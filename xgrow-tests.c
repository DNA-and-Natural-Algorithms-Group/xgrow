/* A test to make sure that xgrow is obeying detailed balance in
 * various situations.  The general idea of these tests is to run some
 * simple assembly reactions until they are approximately at 
 * equilibrium, then test whether the ratio of states at this point
 * look close to detailed balance.
 * 
 * A more quantitative discussion of what we are going to do: 1. We
 * first need to run our simulation until we are "approximately" at
 * equilibrium.  Our approximate test is that designed by Gelman and
 * Rubin (1992) "Inference from Iterative Simulation Using Multiple
 * Sequences" in which we start with an overdispersed set of starting
 * points, and run our simulation until, by a series of T-tests, the
 * variance within the simulation is approximately equal to the
 * variance between simulations.
 * 
 * The test formulated by Gelman and Rubin considers the variance of a
 * variable in order to test for equilibrium.  The indicator variables
 * that we use are k indicator variables that are each 1 if a
 * simulation is in a state s, and 0 otherwise.  The variance
 * within a sequence with respect to state s is the variance of the
 * indicator function I(s), which is 1 if we are in state s at time t,
 * and 0 otherwise.  There is a mean over time of I, and the variance
 * between these means is the other variance that we measure.  We can
 * select these k states by starting the simulation with a random tile
 * and running a simulation for a random period of time (exponentially
 * distributed), and selecting the state that we end up in as one of
 * our k states, throwing out duplicates.
 * 
 * This test requires that the initial state of each chain be drawn
 * from an over-dispersed distribution.  If we have n chains, we can
 * simply doing the same thing as above to generate the n initial
 * states, but by annealing to the desired temperature in order to
 * sample all the possible end states at the desired physical
 * conditions.  We then randomly choose the initial states from all
 * the states visited by a sufficiently slow anneal.  For simplicity,
 * these initial states are also the states used as indicator
 * variables.
 *
 * Once we have determined that we have k chains at stationary state,
 * we evaluate whether the ratio of times spent in each of these k
 * states is compatible with detailed balance.  From this we get a
 * series of k^2 ratios.  "Mathematical Statistics and Data Analysis"
 * by John Rice provides a derivation of computing the expected mean
 * and variance of a ratio of sampled variables.  We calculate this
 * mean and variance, and determine whether we are sufficiently
 * confident that we are obeying detailed balance.
 *
 * Some options (notably fission and chunk_fission) cause a simulation
 * to not obey detailed balance.  It seems possible that at least in
 * the case of the latter, these tests could be used to tune the
 * simulation to be closer to detailed balance.
 *
 * Rebecca Schulman, May 13, 2005
 */

#include <assert.h>
#include <glib-2.0/glib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xgrow-tests.h"

#define UPDATE_RATE 10000
#define SAMPLING_RATE 0.001
#define CHAIN_COUNT 20
#define STATES_TO_ADD_PER_ANNEAL 2
#define BLOCK_TIME 10
#define SCALE_REDUCTION_LIMIT 1.1

#define FLOAT_TOLERANCE 1e-7
#define CONFIDENCE_CONSTANT 1.96
#define MINIMUM_STATES_SEEN 35

double ok_seen_states_ratio = 1.01;
int time_constants_to_run = 10;

double delta = 1;


typedef struct interval_list {
  struct interval_list *next;
  double start_t;
  double end_t;
}  interval_list;

typedef struct chain_state_record {
  double start_t;
  int intervals;
} chain_state_record;

typedef struct indicator_data {
  int assembly;               /* The hash code of the assembly */
  double *means;              /* Means of the indicator variance, for each chain */
  int n;                   /* The number of sampling iterations */
  double mean_of_means;       /* Mean of the means - the target mean */
  double mean_of_square_means; /* Mean of the squares means (needed for covariance calculations) */
  double B;                   /* Variance of the means */
  double *variances;          /* Variances of the indicator variance, for each chain */
  double W;                   /* Mean of the indicator variances  */
  double var_si_squared;      /* Variance of the indicator variances */
  double sigma_hat_squared;    /* Target variance */
  double V_hat;               /* The approximate scale of the t distribution of the means and variances */
  double cov_one;             /* The covariance between the variances and the means squared */
  double cov_two;             /* The covariance between the variances and the means */
  double var_V;               /* The variance of V_hat */
  double df;                  /* The degrees of freedom of the distribution */
  double R_hat;               /* The potential scale reduction if we kept sampling */

} indicator_data;

typedef struct detailed_balance_data {
  double *means;
  double *variances;
  double **covariances;
} detailed_balance_data;

void *malloc_err(size_t size) {
  void *a;
  a = malloc(size);
  if (a == NULL) {
    fprintf(stderr,"Malloc error.  Aborting.\n");
    exit(1);
  }

  return a;
}

int assemblies_are_equal(Assembly a, Assembly b, int size) {
  int j, copy = 1;

  for (j = 0; j < size; j++) {
    if (memcmp(a[j],b[j],size)) {
      copy = 0;
      break;
    }
  }
  return copy;
}

int assembly_is_a_duplicate_in_array (Assembly *assemblies,
				      int number_of_assemblies,
				      Assembly a,
				      int size) {
  int i;

  for (i = 0; i < number_of_assemblies; i++) {
    if (assemblies_are_equal(assemblies[i],a,size)) {
      return 1;
    }
  }
  return 0;
}

const int primes[] = { 503653,503663,503707,503717,503743,503753,
		       503771,503777,503779,503791,503803,
		       24523987,24524009,24524021,24524023,
		       24524033,24524057,24524063,24524089,24524099,24524119,
		       24524173,24524189,24524221,24524243,93851,
		       93871,93887,93889,93893,93901,93911,93913,93923,
		       93937,93941,93949,93967,93971,93979,93983,93997,
		       94007,94009,94033,94049,94057,94063,94079,94099,94109,
		       94111,94117,94121,94151,94153,94169,94201,94207,94219,
		       94229,94253,94261,94273,94291,94307,94309,94321,94327,94331,94343,
		       88747,88771,88789,88793,88799,88801,88807,88811,88813,88817,88819,
		       88843,88853,88861,88867,88873,88883,88897,88903, 97423,97429,97441,
		       97453,97459,97463,97499,97501,97511,97523,97547,97549,97553,97561,
		       97571,97577,97579,97583,97607,97609,97613,97649,97651,97673,97687,
		       97711,97729,97771,97777,97787,97789,97813,97829,97841,97843,97847,
		       97849,97859,97861,97871,97879,97883,97919,97927,97931,97943, 2736637,
		       2736673,2736689,2736691,2736707,2736733,2736737,2736757,2736787,
		       2736863,2736869,2736889,2736911,2736917,2736941,2736947,2736961,2736967,
		       8237461,8237477,8237497,8237513,8237527,8237533,8237543,8237561,
		       8237617,8237627,8237629,8237681,8237683,8237687,8237693,8237699,
		       8237701,8237729,8237753,8237767,8237777,8237783,8237821,8237839,
		       8237869,8237881,8237891,8237899,8237903,8237921,8237923,8237951,
		       8237969,8237981,8238011,8238023,8238029,8238053,8238067,8238101,
		       8238121,8238151,8238169,8238179,823757989,823757993,823758007,
		       823758037,823758077,823758083,823758107,823758109,823758139,823758157,
		       823758163,823758197,823758203,823758211,823758227,823758251,823758263,
		       823758283,823758293,823758337,823758373,823758379,823758421,823758461,
		       823758469,823758473,823758487,823758497,823758503,823758521,823758541,
		       823758557,823758589,823758599,823758641,823758643,78700087,78700091,78700099,
		       78700129,78700133,78700151,78700183,
		       78700199,78700207,78700211,78700217,
		       78700247,78700267,78700301,78700313,
		       78700331,78700373,78700381,78700387,
		       78700397,78700399,78700421,78700423,
		       78700471,78700481,78700483,78700493,
		       78700499,78700529,78700541,78700543,
		       78700571,78700579,78700597,78700613,
		       78700619,78700649,78700681,78700691,
		       78700709,78700711,78700753,78700759,
		       78700763,78700777,78700807,78700819,
		       78700823,78700829,78700841,78700861,
		       78700891,78700903,78700927,78700931,
		       78700933,78700949,78700967,78700969,
		       78700981,78701023,78701039,78701083,
		       78701087,78701089,78701093,78701107,
		       78701111,78701113,78701141,78701149,
		       78701179,78701191,78701197,78701213,
		       78701243,78701251,78701261,78701291,
		       78701317,78701369,78701411,78701419,
		       78701429,78701437,78701443,78701489,
		       78701503,78701527,78701531,78701533,
		       78701549,78701569,78701593,78701611,
		       78701633,78701657,78701663,78701669,
		       78701671,78701681,78701723,78701729,
		       78701731,78701747,78701759,78701761,
		       78701771,78701783,78701807,78701827,
		       78701863,78701869,78701879,78701897,
		       78701921,78701941,78701977,78701981,
		       78701999,78702007,78702053,78702073,
		       78702101,78702131,78702133,78702157,
		       78702181,78702187,78702203,78702227,
		       78702229,78702233,78702263,78702269,
		       78702277,78702301,78702307,78702347,
		       78702391,78702413,78702419,78702443,
		       78702451,78702467,78702479,78702497,
		       78702499,78702511,78702523,78702539,
		       78702553,78702649,78702653,78702661,
		       78702671,78702677,78702683,78702697,
		       78702713,78702751,78702763,78702797,
		       78702817,78702829,78702847,78702853,
		       78702859,78702917,78702919,78702929,
		       78702931,78702947,78702961,78702971,
		       78702977,78703061,78703087,78703103,
		       78703127,78703133,78703147,78703151,
		       78703153,78703159,78703241,78703243,
		       78703253,78703267,78703283,78703321,
		       78703337,78703351,78703363,78703367,
		       78703393,78703399,78703487,78703507,
		       78703561,78703579,78703589,78703657,
		       78703673,78703721,78703753,78703759,
		       78703783,78703787,78703789,78703799,
		       78703813,78703831,78703841,78703873,
		       78703901,78703913,78703939,78703949,
		       78703967,78704027,78704039,78704063,
		       78704089,78704107,78704113,78704123,
		       78704159,78704167,78704191,78704207,
		       78704221,78704267,78704317,78704357,
		       78704363,78704377,78704389,78704393,
		       78704399,78704401,78704419,78704429};

/* Don't use the edges of the assembly for hashing */
int hash_assembly (Assembly a, int size) {
  int i,*j,k;
  int h, hash = 0;
  for (i = 1; i < (size + 1); i++) {
    h = 0;
    for (k = 1; k < (size + 1); k+=4) {
      j = (int *) (a[i] + k);
      h += (*j << ((k - i) % 32)) + (*j >> (32 - (k-i)%32));
    }
    hash += h % primes[i];
  }
  return hash;
}


int assembly_is_a_duplicate (void *states_seen,
			     Assembly a, int size) {
  int key;
  int dup;

  key = hash_assembly (a, size);
  dup = (g_hash_table_lookup ((GHashTable *) states_seen, &key) != NULL);
  //printf ("Is duplicate : %d.\n",dup);
  return dup;
}

void reset_tube (tube *tp) {
  double max_bond_strength;
  int i,j;

  tp->states_seen_count = 0;
  tp->states_seen_hash = g_hash_table_new (g_int_hash, g_int_equal);
  tp->t = 0;
  tp->events = 0;
  tp->stat_a = 0;
  tp->stat_d = 0;
  tp->num_flakes = 0;
  // TODO:  If Gfc is nonzero, restore original concentrations
  tp->flake_list=NULL;
  tp->flake_tree=NULL;

  // Set annealing start point so that tau = max_bond_strength + 1
  max_bond_strength = 0;
  for (i = 0; i < tp->num_bindings; i++) {
    max_bond_strength = MAX(tp->strength[i],max_bond_strength);
    for (j = 0; j < tp->num_bindings; j++) {
      max_bond_strength = MAX(tp->glue[i][j],max_bond_strength);
    }
  }
  tp->anneal_g = tp->Gmc / (max_bond_strength + 1);
}


void print_assembly (Assembly cell, int size) {
  int row, column = 0; 
  int i,j;

  for (i = 1; i < size + 1; i++) {
    row = 0;
    column = 0;
    for (j = 1; j < size + 1; j++) {
      if (cell[i][j]) {
	row = i;
	if (!column) {
	  column = j;
	}
      }
    }
    if (row) {
      printf("Row=%d Column=%d ",row, column);
      for (j = 1; j < size; j++) {
	if (cell[i][j] != 0) {
	  printf("%d ",cell[i][j]);
	}
	else {
	  printf("  ");
	}
      }
      if (cell[i][size]) {
	printf("%d\n",cell[i][size]);
      }
      else {
	printf("\n");
      }
    }
  }
  printf("\n");
}

Assembly copy_assembly (Assembly cell, int size) {
  Assembly new;
  int i;
  
  new = (Assembly) malloc_err ((size+2) * sizeof(unsigned char *));
  for (i = 0; i < size + 2; i++) {
    new[i] = (unsigned char *) malloc_err((size+2)*sizeof(unsigned char));
    memcpy (new[i],&cell[i][0], (size+2)*sizeof(unsigned char));
  }
  return new;
}


void add_assembly_to_seen (tube *tp) {
  Assembly cur;
  int size, *hash;


  cur = tp->flake_list->cell;
  size = (1<<(tp->flake_list->P));
  //printf("Adding.\n");
  //print_assembly (cur,size);  
  hash = malloc(sizeof(int));
  *hash = hash_assembly (cur, size);
  g_hash_table_insert ((GHashTable *) tp->states_seen_hash, hash, copy_assembly (cur, size));
  tp->states_seen_count++;
}


void free_assembly (Assembly x, int size) {
  int i;
  
  for (i = 0; i < size + 2; i++) {
    free(x[i]);
  }
  free (x);
}



/* Sometimes, change cell "adds" an assembly that it then removes,
   because it is disconnected.  We should not consider a seen state
   either */
void remove_assembly_from_seen (tube *tp) {
  Assembly cur;
  int size, *hash;
  Assembly a;


  cur = tp->flake_list->cell;
  size = (1<<(tp->flake_list->P));
  //printf("Removing.\n");
  //print_assembly (cur,size);  
  hash = malloc(sizeof(int));
  *hash = hash_assembly (cur, size);
  a = g_hash_table_lookup ((GHashTable *) tp->states_seen_hash, hash);
  free_assembly (a, size);
  g_hash_table_remove ((GHashTable *) tp->states_seen_hash, hash);
  tp->states_seen_count--;
} 
  


void free_key (gpointer key,
	       gpointer value,  
	       gpointer user_data) {
  free (key);
}

void free_assembly_value (gpointer key,
			  gpointer value,  
			  gpointer user_data) {
  int size;
  size = (int) user_data;
  free_assembly (value, size);
}

void clear_seen_states (tube *tp) {
  int size;

  size = 1<<(tp->P);
  g_hash_table_foreach(tp->states_seen_hash, free_key,NULL);
  g_hash_table_foreach(tp->states_seen_hash, free_assembly_value,(gpointer) size);
  g_hash_table_destroy(tp->states_seen_hash);
		      
}

typedef struct maybe_add_data {
  int total_states_added;
  int states_added_this_anneal;
  tube *tp;
} maybe_add_data;

void maybe_add_to_chain_states (gpointer key,
				gpointer value,
				gpointer user_data) {
  int *assembly_hash;
  int *assembly_code;
  maybe_add_data *data;
  tube *tp;
  int size;
  Assembly a;

  assembly_hash = (int *) key;
  data = (maybe_add_data *) user_data;
  tp = data->tp;
  size = 1<<(tp->P);
  
  if (data->total_states_added < tp->chains &&
      data->states_added_this_anneal < STATES_TO_ADD_PER_ANNEAL &&
      random () > ((tp->states_seen_count - tp->chains) / 
		   tp->states_seen_count) * RAND_MAX &&
      !g_hash_table_lookup (tp->chain_states, assembly_hash)) {
    
    assembly_code = (int *) malloc (sizeof (int));
    assembly_code = assembly_hash;
    g_hash_table_insert (tp->chain_states, assembly_code, assembly_code);
    a = g_hash_table_lookup (tp->states_seen_hash, assembly_code);
    assert (a);
    tp->start_states[data->total_states_added] = 
      copy_assembly(a,size);

    printf("Assembly %d:\n",data->total_states_added);
    print_assembly (tp->start_states[data->total_states_added],size);
    //printf("Assembly hash code is %d.\n",*assembly_code);
    data->total_states_added++;
    data->states_added_this_anneal++;
  }
  
}

void generate_initial_chain_states (tube *tp, int seed_i, int seed_j, int seed_n) {
  unsigned char tile_to_start_with;
  double time_to_run;
  int last_seen_states = 0;
  flake *fp;
  int size;
  maybe_add_data data_s;
  int total_states_added, states_added_this_anneal;
  
  tp->chain_states = g_hash_table_new (g_int_hash, g_int_equal);
  tp->start_states = (Assembly *) malloc_err(tp->chains * sizeof(Assembly));
  tp->start_state_Gs = (double *) malloc_err(tp->chains * sizeof(double));

  printf("Finding indicator variables.\n");
  last_seen_states = 0;
  time_constants_to_run = 1;
  tp->anneal_t = 0.001;
  total_states_added = 0;
  while (total_states_added < tp->chains) {
    while (1) {
      //tile_to_start_with = (unsigned char) ((random () / ((double) RAND_MAX)) * ((double) tp->N)) + 1;
      tile_to_start_with = seed_n;
      //printf("Starting with tile %d.\n",tile_to_start_with);
      reset_tube (tp);
      fp=init_flake(tp->P,tp->N,
		    seed_i,seed_j,tile_to_start_with,0);
      if (tp->dt_right[tile_to_start_with] || tp->dt_left[tile_to_start_with])
	fp->seed_is_double_tile = 1;
      if (tp->dt_right[tile_to_start_with])
	change_cell(fp,seed_i,seed_j+1,tp->dt_right[tile_to_start_with]);
      if (tp->dt_left[tile_to_start_with])
	change_cell(fp,seed_i,seed_j-1,tp->dt_left[tile_to_start_with]);
      insert_flake(fp, tp);    
      size = 1<<(fp->P);
      time_to_run = tp->anneal_t*log(2.0)*time_constants_to_run;
      printf("Time to run is %e.\n",time_to_run);
      time_to_run = tp->anneal_t*time_constants_to_run;
      tp->tracking_seen_states = 1;
      while (tp->t < time_to_run) {
	simulate (tp, UPDATE_RATE, time_to_run, 0, 0, 0);
      }
      if (((double) tp->states_seen_count / ok_seen_states_ratio) <= last_seen_states && 
	  tp->states_seen_count - STATES_TO_ADD_PER_ANNEAL > total_states_added  &&
	  tp->states_seen_count > MINIMUM_STATES_SEEN) {
	printf("Saw %d states in the final round.\n",tp->states_seen_count);
	break;
      }
      else {
	last_seen_states = tp->states_seen_count;
	clear_seen_states (tp);
	free_flake(fp);
	tp->flake_list = NULL;
	tp->anneal_t *= 1.5;
	printf("Increasing time constant to %lf.\n",tp->anneal_t);
	printf("States seen this round were %d.\n",last_seen_states);
      }
    }
    states_added_this_anneal = 0;
    while (total_states_added < tp->chains &&
	   states_added_this_anneal < STATES_TO_ADD_PER_ANNEAL) {
      data_s.states_added_this_anneal = states_added_this_anneal;
      data_s.total_states_added = total_states_added;
      data_s.tp = tp;
      g_hash_table_foreach (tp->states_seen_hash, maybe_add_to_chain_states, &data_s);
      total_states_added = data_s.total_states_added;
      states_added_this_anneal = data_s.states_added_this_anneal;
    }
  }
  tp->tracking_seen_states = 0;
}

static void free_start_states (Assembly *a, int size, int assemblies) {
  int i, j;
  for (i = 0; i < assemblies; i++) {
    for (j = 0; j < size + 2; j++) {
      free(a[i][j]);
    }
    free(a[i]);
  }
  free(a);
}

void print_key (gpointer key, gpointer value, gpointer user_data) {
  int *key_int;
  key_int = (int *) key;
  printf("Key is %d.\n",*key_int);
}

int sample_count (double start_time, double end_time) {
  int contains_an_interval_sample;
  int extra_samples;

  contains_an_interval_sample = 
    floor (start_time / SAMPLING_RATE) != floor (end_time / SAMPLING_RATE);
  extra_samples = MAX(0,floor((end_time - start_time)/SAMPLING_RATE));
  return contains_an_interval_sample + extra_samples;
}

double variance_total (flake *flake, double mean, int assembly_code, double cur_time) {
  double v;
  chain_state_record *c;

  c = g_hash_table_lookup (flake->chain_hash, &assembly_code); 
  if (c) {
    v = c->intervals * pow(1 - mean,2) + 
      (sample_count (0, cur_time) - c->intervals) * pow(mean,2);
  }
  else {
    v = 0;
  }
  return v;
}

//#define DEBUG_CONVERGED 1

int converged (tube *tp, indicator_data *data) {
  int size, i, j, count_total;
  flake *flake;
  double n,m;
  chain_state_record *c;

  for (j = 0; j < tp->chains; j++) {
#ifdef DEBUG_CONVERGED
    printf("For indicator variable %d:\n",j);
#endif
    size = 1<<(tp->P);
    data[j].assembly = hash_assembly(tp->start_states[j],size);
    data[j].n = sample_count (0, tp->t);

    n = (double) data[j].n;
    m = (double) tp->chains;
    
    /* Calculate mean of the indicator variable for each chain */
    i = 0;
    data[j].mean_of_means = 0;
    data[j].mean_of_square_means = 0;
    for (flake = tp->flake_list; flake != NULL; flake = flake->next_flake) {
      count_total = 0;
      c = g_hash_table_lookup (flake->chain_hash, &(data[j].assembly)); 
      if (c) {
	data[j].means[i] = c->intervals / n;
      }
      else {
	data[j].means[i] = 0;
      }
#ifdef DEBUG_CONVERGED
      printf("Mean %d is %e.\n",i,data[j].means[i]);
#endif
      data[j].mean_of_means += data[j].means[i];
      data[j].mean_of_square_means += pow(data[j].means[i],2);
      i++;
    }
    data[j].mean_of_means /= tp->chains;    
    data[j].mean_of_square_means /= tp->chains;
#ifdef DEBUG_CONVERGED
    printf("Time now is %e.  Mean percentage time for all flakes is %e.\n",tp->t,data[j].mean_of_means);
#endif

    /* Calculate variance of the indicator variables means for each chain */
    i = 0;
    for (flake = tp->flake_list; flake != NULL; flake = flake->next_flake) {
      data[j].variances[i] = variance_total (flake, data[j].means[i], data[j].assembly, tp->t);
      data[j].variances[i] /= n - 1;
#ifdef DEBUG_CONVERGED
      printf("Variance of chain %d is %e.\n",i,data[j].variances[i]);
#endif
      i++;
    }
    
    data[j].B = 0;
    for (i = 0; i < m; i++) {
#ifdef DEBUG_CONVERGED
      printf("Variance term for chain %d is %e.\n",i,pow (data[j].mean_of_means - data[j].means[i],2));
#endif
      data[j].B += pow (data[j].mean_of_means - data[j].means[i],2);
    }
    data[j].B /= (m - 1);
    data[j].B *= n;
#ifdef DEBUG_CONVERGED
    printf("Variance of the means (B) is %e.\n",data[j].B);
#endif

    /* Calculate variance of the indicator variable for each chain */
    data[j].W = 0;
    for (i = 0; i < m; i++) {
      data[j].W += data[j].variances[i];
    }
    /* Calculate the mean of the variances */
    data[j].W /= m;
#ifdef DEBUG_CONVERGED
    printf("Mean variance (W) is %e.\n",data[j].W);
#endif
  
    /* Calculate the target variance */
    data[j].sigma_hat_squared = (n - 1)/n * data[j].W + data[j].B/n;
#ifdef DEBUG_CONVERGED
    printf("Sigma hat squared is %e.\n", data[j].sigma_hat_squared);
#endif
    /* Now we are going to calculate a student's T distribution to
       calculate the likely mean and variance of the target. To do this
       we calculate the scale of the distribution */
    data[j].V_hat = data[j].sigma_hat_squared + data[j].B/(m*n);
#ifdef DEBUG_CONVERGED
    printf("Target scaling (V_hat) is %e.\n",data[j].V_hat);
#endif
    data[j].var_si_squared = 0;
    data[j].cov_one = 0;
    data[j].cov_two = 0;
    for (i = 0; i < m; i++) {
      data[j].var_si_squared += (data[j].W - data[j].variances[i]) * (data[j].W - data[j].variances[i]);
      data[j].cov_one += (data[j].mean_of_square_means - data[j].means[i]*data[j].means[i]) * (data[j].W - data[j].variances[i]);
      data[j].cov_two += (data[j].mean_of_means - data[j].means[i]) * (data[j].W - data[j].variances[i]);
    }
    data[j].var_si_squared /= (m - 1);
#ifdef DEBUG_CONVERGED
    printf("Variance of si squared is %e.\n",data[j].var_si_squared);
#endif
    data[j].cov_one /= m;
    data[j].cov_two /= m;
#ifdef DEBUG_CONVERGED
    printf("Covariance one is %e, covariance two is %e.\n",data[j].cov_one, data[j].cov_two);
#endif
    data[j].var_V = ((n - 1) / n) * ((n - 1) / n) / m * data[j].var_si_squared + 
      ((m + 1) / (m*n)) * ((m + 1) / (m*n)) * 2 / (m - 1) * data[j].B * data[j].B + 
      2 * ((m + 1) * (n - 1)) / (m * m * n) *
      (data[j].cov_one - 2*data[j].mean_of_means*data[j].cov_two);

#ifdef DEBUG_CONVERGED
    printf("The variance on V is %e.\n",data[j].var_V);
#endif
    data[j].df = 2*pow (data[j].V_hat,2) / data[j].var_V;
#ifdef DEBUG_CONVERGED
    printf("Degrees of freedom in variance are %e.\n",data[j].df);
#endif
    data[j].R_hat = data[j].V_hat/data[j].W*data[j].df/(data[j].df - 2);
#ifdef DEBUG_CONVERGED
    data[j].R_hat = data[j].V_hat/data[j].W;
#endif
    printf("Potential scale reduction for indicator variable %d is %e.\n",j,data[j].R_hat);
    if (isnan (data[j].R_hat) || data[j].R_hat > SCALE_REDUCTION_LIMIT) {
      printf("Scale factor is too high, continuing.\n");
      return 0;
    }
    if (data[j].R_hat < 0) {
      printf("Not enough data to compute scale factor. Continuing.\n");
      return 0;
    }
}
  //return return_value;
  return 1;
}



static 
indicator_data *run_flakes_past_burn(tube *tp, int size) {
  int i,j,k;
  flake *fp;
  int l,m;
  int not_empty;
  indicator_data *data;

  reset_tube(tp);
  tp->anneal_t = 0;
  tp->Gse = tp->Gse_final;
  data = (indicator_data *) malloc (tp->chains*sizeof(indicator_data));
  for (j = 0; j < tp->chains; j++) {
    data[j].means = (double *) malloc(tp->chains*sizeof(double));
    data[j].variances = (double *) malloc(tp->chains*sizeof(double));
  }

  for (i = 0; i < tp->chains; i++) {
    //printf("Setting up flake %d.\n",i);
    fp=init_flake(tp->P,tp->N,1,1,1,0);
    insert_flake(fp, tp);    
    not_empty = 0;
    for (j = 1; j < size + 1; j++) {
      for (k = 1; k < size + 1; k++) {
	assert (tp->start_states[i]);
	change_cell (fp, j-1,k-1, tp->start_states[i][j][k]);
	if (tp->start_states[i][j][k] > 0) {
	  not_empty = 1;
	}
      }
    }
    assert (not_empty);
    l = size * (((double)random()) / ((double)RAND_MAX));
    m = size * (((double)random()) / ((double)RAND_MAX));
    while ((fp->Cell(l,m)) == 0) {
      l = size * (((double)random()) / ((double)RAND_MAX));
      m = size * (((double)random()) / ((double)RAND_MAX));
    }
    fp->seed_i = l;
    fp->seed_j = m;
    fp->seed_n = fp->Cell(l,m);
    recalc_G(fp);
    tp->start_state_Gs[i] = fp->G;
    fp->chain_hash = g_hash_table_new_full (g_int_hash, g_int_equal, 
					    free, free);
    update_state_on_indicator(fp,tp->start_states[i], size);
  }
  tp->watching_states = 1;
  printf("Starting simulation to pass burn.\n");
  i = 0;
  while (1) {
    printf("Simulating block %d.\n",i);
    while (tp->t < BLOCK_TIME*i) {
      simulate (tp, UPDATE_RATE, BLOCK_TIME*i, 0, 0, 0);
      //printf("Time is %e.\n",tp->t);
    }
    /* After running a block, recalculate parameters to see if we are past burn */
    if (converged (tp, data)) {
      break;
    }
    i++;
  }
  printf("Time at end is %e.\n",tp->t);
  printf("Freeing states.\n");
  free_start_states (tp->start_states, size, tp->chains);
  return data;
}

void update_state_on_indicator(flake *fp, Assembly a, int size) {
  int h, *k;
  chain_state_record *j;
  tube *tp;

  tp = fp->tube;
  h = hash_assembly (a, size);
  if (g_hash_table_lookup(tp->chain_states, &h)) {
    if ((j = g_hash_table_lookup (fp->chain_hash, &h))) {
      j->start_t = tp->t;
    }
    else {
      j = (chain_state_record *) malloc(sizeof (chain_state_record));
      j->start_t = tp->t;
      j->intervals = 0;
      k = (int *) malloc(sizeof(int));
      *k = h;
      g_hash_table_insert(fp->chain_hash, k, j);
    }
    fp->chain_state = h;
    //printf("Entering chain state %d for flake %p.\n",fp->chain_state,fp);
  }
  else {
    //printf("Rejecting state %d:\n",h);
    //print_assembly (a, size);
  }
}

void update_state_off_indicator(flake *fp) {
  chain_state_record *i;
  //printf("Turning off chain state %d for flake %p.\n",fp->chain_state,fp);
  i = g_hash_table_lookup(fp->chain_hash, &fp->chain_state);
  assert (i);
  i->intervals += sample_count (i->start_t, fp->tube->t);
  fp->chain_state = 0;
}

int test_detailed_balance (tube *tp, indicator_data *data) {
  int i,j;
  flake *flake;
  detailed_balance_data d;
  double correct_state_ratio;
  double r, est_mean, est_variance, confidence_interval;
  double confidence_ratio;
  int n;

  printf("Testing whether detailed balance has been achieved.\n");
  d.means = (double *) malloc (sizeof(double) * tp->chains);
  d.variances = (double *) malloc (sizeof(double) * tp->chains);
  d.covariances = (double **) malloc (sizeof(double *) * tp->chains);
  for (i = 0; i < tp->chains; i++) {
    d.covariances[i] = (double *) malloc (sizeof(double) * tp->chains);
  }

  for (i = 0; i < tp->chains; i++) {
    d.means[i] = data[i].mean_of_means;
  }
 
  n = (floor(tp->t / SAMPLING_RATE)  - 1) * tp->chains;  
  for (i = 0; i < tp->chains; i++) {
    d.variances[i] = 0;
    for (flake = tp->flake_list; flake != NULL; flake = flake->next_flake) {
      d.variances[i] += variance_total (flake, d.means[i], data[i].assembly, tp->t);
    }
    d.variances[i] /= n - 1;
  }

  /* Calculate covariance of indicator variables i and j */
  for (i = 0; i < tp->chains; i++) {
    for (j = 0; j < tp->chains; j++) {
      d.covariances[i][j] = 0;
      d.covariances[i][j] += (1 - d.means[i]) * d.means[j] * d.means[i];
      d.covariances[i][j] += (1 - d.means[j]) * d.means[i] * d.means[j];
      d.covariances[i][j] += d.means[i] * d.means[j] * (1 - d.means[i] - d.means[j]);
    }
  }
      
  for (i = 0; i < tp->chains; i++) {
    for (j = 0; j < tp->chains; j++) {
      if (j == i) {
	continue;
      }
      correct_state_ratio = exp(tp->start_state_Gs[j]-tp->start_state_Gs[i]);
      r = d.means[i] / d.means[j];
      est_variance = 1/(n * pow(d.means[j],2)) * 
	(r*r*d.variances[j] + d.variances[i] - 2*r*d.covariances[i][j]);
      //printf("Estimated variance is %e.\n",est_variance);
      est_mean = r + 1/(n * pow(d.means[j],2)) * 
	(r*pow(d.variances[j],2) - 2*r*d.covariances[i][j]);
      if (est_variance >= 0) {
	confidence_interval = CONFIDENCE_CONSTANT * sqrt(est_variance);
      }
      else {
	confidence_interval = CONFIDENCE_CONSTANT * sqrt(-est_variance);
      }
      confidence_ratio = confidence_interval/correct_state_ratio;
      printf("\nTrue ratio of variables %d and %d is %f.\n",j,i,correct_state_ratio);
      printf("Simulated ratio is %f.  Confidence interval is %f, or %2.1f%%.\n",
	     est_mean, confidence_interval,100*confidence_ratio);
      if (est_mean - 3*confidence_interval > correct_state_ratio) {
	printf("Estimated ratio is too high.\n");
	//return 0;
      }
      if (est_mean + 3*confidence_interval < correct_state_ratio) {
	printf("Estimated ratio is too low.\n");
	//return 0;
      }
    }
  }
  return 1;
}

void run_xgrow_tests (tube *tp,double Gmc, double Gse, int seed_i, int seed_j, int seed_n, int size) {
  indicator_data *data;

  tp->chains = CHAIN_COUNT;
  /* First, generate the initial set of random states, by starting
   * with a random tile and running until we've seen most of the
   * states we will see locally */
  generate_initial_chain_states(tp, seed_i, seed_j, seed_n);

  /* Now, we run the chains with the initial states until we pass our
     t-test on the variance of the random states we're sampling occurring */
  data = run_flakes_past_burn(tp, size);

  /* Now that we're past the initial "burn" in which our state is
     highly dependent on the starting state, we can test whether we
     obey detailed balance, by comparing the time spent in different
     states by our different processes. */
  if (test_detailed_balance (tp, data)) {
    printf("Detailed balance seems to have been achieved.\n");
  }
  else {
    printf("There was a problem.  This may be a statistical error, but if it repeats, there may be a bug.\n");
  }
  
}
