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
#include <openssl/md4.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xgrow-tests.h"

#define CHOP 10

#define UPDATE_RATE 10000
#define CHAIN_COUNT 20
#define STATES_TO_ADD_PER_ANNEAL 2
#define SCALE_REDUCTION_LIMIT 1.01

#define FLOAT_TOLERANCE 1e-7
#define CONFIDENCE_CONSTANT 3.30
#define MINIMUM_STATES_SEEN 20

double ok_seen_states_ratio = 1.01;
int time_constants_to_run = 10;

double block_time, sampling_rate;
int seed_i, seed_j, seed_n;

typedef struct interval_list {
  struct interval_list *next;
  double start_t;
  double end_t;
}  interval_list;

typedef struct chain_state_record {
  double start_t;
  double old_start_t;
  int intervals;
  double times;
} chain_state_record;

typedef struct indicator_data {
  unsigned char *assembly;    /* The hash code of the assembly */
  double *means;              /* Means of the indicator variance, for each chain */
  int n;                   /* The number of sampling iterations */
  double mean_of_means;       /* Mean of the means - the target mean */
  double total_time;         /* The total amount of time spent by all chains in this state */
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
unsigned char *hash_assembly(Assembly a, int size) {
  int blocks=0, in_block=0, non_zero_cells=0;
  int i,j, pos, len;
  unsigned char *result;

  for (i = 0; i < (size+1); i++) {
    for (j = 0; j < (size+1); j++) {
      if (a[i][j]) {
	non_zero_cells++;
	if (!in_block) {
	  in_block = 1;
	  blocks++;
	}
      }
      else if (in_block) {
	in_block = 0;
      }
    }
  }
  len = non_zero_cells + 3*blocks;
  result = (unsigned char *) malloc_err ((len + 1)*sizeof (unsigned char)); 
  pos = 0;
  in_block = 0;
  for (i = 0; i < (size+1); i++) {
    for (j = 0; j < (size+1); j++) {
      assert (pos <= len);
      if (a[i][j]) {
	if (in_block) {
	  result[pos++] = a[i][j];
	}
	else {
	  result[pos++] = 255;
	  result[pos++] = i + 1;
	  result[pos++] = j + 1;
	  result[pos++] = a[i][j];
	  in_block = 1;
	}
      }
      else {
	in_block = 0;
      }
    }
  }
  result[pos] = 0;
  return result;
}

int assembly_is_a_duplicate (void *states_seen,
			     Assembly a, int size) {
   int dup;
  unsigned char *key;
  
  key = hash_assembly (a, size);
  dup = (g_hash_table_lookup ((GHashTable *) states_seen, key) != NULL);
  free (key);
  //printf ("Is duplicate : %d.\n",dup);
  return dup;
}

void reset_tube (tube *tp) {
  double max_bond_strength;
  int i,j;

  tp->states_seen_count = 0;
  tp->states_seen_hash = g_hash_table_new (g_str_hash, g_str_equal);
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
  int size;
  unsigned char *hash;

  cur = tp->flake_list->cell;
  size = (1<<(tp->flake_list->P));
  //printf("Adding.\n");
  //print_assembly (cur,size);  
  hash = hash_assembly (cur, size);
  //printf("hash %s.\n",hash);
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
  int size;
  Assembly a;
  unsigned char *hash;

  cur = tp->flake_list->cell;
  size = (1<<(tp->flake_list->P));
  //printf("Removing.\n");
  //print_assembly (cur,size);  
  hash = hash_assembly (cur, size);
  a = g_hash_table_lookup ((GHashTable *) tp->states_seen_hash, hash);
  free_assembly (a, size);
  g_hash_table_remove ((GHashTable *) tp->states_seen_hash, hash);
  free (hash);
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
  unsigned char *assembly_hash;
  unsigned char *assembly_code;
  maybe_add_data *data;
  tube *tp;
  int size;
  Assembly a;

  assembly_hash = (unsigned char *) key;
  data = (maybe_add_data *) user_data;
  tp = data->tp;
  size = 1<<(tp->P);
  if (data->total_states_added < tp->chains &&
      data->states_added_this_anneal < STATES_TO_ADD_PER_ANNEAL &&
      drand48() < ((double) STATES_TO_ADD_PER_ANNEAL  /
		   (double) tp->states_seen_count) &&
      !g_hash_table_lookup (tp->chain_states, assembly_hash)) {
    assembly_code = (unsigned char *) malloc_err((strlen((char *) assembly_hash)+1)*sizeof (unsigned char));
    memcpy(assembly_code, assembly_hash, strlen((char *) assembly_hash) + 1);
    g_hash_table_insert (tp->chain_states, assembly_code, assembly_code);
    a = g_hash_table_lookup (tp->states_seen_hash, assembly_code);
    assert (a);
    tp->start_states[data->total_states_added] = 
      copy_assembly(a,size);
    printf("Choosing representative assembly %d:\n",data->total_states_added);
    print_assembly (tp->start_states[data->total_states_added],size);
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
  
  tp->chain_states = g_hash_table_new (g_str_hash, g_str_equal);
  tp->start_states = (Assembly *) malloc_err(tp->chains * sizeof(Assembly));
  tp->start_state_Gs = (double *) malloc_err(tp->chains * sizeof(double));

  printf("\n**************************************************\n");
  printf("Entering stage one : Finding representative states\n\n");
  printf("**************************************************\n\n");
  printf("We'll start with a fast anneal, and slow down until we've sampled \n"
	 "enough states to test.\n\n");
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
      printf("Current annealing time constant is %1.1e.\n",time_to_run);
      time_to_run = tp->anneal_t*time_constants_to_run;
      tp->tracking_seen_states = 1;
      add_assembly_to_seen(tp);
      while (tp->t < time_to_run) {
	simulate (tp, UPDATE_RATE, time_to_run, 0, 0, 0, -1);
      }
      if (((double) tp->states_seen_count / ok_seen_states_ratio) <= last_seen_states && 
	  tp->states_seen_count - STATES_TO_ADD_PER_ANNEAL > total_states_added  &&
	  tp->states_seen_count > MINIMUM_STATES_SEEN) {
	printf("\nSaw %d states in this anneal.  \n"
	       "Sampling representative states from this anneal:\n",
	       tp->states_seen_count);
	break;
      }
      else {
	last_seen_states = tp->states_seen_count;
	clear_seen_states (tp);
	free_flake(fp);
	tp->flake_list = NULL;
	tp->anneal_t *= 1.5;
      }
    }
    states_added_this_anneal = 0;
    printf("\n");
    while (total_states_added < tp->chains &&
	   states_added_this_anneal < STATES_TO_ADD_PER_ANNEAL) {
      data_s.states_added_this_anneal = states_added_this_anneal;
      data_s.total_states_added = total_states_added;
      data_s.tp = tp;
      g_hash_table_foreach (tp->states_seen_hash, maybe_add_to_chain_states, &data_s);
      total_states_added = data_s.total_states_added;
      states_added_this_anneal = data_s.states_added_this_anneal;
    }
    printf("--------------------\n\n");
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
/*
int sample_count (double start_time, double end_time) {
  int contains_an_interval_sample;
  int extra_samples;

  contains_an_interval_sample = 
    floor (start_time / sampling_rate) != floor (end_time / sampling_rate);
  extra_samples = MAX(0,floor((end_time - start_time)/sampling_rate));
  return contains_an_interval_sample + extra_samples;
}
*/

int sample_count (double start_time, double end_time) {
  int contains_an_interval_sample;
  int extra_samples;

  contains_an_interval_sample = 
    (((int) floor (start_time / sampling_rate))/CHOP) != 
    (((int) floor (end_time / sampling_rate))/CHOP);
  /*if (contains_an_interval_sample) {
  printf("without chop %d to %d.  With chop %d to %d\n",
	   (int) floor (start_time / sampling_rate), (int) floor (end_time / sampling_rate),
	   ((int) floor (start_time / sampling_rate))/CHOP, 
	   ((int) floor (end_time / sampling_rate))/CHOP);
	   }*/
  extra_samples = MAX(0,floor((end_time - start_time)/sampling_rate)/CHOP);
  return contains_an_interval_sample + extra_samples;
}


double variance_total (flake *flake, double mean, unsigned char *assembly_code, double cur_time) {
  double v;
  chain_state_record *c;

  c = g_hash_table_lookup (flake->chain_hash, assembly_code); 
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
    data[j].total_time = 0;
    data[j].mean_of_square_means = 0;
    for (flake = tp->flake_list; flake != NULL; flake = flake->next_flake) {
      count_total = 0;
      c = g_hash_table_lookup (flake->chain_hash, (data[j].assembly)); 
      if (c) {
	data[j].means[i] = c->intervals / n;
	data[j].total_time += c->times;
      }
      else {
	data[j].means[i] = 0;
	data[j].total_time = 0;
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
    if (isnan (data[j].R_hat) || data[j].R_hat > SCALE_REDUCTION_LIMIT) {
      printf("Potential scale reduction for representative state %d is %1.4f.\n",j,data[j].R_hat);
      printf("Scale factor is too high, continuing.\n");
      return 0;
    }
    if (data[j].R_hat < 0) {
      printf("Not enough data to compute scale factor for representative state %d. \n"
	     "Continuing.\n",j);
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
  set_Gses(tp,tp->Gse,0);
  data = (indicator_data *) malloc_err (tp->chains*sizeof(indicator_data));
  for (j = 0; j < tp->chains; j++) {
    data[j].means = (double *) malloc_err(tp->chains*sizeof(double));
    data[j].variances = (double *) malloc_err(tp->chains*sizeof(double));
  }

  for (i = 0; i < tp->chains; i++) {
    //printf("State %d followed:\n",i);
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
    if (!wander) {
      assert (fp->Cell(seed_i,seed_j) == seed_n);
    }

    if (wander) {
      l = size * (((double)random()) / ((double)RAND_MAX));
      m = size * (((double)random()) / ((double)RAND_MAX));
      while (fp->Cell(l,m) == 0 || tp->dt_left[fp->Cell(l,m)]) {
	l = size * (((double)random()) / ((double)RAND_MAX));
	m = size * (((double)random()) / ((double)RAND_MAX));
      }
      fp->seed_i = l;
      fp->seed_j = m;
      fp->seed_n = fp->Cell(l,m);
      fp->seed_is_double_tile = tp->dt_right[fp->Cell(l,m)];
    }
    else {
      fp->seed_i = seed_i;
      fp->seed_j = seed_j;
      fp->seed_n = seed_n;
      fp->seed_is_double_tile = tp->dt_right[seed_n];
    }
    //print_assembly(fp->cell,1<<(tp->P));
    //printf("seed is %d.\n",fp->seed_n);
    recalc_G(fp);
    tp->start_state_Gs[i] = fp->G;
    fp->chain_hash = g_hash_table_new_full (g_str_hash, g_str_equal, 
					    free, free);
    fp->chain_state = NULL;
    update_state_on_indicator(fp,tp->start_states[i], size);
  }
  tp->watching_states = 1;
  printf("\n**************************************************************\n");
  printf("Entering stage two: Running simulation until chains have seen \n"
	 "representative states about equally.\n");
  printf("**************************************************************\n");

  i = 0;
  while (1) {
    printf("\nTotal simulated time is %f seconds.\n",((double) i)*((double) block_time));
    printf("Simulating block %d:\n",i);
    while (tp->t < block_time*i) {
      simulate (tp, UPDATE_RATE, block_time*i, 0, 0, 0, -1);
      //printf("Time is %e.\n",tp->t);
    }
    /* After running a block, recalculate parameters to see if we are past burn */
    if (converged (tp, data)) {
      break;
    }
    i++;
  }
  free_start_states (tp->start_states, size, tp->chains);
  return data;
}

void update_state_on_indicator(flake *fp, Assembly a, int size) {
  chain_state_record *j;
  tube *tp;
  unsigned char *h, *c, *old_key;

  tp = fp->tube;
  h = hash_assembly (a, size);
  assert (fp->chain_state == NULL);
  if (g_hash_table_lookup(tp->chain_states, h)) {
    if (g_hash_table_lookup_extended (fp->chain_hash, h, 
				      (gpointer *) &old_key,
				      (gpointer *) &j)) {
      j->old_start_t = j->start_t;
      j->start_t = tp->t;
      fp->chain_state = old_key;
    }
    else {
      j = (chain_state_record *) malloc_err(sizeof (chain_state_record));
      j->old_start_t = 0;
      j->start_t = tp->t;
      j->intervals = 0;
      j->times = 0;
      c = (unsigned char *) malloc_err((strlen((char *) h) + 1)*sizeof (unsigned char));
      memcpy(c,h,strlen((char *)h) + 1);
      g_hash_table_insert(fp->chain_hash, c, j);
      fp->chain_state = c;
    }

    //printf("Entering chain state %d for flake %p.\n",fp->chain_state,fp);
  }
  else {
    //printf("Rejecting state %d:\n",h);
    //print_assembly (a, size);
  }
  free (h);
}

void update_state_off_indicator(flake *fp) {
  chain_state_record *i;
  //printf("Turning off chain state %d for flake %p.\n",fp->chain_state,fp);
  i = g_hash_table_lookup(fp->chain_hash, fp->chain_state);
  assert (i);
  i->intervals += sample_count (i->start_t, fp->tube->t);
  i->times += (fp->tube->t - i->start_t);
  fp->chain_state = NULL;
}

void undo_state_off_indicator(flake *fp) {
  unsigned char *h;
  chain_state_record *i;

  //printf("Turning off chain state %d for flake %p.\n",fp->chain_state,fp);
  h = hash_assembly (fp->cell, 1<<(fp->tube->P));
  i = g_hash_table_lookup(fp->chain_hash, h);
  assert (i);
  i->intervals -= sample_count (i->old_start_t,i->start_t);
  i->times -= i->start_t - i->old_start_t;
  i->start_t = i->old_start_t;
  free (h);
}

int test_detailed_balance (tube *tp, indicator_data *data) {
  int i,j;
  flake *flake;
  detailed_balance_data d;
  double correct_state_ratio;
  double r, est_mean, est_variance, confidence_interval;
  double confidence_ratio;
  int n;
  int unbalanced_states = 0;

  printf("\n********************************************************\n");
  printf("Entering stage three: testing whether detailed balance \n"
	 "has been achieved.\n");
  printf("*******************************************************\n\n");
  d.means = (double *) malloc_err (sizeof(double) * tp->chains);
  d.variances = (double *) malloc_err (sizeof(double) * tp->chains);
  d.covariances = (double **) malloc_err (sizeof(double *) * tp->chains);
  for (i = 0; i < tp->chains; i++) {
    d.covariances[i] = (double *) malloc_err (sizeof(double) * tp->chains);
  }

  for (i = 0; i < tp->chains; i++) {
    d.means[i] = data[i].mean_of_means;
    printf("Intervals spent in state %d is %d out of %d intervals total.\n",
	 i,(int) (d.means[i]*(tp->t/sampling_rate)/CHOP),(int) ((tp->t/sampling_rate*tp->chains)/CHOP));
    printf("Time spent in state %d is %f seconds out of %f total.\n",
	   i, data[i].total_time, tp->t*tp->chains);
  }

  n = (floor((tp->t / sampling_rate)/CHOP)  - 1) * tp->chains;  
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
      if (r - confidence_interval > correct_state_ratio) {
	printf("\nTrue ratio of variables %d and %d is %1.2e.  ",j,i,correct_state_ratio);
	printf("Simulated ratio is %1.2e.\nThe computed confidence interval is %2.1f%%.\n",
	       r,100*confidence_ratio);
	printf("Estimated ratio is too high.\n");
	unbalanced_states++;
      }
      if (r + confidence_interval < correct_state_ratio) {
	printf("\nTrue ratio of variables %d and %d is %1.2e.  ",j,i,correct_state_ratio);
	printf("Simulated ratio is %1.2e.\nThe computed confidence interval is %2.1f%%.\n",
	       r,100*confidence_ratio);
	printf("Estimated ratio is too low.\n");
	unbalanced_states++;
      }
      /*
      if (data[i].total_time/data[j].total_time - confidence_interval > correct_state_ratio) {
	printf("\nTrue ratio of variables %d and %d is %1.2e.  ",j,i,correct_state_ratio);
	printf("Simulated ratio is %1.2e.\nThe computed confidence interval is %2.1f%%.\n",
	       data[i].total_time/data[j].total_time,100*confidence_ratio);
	printf("Estimated ratio is too high.\n");
	unbalanced_states++;
      }
      if (data[i].total_time/data[j].total_time + confidence_interval < correct_state_ratio) {
	printf("\nTrue ratio of variables %d and %d is %1.2e.  ",j,i,correct_state_ratio);
	printf("Simulated ratio is %1.2e.\nThe computed confidence interval is %2.1f%%.\n",
	       data[i].total_time/data[j].total_time,100*confidence_ratio);
	printf("Estimated ratio is too low.\n");
	unbalanced_states++;
      }
      */
    }
  }
  return unbalanced_states;
}

void run_xgrow_tests (tube *tp,double Gmc, double Gse, int si, int sj, int sn, int size) {
  indicator_data *data;
  int unbalanced_states; 

  block_time = 0.005 * exp(Gmc);
  sampling_rate = 5e-7 * exp(Gmc);
  printf("\nAttention!\n\n");
  printf("xgrow is being run in testing mode.  In testing mode, xgrow will \n"
	 "attempt to see if xgrow's simulation algorithm, when it reaches a \n"
	 "state close to equilibrium, satisfies detailed balance on the tileset \n"
	 "given as the argument.  These tests will allow you to identify any bugs \n"
	 "in which some states are incorrectly weighted over other states, skewing \n"
	 "the results of your simulations.  However, since reaching a state close to \n"
	 "equilibrium is required, you should only run tests on a \n"
	 "tile set and physical conditions where this can be achieved in \n"
	 "a realistic simulation.\n\n");
  printf("Testing consists of three stages.\n\n");
  printf("In the first stage, xgrow will run a series of anneals from above the \n"
	 "melting temperature of the tile set down to the physical conditions \n"
	 "set by the tile set or xgrow defaults.  At the end of each \n"
	 "anneal, xgrow will choose %d representative states to check for \n"
	 "detailed balance. It will run %d different anneals in order to find a \n"
	 "a representative set of states to follow in the next stage.  \n"
	 "The goal of the next stages will be to check whether the ratio \n"
	 "between the amount of time spent in these states during simulations \n"
	 "satisfies detailed balance.\n\n",STATES_TO_ADD_PER_ANNEAL,
	 CHAIN_COUNT/STATES_TO_ADD_PER_ANNEAL);
  printf("In the second stage, xgrow will run several simulations in parallel \n"
	 "until the time spent in each of the representative samples becomes \n"
	 "approximately equal.  Specifically, xgrow will follow %d simulations.\n"
	 "It will run them for a block of %f seconds, then compute an approximate \n"
	 "distance to equilibrium by computing a scaling factor (Gelman, Rubin 1992) \n"
	 "for each state.  A scaling factor is the approximate difference in the \n"
	 "amount of time spent in a state we are following between chains.  At \n"
	 "equilibrium, we expect the sampled state to appear equally in all of \n"
	 "the chains. Thus, we run until the sampled states are seen with \n"
	 "approximately the same frequency.  At this point, xgrow will enter the \n"
	 "third stage.\n\n",CHAIN_COUNT,block_time);
  printf("In the third stage, xgrow will compare the ratios of times spent in \n"
	 "each representative state, recorded in the second stage, and compute \n"
	 "the delta G of each of the representative states.  It will use this \n"
	 "information to check whether these ratios satisfy detailed balance. \n"
	 "If any of the ratios determined by simulation are not a reasonable \n"
	 "statistical approximation for the ratios of the delta G's, xgrow will \n"
	 "notify you.  Occasional statistical outliers are to be expected, but a \n"
	 "persistent problem, especially among a particular set of states, may \n"
	 "indicate a bug in xgrow.\n\n");
  printf("Press enter to start testing.\n");
  fscanf(stdin,"%*c");

  seed_i = si;
  seed_j = sj;
  seed_n = sn;
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
  unbalanced_states = test_detailed_balance (tp, data);
  if (!unbalanced_states) {
    printf("\n\n***************************************************\n");
    printf("Detailed balance seems to have been acheived.\n");
    printf("***************************************************\n");
  }
  else {
    printf("\n\n************************************************************\n");
    printf("There was a problem.  This may be a statistical error, \n"
	   "but if it repeats, there may be a bug.\n");
    printf("************************************************************\n");
  }
  
}
