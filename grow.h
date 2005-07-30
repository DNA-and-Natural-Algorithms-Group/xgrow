/* grow.h

  This code is freely distributable.

  by Erik Winfree
*/


#ifndef __GROW_H__
#define __GROW_H__


#define DEBUG 1
#define dprintf if (DEBUG) printf
#define d2printf if (DEBUG==2) printf

double drand48(); long lrand48(); 
double exp(); double log();

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif
/* make off-by-one error less likely : include a boundary of empty */
/* cells that will never be modified                               */
#define Cell(i,j) cell[(i)+1][(j)+1]
/* note i,j here are indexed 0 <= i,j < (1<<P)                     */
/* 
   for periodic boundary conditions, where size=2^P,  
    cell[0][*] cell[size+1][*] cell[*][0] cell[*][size+1]
   are maintained with the contents of the opposite wall;
   for non-periodic conditions, they are maintained as 0=empty

  Thus, it is meaningful to request 
    Cell(-1,j), Cell(size,j), Cell(i,-1), Cell(i,size)
*/

/* for times when it's inconvenience to know if i,j are within bounds */
#define CellM(i,j) cell[periodic?((i+size)%size):MAX(0,MIN((i)+1,size+1))][periodic?((j+size)%size):MAX(0,MIN((j)+1,size+1))]

/* macro definition of summed sticky end bond energy                    */
/* computes energy IF Cell(i,j) were n, given its current neighbors     */
/* assumes "fp" arg is a simple variable, but others can be expressions */
/* note that n != 0  and assumes 0 <= i,j < (1<<fp->P)                  */
#define Gse(fp,i,j,n) (                                \
        fp->tube->Gse_EW[ n ] [ fp->Cell(i,(j)-1) ] +  \
        fp->tube->Gse_EW[ fp->Cell(i,(j)+1) ] [ n ] +  \
        fp->tube->Gse_NS[ n ] [ fp->Cell((i)+1,j) ] +  \
        fp->tube->Gse_NS[ fp->Cell((i)-1,j) ] [ n ] )

/* definition for total sticky-end strength around a pair or a 2x2 chunk */
/* -- note that if some site is empty, this still gives the correct      */
/*    energy for dissociation                                            */
/* Here, 0 <= i,j < 2^P if periodic, and 0 <= i,j < 2^P-1 otherwise.     */
#define chunk_Gse_EW(fp,i,j,n) ( \
        Gse(fp,i,j,n) +                                                         \
        Gse(fp,i,((j)+1)%size,fp->Cell(i,((j)+1)%size)) -                       \
        2*fp->tube->Gse_EW[ fp->Cell(i,(j)+1) ] [ n ] ) 
#define chunk_Gse_NS(fp,i,j,n) ( \
        Gse(fp,i,j,n) +                                                         \
        Gse(fp,((i)+1)%size,j,fp->Cell(((i)+1)%size,j)) -                       \
        2*fp->tube->Gse_NS[ n ] [ fp->Cell((i)+1,j) ] ) 
#define chunk_Gse_2x2(fp,i,j,n) ( \
        Gse(fp,i,j,n) +                                                         \
        Gse(fp,i,((j)+1)%size,fp->Cell(i,((j)+1)%size)) +                       \
        Gse(fp,((i)+1)%size,j,fp->Cell(((i)+1)%size,j)) +                       \
        Gse(fp,((i)+1)%size,((j)+1)%size,fp->Cell(((i)+1)%size,((j)+1)%size)) - \
        2*fp->tube->Gse_EW[ fp->Cell(i,(j)+1) ] [ fp->Cell(i,j) ] -             \
        2*fp->tube->Gse_NS[ fp->Cell(i,j) ] [ fp->Cell((i)+1,j) ] -             \
        2*fp->tube->Gse_EW[ fp->Cell((i)+1,(j)+1) ] [ fp->Cell((i)+1,j) ] -     \
        2*fp->tube->Gse_NS[ fp->Cell(i,(j)+1) ] [ fp->Cell((i)+1,(j)+1) ]  ) 


/* similar definition to count the number of sides that are mismatched   */
/* also,  n != 0   and assumes 0 <= i,j < (1<<fp->P).                    */
/* this gives the number of mismatched bonds (not null bonds, not equal) */
/* NOT the number of tiles with mismatches.                              */
/* note that this counts twice: sum_i,j Mism(i,j) == 2ce # mism. bonds.  */
/* however, *if* the sum is accumulated during assembly, exactly when    */
/* the tile at i,j is being added, then it counts each mismatch ONCE.    */

#define Mism(fp,i,j,n) (                                                   \
 ((fp->tube->tileb)[n][1] != (fp->tube->tileb)[fp->Cell(i,(j)+1)][3] &&    \
  (fp->tube->tileb)[n][1]*(fp->tube->tileb)[fp->Cell(i,(j)+1)][3] > 0) +   \
 ((fp->tube->tileb)[n][3] != (fp->tube->tileb)[fp->Cell(i,(j)-1)][1] &&    \
  (fp->tube->tileb)[n][3]*(fp->tube->tileb)[fp->Cell(i,(j)-1)][1] > 0) +   \
 ((fp->tube->tileb)[n][2] != (fp->tube->tileb)[fp->Cell((i)+1,j)][0] &&    \
  (fp->tube->tileb)[n][2]*(fp->tube->tileb)[fp->Cell((i)+1,j)][0] > 0) +   \
 ((fp->tube->tileb)[n][0] != (fp->tube->tileb)[fp->Cell((i)-1,j)][2] &&    \
  (fp->tube->tileb)[n][0]*(fp->tube->tileb)[fp->Cell((i)-1,j)][2] > 0) )



typedef struct flake_struct {
  struct tube_struct *tube; /* contains tile set, reaction conditions,     */
                       /* time, event stats, scratch space...              */
                       /* all flakes are the same size, 2^(tube->P)        */
  unsigned char N, P;  /* # non-empty tile types; 2^P active cell grid     */

  unsigned char **cell;/* tile type at [i][j]; array of arrays             */
                       /* note 0 <= i,j <= 2^P+1, allowing for borders     */
  double ***rate;      /* hierarchical rates for events in non-empty cells */
                       /* rate[p][i][j] has 0 <= i,j < 2^p                 */
                       /* rate[P][i][j] = sum rates for Cell(i,j)          */
                       /* rate[p][i][j] =   sum rate[p+1][2*i+di][2*j+dj]  */
                       /*                 (di,dj in {0,1})                 */
                       /* rate[0][0][0] + k * sum conc[] = net event rate  */
  int ***empty;        /* hierarchical tally of number of empty cells      */
                       /* adjacent to some non-empty cell.                 */
                       /* for irreversible model, counts only if there is  */
                       /* a tile type that could stick here.               */  
  double flake_conc;   /* for depleting concentrations of monomers as they */
                       /* are incorporated; if 0, then no depletion occurs */
  double G;            /* cumulative energy of tile flake                  */
  int seed_i,seed_j;   /* special site which won't change                  */
  unsigned char seed_n; 
  long int events;     /* total on, off, hydrolysis events in this flake   */
  int tiles;           /* total number of tiles in this flake              */
  int seed_is_double_tile;  /* If the seed is a double tile, it will be a monomer,
				but the number of tiles will be reported as 2.  */
  int mismatches;      /* number of se edges that don't agree              */
  struct flake_struct *next_flake;  /* for NULL-terminated linked list     */
  struct flake_tree_struct *tree_node;  /* for tree of flakes              */
  int flake_ID;        /* which flake is this (for display use only)       */
  void *chain_hash;    /* When flake has visited particular states;        */
		       /* used for testing purposes                        */
  int chain_state;     /* If we're currently at a configuration that has an*/
		       /* indicator variable, the hash code for that state */

} flake;          

typedef struct flake_tree_struct {
  struct flake_tree_struct *left, *right, *up;    /* up==NULL iff root     */
  flake *fp;                         /* if fp != NULL, left==right==NULL   */
  int empty; double rate;            /* analogous to empty & rate in flake */
} flake_tree;

typedef struct assembly_list_struct {
  struct assembly_list_struct *next;
  unsigned char **assembly;
} assembly_list;

typedef struct tube_struct {
  int **tileb;         /* {N E S W} bond types for each tile type */
  double *strength;    /* bond strengths.  assumes tile types stick to
			  each other and not other types of tiles          */
  double **glue;       /* A generalized version of the strength function:
			  The level of "glue" between various bond types   */
  int num_bindings;    /* length of strength */
  int *dt_right;     /* The right half of a double tile.  
			If it doesn't have one, or the tile is the  
			left half of a double tile, the value 
			here is 0 */
  int *dt_left;     /* The left half of a double tile.  
			If it doesn't have one, or the tile is the  
			right half of a double tile, the value 
			here is 0 */
  int tinybox;         /* If this value is nonzero, indicates that
			  each kind of two tile flake should be
			  dynamically created at a rate
			  tinybox*k_f*e^{-Gmc} */
  double anneal_g,     /* Used if annealing is on.  If so, adjust Gse over */
     anneal_t;         /* time, with time constant anneal_t                */  
  double Gse_final;    /* Gse to approach asymptotically in anneal         */
  double update_freq;     /* Number of times to update the interval per time 
			  constant                                         */
  int updates;         /* Number of updates that have taken place          */
  double Gse;          /* Current Gse                                      */
  double Gmc;          /* Gmc                                              */
  double next_update_t;   /* Precompute next update time                   */
  unsigned char N, P;  /* # non-empty tile types; 2^P active cell grid     */

  int num_flakes;      /* how many flakes do we have here?                 */
  flake *flake_list;   /* for NULL-terminated linked list                  */
  flake_tree *flake_tree; /* binary tree for fast event selection          */

  unsigned char hydro; /* does this tile set use hydrolysis rules?         */
                       /* in this case, N must be even, and tiles          */
                       /* 1...N/2 are non-hydrolized; tiles N/2+1...N are  */
                       /* hydrolized, with weaker sticky end strengths     */
  double T;            /* threshold T for irreversible Tile Assembly Model */
                       /* used to prevent dissociation only if T>0         */
                       /* also, prevents incorrect association             */
                       /* on-rates & off-rates are calculated as usual,    */
                       /* but events violating the model are discarded     */
  double k;            /* forward rate constant for on-events.             */
  double kas,kao,      /* f.r.c (ratio to k) for "hydrolysis" spontaneous, */
         kam,kae,kah;  /* and when input se are mismatched, empty, or      */
                       /* input tiles are "hydrolized". kao is the ratio   */
                       /* of output-triggered/input-triggered rates.       */
  double *conc;        /* concentration of each tile type n = 1...N        */
                       /* conc[0] = sum conc[n] from n = 1...N             */
                       /* conc[n] =def= exp(-Gmc[n])                       */
  double *Gcb;         /* neg std free energy of chemical bonds in tile    */
                       /* (zero typically; pos for hydrolysis models)      */
  double **Gse_EW;     /* bond energy between east-west oriented tiles     */
                       /* Gse_EW[n1][n2] =  n1's west se  : n2's east se   */
  double **Gse_NS;     /* Gse_NS[n1][n2] =  n1's south se : n2's north se  */
                       /* NOTE n could be empty tile, for which Gse = 0    */
                       /* coordinate system is: i+,j+ moves S,E            */
  double t;            /* cumulative time in seconds                       */
  long int events;     /* cumulative number of events                      */
  long int stat_a,stat_d,/* tally of number of association, dissociation,  */
      stat_h,stat_f;   /* "hydrolysis", and "fission" events               */
  int ewrapped;        /* has the event counter wrapped around?            */
  double *rv;          /* scratch space, size fp->1+N+4 (for chunk_fission)*/
  int *Fnext, *Fgroup; /* size x size array for fill scratch space         */
  /* Testing variables */
  int watching_states; /* true if we are testing xgrow, false otherwise */
  int chains;          /* Number of chains that we are going to follow for 
			  testing purposes */
  void *chain_states;  /* Hash of assemblies that are indicator variables; 
			  use for testing purposes */
  int tracking_seen_states;     /* tracking states that are seen             */
  int states_seen_count;    /* For use in testing --find a time at which */
  void  *states_seen_hash;   /* hash that records which states we have been to. */
  unsigned char ***start_states;
  double *start_state_Gs;

} tube;          

extern int periodic;    /* simulation on torus */
extern int wander;      /* of seed tile designation */
extern int fission_allowed; /* allow dissociation that breaks flake in two? */
extern int zero_bonds_allowed; /* allow association of tiles that make only 0 strength bonds to flake? */
extern double blast_rate_alpha;
extern double blast_rate_beta;
extern double blast_rate_gamma;
extern double blast_rate;


tube *init_tube(unsigned char P, unsigned char N, int num_bindings);
flake *init_flake(unsigned char P, unsigned char N,
 int seed_i, int seed_j, int seed_n, double Gfc);
flake *free_flake(flake *fp);
void free_tube(tube *tp);
void insert_flake(flake *fp, tube *tp);
void print_tree(flake_tree *ftp, int L, char s);
void clean_flake(flake *fp, double X, int iters);
void fill_flake(flake *fp, double X, int iters);
void error_radius_flake(flake *fp, double rad);
void repair_flake(flake *fp, double T, double Gse);
void set_params(tube *tp, int** tileb, double* strength, double **glue, 
 double* stoic, double anneal_g, double anneal_t, int updates_per_RC,int *dt_right, int *dt_left, int hydro, double k, double Gmc, double Gse,
 double Gmch, double Gseh, double Ghyd, 
 double Gas, double Gam, double Gae, double Gah, double Gao, double T, int tinybox);
void reset_params(tube *tp, double old_Gmc, double old_Gse,
 double new_Gmc, double new_Gse, double Gseh);
void recalc_G(flake *fp);
double calc_dG_bonds(flake *fp);
int calc_perimeter(flake *fp);
void update_rates(flake *fp, int ii, int jj);
void update_tube_rates(flake *fp);
void change_cell(flake *fp, int i, int j, unsigned char n);
void change_seed(flake *fp, int new_i, int new_j);
int flake_fission(flake *fp, int i, int j);
void simulate(tube *tp, int events, double tmax, int emax, int smax);
void linear_simulate( double ratek, double Gmc, double Gse,
                      double tmax, int emax, int smax);

#endif /* ifdef __GROW_H__ */



