

#define DEBUG 1
#define dprintf if (DEBUG) printf
#define d2printf if (DEBUG==2) printf

double drand48(); long lrand48(); 
double exp(); double log();

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* make off-by-one error less likely : include a boundary of empty */
/* cells that will never be modified                               */
#define Cell(i,j) cell[(i)+1][(j)+1]
/* note i,j here are indexed 0 <= i,j < (1<<P)                     */
/* 
   for periodic boundary conditions,  
   cell[0][*] cell[size+1][*] cell[*][0] cell[*][size+1]
   are maintained with the contents of the opposite wall;
   for non-periodic conditions, they are maintained as 0=empty
*/

typedef struct flake_struct {
  int **units;         /* tile types */
  double *strength;    /* bond strengths */
  int num_bindings;    /* length of strength */
  
  unsigned char N, P;  /* # non-empty tile types; 2^P active cell grid     */
  unsigned char hydro; /* does this tile set use hydrolysis rules?         */
                       /* in this case, N must be even, and tiles          */
                       /* 1...N/2 are non-hydrolized; tiles N/2+1...N are  */
                       /* hydrolized, with weaker sticky end strengths     */
  unsigned char **cell;/* tile type at [i][j]; array of arrays             */
                       /* note 0 <= i,j <= 2^P+1, allowing for borders     */
  double ***rate;      /* hierarchical rates for events in non-empty cells */
                       /* rate[p][i][j] has 0 <= i,j <= 2^p                */
                       /* rate[P][i][j] = sum rates for Cell(i,j)          */
                       /* rate[p][i][j] =   sum rate[p+1][2*i+di][2*j+dj]  */
                       /*                 (di,dj in {0,1})                 */
                       /* rate[0][0][0] + k * sum conc[] = net event rate  */
  int ***empty;        /* hierarchical tally of number of empty cells      */
                       /* adjacent to some non-empty cell                  */
  double k;            /* forward rate constant for on-events.             */
  double kas,kao,      /* f.r.c (ratio to k) for "hydrolysis" spontaneous, */
         kam,kae,kah;  /* and when input se are mismatched, empty, or      */
                       /* input tiles are "hydrolyized". kao is the ratio  */
                       /* of output-triggered/input-triggered rates.       */
  double *conc;        /* concentration of each tile type n = 1...N        */
                       /* conc[0] = sum conc[n] from n = 1...N             */
                       /* conc[n] =def= exp(-Gmc[n])                       */
  double flake_conc;   /* for depleting concentrations of monomers as they */
                       /* are incorporated; if 0, then no depletion occurs */
  double *Gcb;         /* neg std free energy of chemical bonds in tile    */
                       /* (zero typically; pos for hydrolysis models)      */
  double **Gse_EW;     /* bond energy between east-west oriented tiles     */
                       /* Gse_EW[n1][n2] =  n1's west se  : n2's east se   */
  double **Gse_NS;     /* Gse_NS[n1][n2] =  n1's south se : n2's north se  */
                       /* NOTE n could be empty tile, for which Gse = 0    */
                       /* coordinate system is: i+,j+ moves S,E            */
  double G;            /* cumulative energy of system                      */
  double t;            /* cumulative time in seconds                       */
  int seed_i,seed_j;   /* special site which won't change                  */
  unsigned char seed_n; 
  int events;          /* cumulative number of events                      */
  int stat_a,stat_d,   /* tally of number of association, dissociation,    */
      stat_h,stat_f;   /* "hydrolysis", and "fission" events               */
  int mismatches;      /* number of se edges that don't agree              */
  double *rv;          /* scratch space, size fp->N+1                      */
  int *Fnext, *Fgroup; /* size x size array for fill scratch space         */
} flake;          

extern int periodic;    /* simulation on torus */
extern int wander;      /* of seed tile designation */

flake *init_flake(unsigned char P, unsigned char N, int num_bindings);
void free_flake(flake *fp);
void clean_flake(flake *fp);
void set_params(flake *fp, int** units, double* strength, double* relconc,
 int hydro, double k, double Gmc, double Gse,
 double Gmch, double Gseh, double Ghyd, 
 double Gas, double Gam, double Gae, double Gah, double Gao,
 double Gfc);
void change_cell(flake *fp, int i, int j, unsigned char n);
void simulate(flake *fp, int events, double tmax, int emax, int smax);
void linear_simulate( double ratek, double Gmc, double Gse,
                      double tmax, int emax, int smax);


