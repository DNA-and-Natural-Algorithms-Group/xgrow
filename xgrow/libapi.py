import cffi
import numpy as np
import os, re

_XGROW_LIB = re.sub(r" ", r"\ ", os.path.join(os.path.dirname(__file__), "libxgrow.so"))

ffi = cffi.FFI()

ffi.cdef(
    """ /* CDEF copy-paste on 2021-07-06 */
typedef struct flake_struct {
   struct tube_struct *tube; /* contains tile set, reaction conditions,     */
   /* time, event stats, scratch space...              */
   /* all flakes are the same size, 2^(tube->P)        */
   unsigned int N, P; /* # non-empty tile types; 2^P active cell grid     */

   unsigned int *cell; /* tile type at [i][j]; array of arrays             */
   /* note 0 <= i,j <= 2^P+1, allowing for borders     */
   double ***rate; /* hierarchical rates for events in non-empty cells */
   /* rate[p][i][j] has 0 <= i,j < 2^p                 */
   /* rate[P][i][j] = sum rates for Cell(i,j)          */
   /* rate[p][i][j] =   sum rate[p+1][2*i+di][2*j+dj]  */
   /*                 (di,dj in {0,1})                 */
   /* rate[0][0][0] + k * sum conc[] = net event rate  */
   int ***empty; /* hierarchical tally of number of empty cells      */
   /* adjacent to some non-empty cell.                 */
   /* for irreversible model, counts only if there is  */
   /* a tile type that could stick here.               */
   double flake_conc; /* for depleting concentrations of monomers as they */
   /* are incorporated; if 0, then no depletion occurs */
   double G;           /* cumulative energy of tile flake                  */
   int seed_i, seed_j; /* special site which won't change                  */
   unsigned int seed_n;
   unsigned long long events;                    /* total on, off, hydrolysis events in this flake   */
   int tiles;                       /* total number of tiles in this flake              */
   int seed_is_double_tile;         /* If the seed is a double tile, it will be a monomer,
                                       but the number of tiles will be reported as 2.  */
   int seed_is_vdouble_tile;        /* same for vdoubles */
   int mismatches;                  /* number of se edges that don't agree              */
   struct flake_struct *next_flake; /* for NULL-terminated linked list     */
   struct flake_tree_struct *tree_node; /* for tree of flakes              */
   int *is_present;                     /* records whether each of the watched
                                           tile types are present              */

   int flake_ID;     /* which flake is this (for display use only)       */
   void *chain_hash; /* When flake has visited particular states;        */
   /* used for testing purposes                        */
   unsigned char *chain_state; /* If we're currently at a configuration that has an*/
   /* indicator variable, the hash code for that state */
   int periodic; /* whether canvasa is periodic: don't want to refer to tube */

} flake;

typedef struct flake_tree_struct {
   struct flake_tree_struct *left, *right, *up; /* up==NULL iff root     */
   flake *fp;                                   /* if fp != NULL, left==right==NULL   */
   int empty;
   double rate; /* analogous to empty & rate in flake */
} flake_tree;

typedef struct assembly_list_struct {
   struct assembly_list_struct *next;
   unsigned int **assembly;
} assembly_list;

typedef struct tube_struct {
   int *tileb;       /* {N E S W} bond types for each tile type */
   double *strength; /* bond strengths.  assumes tile types stick to
                        each other and not other types of tiles          */
   double *glue;     /* A generalized version of the strength function:
                         The level of "glue" between various bond types   */
   int num_bindings; /* length of strength */
   int *dt_right;    /* The right half of a double tile.
                        If it doesn't have one, or the tile is the
                        left half of a double tile, the value
                        here is 0 */
   int *dt_left;     /* The left half of a double tile.
                        If it doesn't have one, or the tile is the
                        right half of a double tile, the value
                        here is 0 */
   int *dt_down;     /* dt_down and dt_up are analagous to dt_right and dt_left
                      for vertical tiles. */
   int *dt_up;
   int double_tile_count;
   int vdouble_tile_count;
   double tinybox;     /* If this value is nonzero, indicates that
                          each kind of two tile flake should be
                          dynamically created at a rate
                          tinybox*k_f*tp->conc[0]^2 */
   double anneal_g,    /* Used if annealing is on.  If so, adjust Gse over */
       anneal_t;       /* time, with time constant anneal_t                */
   double Gse_final;   /* Gse to approach asymptotically in anneal         */
   double update_freq; /* Number of times to update the interval per time
                          constant                                         */
   int updates;        /* Number of updates that have taken place          */
   double anneal_h;    /* Next 5 are variables for linear anneal           */
   double anneal_s;
   double startC;
   double currentC;
   double endC;
   double seconds_per_C;
   double Gse;           /* Current Gse                                      */
   double Gmc;           /* Gmc                                              */
   double next_update_t; /* Precompute next update time                   */
   unsigned int N, P;            /* # non-empty tile types; 2^P active cell grid     */

   int num_flakes;         /* how many flakes do we have here?                 */
   int total_flakes;       /* how many flakes have we made, total
                              (in tinybox may not be the same as num_flakes)   */
   int largest_flake;      /* id of largest flake                              */
   int largest_flake_size; /* size of largest flake                          */
   flake *flake_list;      /* for NULL-terminated linked list                  */
   flake_tree *flake_tree; /* binary tree for fast event selection          */
   int default_seed_i,     /* The last seed_i stated, which is used in creating
                              new flakes */
       default_seed_j;
   double initial_Gfc;
   unsigned int hydro; /* does this tile set use hydrolysis rules?         */
   /* in this case, N must be even, and tiles          */
   /* 1...N/2 are non-hydrolized; tiles N/2+1...N are  */
   /* hydrolized, with weaker sticky end strengths     */
   double T; /* threshold T for irreversible Tile Assembly Model */
   /* used to prevent dissociation only if T>0         */
   /* also, prevents incorrect association             */
   /* on-rates & off-rates are calculated as usual,    */
   /* but events violating the model are discarded     */
   double Gseh;
   double Gmch;
   double k;          /* forward rate constant for on-events.             */
   double kas, kao,   /* f.r.c (ratio to k) for "hydrolysis" spontaneous, */
       kam, kae, kah; /* and when input se are mismatched, empty, or      */
   /* input tiles are "hydrolized". kao is the ratio   */
   /* of output-triggered/input-triggered rates.       */
   double *conc; /* concentration of each tile type n = 1...N        */
   /* conc[0] = sum conc[n] from n = 1...N             */
   /* conc[n] =def= exp(-Gmc[n])                       */
   double *Gcb; /* neg std free energy of chemical bonds in tile    */
   /* (zero typically; pos for hydrolysis models)      */
   double **Gse_EW; /* bond energy between east-west oriented tiles     */
   /* Gse_EW[n1][n2] =  n1's west se  : n2's east se   */
   double **Gse_NS; /* Gse_NS[n1][n2] =  n1's south se : n2's north se  */
   /* NOTE n could be empty tile, for which Gse = 0    */
   /* coordinate system is: i+,j+ moves S,E            */
   double t;             /* cumulative time in seconds                       */
   unsigned long long events;         /* cumulative number of events                      */
   unsigned long long stat_a, stat_d, /* tally of number of association, dissociation,  */
       stat_h, stat_f;   /* "hydrolysis", and "fission" events               */
   int stat_m;
   int ewrapped;        /* has the event counter wrapped around?            */
   double *rv;          /* scratch space, size fp->1+N+4 (for chunk_fission)*/
   int *Fnext, *Fgroup; /* size x size array for fill scratch space         */
   int all_present;     /* True if all the tiles in untiltiles are in the assembly */
   int untiltilescount;
   /* Testing variables */
   int watching_states;      /* true if we are testing xgrow, false otherwise */
   int chains;               /* Number of chains that we are going to follow for
                                testing purposes */
   void *chain_states;       /* Hash of assemblies that are indicator variables;
                                use for testing purposes */
   int tracking_seen_states; /* tracking states that are seen             */
   int states_seen_count;    /* For use in testing --find a time at which */
   void *states_seen_hash;   /* hash that records which states we have been to. */
   unsigned int ***start_states;
   double *start_state_Gs;
   // Below was previously an extern to xgrow.c
   int periodic;           /* simulation on torus */
   int wander;             /* of seed tile designation */
   int fission_allowed;    /* allow dissociation that breaks flake in two? */
   int zero_bonds_allowed; /* allow association of tiles that make only 0 strength bonds
                              to flake? */
   double blast_rate_alpha;
   double blast_rate_beta;
   double blast_rate_gamma;
   double blast_rate;
   double min_strength;
   int *present_list;
   int present_list_len;
   int untiltiles;
   // End of former externs
} tube;

void setup_tube(tube *tp);
tube *init_tube(unsigned int P, unsigned int N, int num_bindings);
flake *init_flake(unsigned int P, unsigned int N, int seed_i, int seed_j, int seed_n, double Gfc,
                  int present_list_len, int periodic);
flake *free_flake(flake *fp);
void free_tube(tube *tp);
void set_Gses(tube *tp, double Gse, double Gseh);
void insert_flake(flake *fp, tube *tp);
void print_tree(flake_tree *ftp, int L, char s);
void clean_flake(flake *fp, double X, int iters);
void fill_flake(flake *fp, double X, int iters);
void error_radius_flake(flake *fp, double rad);
void repair_flake(flake *fp, double T, double Gse);
void reset_params(tube *tp, double old_Gmc, double old_Gse, double new_Gmc,
                  double new_Gse, double Gseh);
void recalc_G(flake *fp);
double calc_dG_bonds(flake *fp);
int calc_perimeter(flake *fp);
void update_all_rates(tube *tp);
void update_rates(flake *fp, int ii, int jj);
void update_tube_rates(flake *fp);
void change_cell(flake *fp, int i, int j, unsigned int n);
void change_seed(flake *fp, int new_i, int new_j);
int flake_fission(flake *fp, int i, int j);
void simulate(tube *tp, unsigned long long events, double tmax, int emax, int smax, int fsmax,
              int smin, int mmax);
void linear_simulate(double ratek, double Gmc, double Gse, double tmax, int emax,
                     int smax, int mmax);
int Mism(tube *tp, flake *fp, int size, int i, int j, int n);
"""
)

grow = ffi.dlopen(_XGROW_LIB)


class Tube:
    def __init__(self, size, ntiles, nbonds) -> None:
        self.sizeP = int(np.ceil(np.log2(size)))
        self.size = 2 ** self.sizeP
        self._inner = grow.init_tube(self.sizeP, ntiles, nbonds)

        self._tileb = np.frombuffer(
            ffi.buffer(
                self._inner.tileb, (ntiles + 1) * 4 * ffi.sizeof("unsigned int")
            ),
            dtype=np.uint32,
        ).reshape((-1, 4))
        self._glues = np.frombuffer(
            ffi.buffer(
                self._inner.glue, (nbonds + 1) * (nbonds + 1) * ffi.sizeof("double")
            ),
            dtype=np.double,
        ).reshape(((nbonds + 1), (nbonds + 1)))
        self._strengths = np.frombuffer(
            ffi.buffer(self._inner.strength, (nbonds + 1) * ffi.sizeof("double")),
            dtype=np.double,
        )
        self._concs = np.frombuffer(
            ffi.buffer(self._inner.conc, (ntiles + 1) * ffi.sizeof("double")),
            dtype=np.double,
        )

    @property
    def tileb(self) -> np.ndarray:
        return self._tileb

    @tileb.setter
    def tileb(self, v: np.ndarray):
        self._tileb[:] = v

    @property
    def concs(self) -> np.ndarray:
        return self._concs

    @concs.setter
    def concs(self, v: np.ndarray):
        self._concs[:] = v

    @property
    def strengths(self) -> np.ndarray:
        return self._strengths

    @strengths.setter
    def strengths(self, v: np.ndarray):
        self._strengths[:] = v

    @property
    def glues(self) -> np.ndarray:
        return self._glues

    @glues.setter
    def glues(self, v: np.ndarray):
        self._glues[:] = v

    @property
    def Gse(self) -> float:
        return self._inner.Gse

    @Gse.setter
    def Gse(self, v: float):
        self._inner.Gse = v

    @property
    def Gmc(self) -> float:
        return self._inner.Gmc

    @Gmc.setter
    def Gmc(self, v: float):
        self._inner.Gmc = v

    @property
    def k(self) -> float:
        return self._inner.k

    @k.setter
    def k(self, v: float):
        self._inner.k = v

    @property
    def T(self) -> int:
        return self._inner.T

    @T.setter
    def T(self, v: int):
        self._inner.T = v

    def setup(self):
        grow.setup_tube(self._inner)

    # FIXME: events should be long long
    def simulate(
        self,
        events: int,
        tmax: float = 0.0,
        emax: int = 0,
        smax: int = 0,
        fsmax: int = 0,
        smin: int = -1,
        mmax: int = 0,
    ):
        grow.simulate(self._inner, events, tmax, emax, smax, fsmax, smin, mmax)

    def insert_flake(self, flake: "Flake"):
        grow.insert_flake(flake._inner, self._inner)


class Flake:
    def __init__(self, size, ntiles, i, j, n) -> None:
        self.sizeP = int(np.ceil(np.log2(size)))
        self.size = 2 ** self.sizeP
        self._inner = grow.init_flake(self.sizeP, ntiles, i, j, n, 0, 0, 0)
        self._cell = np.frombuffer(
            ffi.buffer(
                self._inner.cell, (2 + self.size) ** 2 * ffi.sizeof("unsigned int")
            ),
            dtype=np.uint32,
        ).reshape((2 + self.size), (2 + self.size))

    @property
    def cell(self) -> np.ndarray:
        return self._cell

    def change_cell(self, i: int, j: int, t: int):
        grow.change_cell(self._inner, i, j, t)

    @property
    def tiles(self) -> int:
        return self._inner.tiles

    @tiles.setter
    def tiles(self, v: int):
        self._inner.tiles = v
