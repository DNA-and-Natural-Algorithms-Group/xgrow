/*  xgrow.c

    This code is freely distributable.

History:
Originally X-Automatalab by Michael Creutz creutz@wind.phy.bnl.gov
modified for excitable media by Erik Winfree Nov 1995
modifications mostly in parallelupdate()
modified again for DNA self-assembly simulations by Erik Winfree March 1996
new non-parallel update, + simulated annealing & energy calculations
3/25/96 changed S.A. model so that in empty space, emptiness is
overwhelmingly favored.  added random number seed control.
11/8/98 split into two modules: xgrow.c and grow.c
simulation completely re-vamped for hierarchical kinetic model
11/27/98 added "error" colors... hydrolysis sim.
1/29/99  added output data
10/28/01 Xwindows didn't work on new 16-bit displays...
incorporating a few minor mods from Cruetz's latest&greatest
But, some problems with fission fragments not disappearing?
7/5/02   Michael DeLorimier adds general tile definitions
7/6/02   EW generalizes hydrolysis rules for any tile set
7/12/02  EW adds GUI for changing Gse and Gmc during simulation
7/16/02  EW adds irreversible Tile Assembly Model
7/17/02  EW adds default param defs to be read from tile file
(can be overwritten by command line)
Allow more field sizes, not just 256, and larger blocks.
(Note: crashes if you try too large for your display)
Allowed variable display rates.
Optimized choose_cell() and update_rates() using gprof.
Put estimated [DX] and temp on display
(from thesis pg 63, "Simulations" (1998) pg 10)
7/19/02  Major revision for multiple flakes.
Created "tube" structure to keep sim params.
5/2/03   minor enhancements by EW:
comments in tile file are allowed
5/3/03   added buttons for producing a "sample" AFM deposition image,
drawn from the distribution of flakes being simulated
and for exporting array information for MATLAB reading.
added 'fission_allowed' option, so as to disable
fission-causing dissociations.  this is now the default!
when considering flake_fission:
if 'wander' is set, and the cell chosen to
dissociate is the seed, then move the seed to its neighbor.
although 'wander' is the "right" model, it sometimes looks
strange, so I won't make it the default (yet).
nice little boxes drawn around blocks, if blocks>4
simulation no longer stops when event counter wraps around.
5/11/03  changed behavior for flakes with []'s:
if the flake is a single tile, then it "absorbs" no []
thus, an excess of single-tile flakes won't decrease the
free monomer concentration experienced by the other flakes.
coding implications:
no monomer depletion when initializing single-tile seed flakes.
when dimer flakes lose a tile, *both* concentrations are restored.
when monomer flake gains a tile, *both* concentrations are depleted.
changed seed tile behavior for rates:
event for removing seed tile is always present (unless tiles==1)
but it may have no effect (except to increment time)
if seed can't wander.
if monomer flake can't grow (because flake_conc > all seed tile conc)
then seed tile type changes to a random value
and this counts as a "hydrolysis" event, oddly enough
sampling now ignores monomer tiles if possible
flake G now reports G_bonds + G_concs *including* the seed tile
(if wander is off, then screen display offsets for seed tile [])
actually fixed the simulation event counter wrap-around.
5/25/03  enhanced tile file definition to allow X named colors.
e.g., {2 3 5 1}[.4](purple)
summer 03 Added command-line option to read in initial tile assembly.  (Shaun Lee)
10/14/03 Changed format of input file to allow a generalized strength function of the form
g(a,b) = c where a and b are tile types, and c is the floating point strength.
(Rebecca Schulman)
   11/10/03 Adding chunk_fission option.  When set, not only can single tiles dissociate,
   but pairs of tiles and 2x2 blocks can dissociate.  Implies fission.
   Should work with periodic and wander options, also.
   Was a bit hairy, especially with boundary conditions.
   Probably a few bugs still lurk here -- I yearn for automated testing.
   chunk_fission is incompatible with hydrolysis rules.
(Erik Winfree)
   11/16/03 Fixed color use for Err & Hyd display.
   Fixed puncture so it fissions if you cut the flake in two.
   Added error messages when tile file read fails.
   Added middle mouse click -> tile # underneath pointer is identified.
   12/29/03
   Note: Mism()      == exists a neighbor for which both se are non-null and different.
   errortile() == exists a neighbor for which se are different.
   (thus it includes growth front and tiles abutting null bonds
    i.e. a "dangling" sticky end)
   Previously: CONNECTED is true for mismatched sticky ends.
   Only not connected if bond 0 is involved.
   But otherwise, zero-strength bonds connect.
   Without this convention, there would be no speckle of random tiles
   around the edge (which is nice, because it visually illustrates that
    any tile can be added anywhere) if all zero-strength interactions
   immediately & instantly dissociate, as they do for speed reasons.
   New:        CONNECTED is true iff there is a non-zero bond (regardless how small).
   Thus, unmatches sticky-ends may dangle in the breeze
   without considering the tiles to be connected.

   This, in combination with an additional quick test to make sure that removing
   a tile won't lead to fission, greatly speeds up the simulation in cases where
   previously there was a deadly slow-down: if strength-0 edges abutted a
   growing region (as happens when sierpinski boundary erroneously grows inside the
    flake) the dissociations are frequent, and required flake_fission fills to
   be performed frequently.
   locally_fission_proof() was added to augment the ring[] test (which verifies
    that all neighbors of a removed tile are connected to each other) with a slightly
   slower but still fast (compared to the full fission fill test) test that verifies
   that the local connectivity of neighbors is *unchanged* with and without the tile.

   12/30/03 Added fill_cycles and fill_X to fill-in small holes after simulation.
   Added zero_bonds option to allow accretion of tiles that "don't stick".
   1/4/04  Added fill / clean button to help clarify exactly what this processing
   does, so error rate stats can be interpretted.
(Also fixed Quit so that it exits normally, outputting data if requested.)
   1/6/04  Added error_radius option.  I am worried that major counted error rates
   are temporary growth front mishaps, when internal errors are very very small.
   Testing...
   1/9/04  Added display option to show tile sides in null/weak/strong colors (if block>4)
   null:   Gse_EW with EW neighbor is < .5 Gse            (e.g. mismatch or null neighbor)
   weak:   Gse_EW with EW neighbor is > .5 and < 1.5 Gse  (e.g. regular strength-1 bond)
strong: Gse_EW with EW neighbor is > 1.5 Gse           (e.g. strength-2 bond)
   1/18/04 Added Rx to clean/fill button, and repair_unique_T to command-line.
   This provides an improved (I think) mismatch measurement, as it is guaranteed to
   fill in in all interior holes (thus not undercounting), and it removes all mismatches
   and fills in with unique correct aTAM tiles if possible (thus not overcounting).
   But in cases of big erroneous holes, it can (a) result in many mismatches, since
non-unique filling proceeds somewhat arbitrarily to completion, and (b) result in a
non-connected flake, which is weird if the simulation is continued. (So far, non-connected
flakes always get cleaned up well, but I'd still call this a BUG.)   Also, since only
interior holes get filled in, any error that stalls growth or leaves an inlet will not get
counted. Also added command-line "pause", primarily for examining saved files w/o starting
simulator. 2/25/04 Fixed bug in locally_fission_proof() call when chunk_fission is set.
(Sometimes bogus value resulted in non-fission when fission should have occurred.)
   6/4/04  Bond types can now be given as names.  See spiral.tiles for example.
   6/20/04 Added command-line option for setting tile stoichiometry.
   Made sure file is flushed immediately after export.
   3/8/05  Fixed bug in chunk_fission that caused bogus "impossible chunk chosen" errors
(RS)

   3/8/05 Implemented a "doubletile" function in configuration files
   that allows two tile types to be considered a connected,
   horizontal rectangle for purposes of tile additions and removals.  This works as
follows:
   - Tile types are entered as usual.  Double tiles are submitted
   as two single tiles.  The bond between the tile types is ignored
   - The option "doubletile(4,5)" listed in the options section
   creates a double tile from tiles 4 and 5, with 4 to the left of 5.
   - Double tiles are added as follows:
   - the first tile can be added anywhere it is adjacent to
   the assembly, provided that the second tile will fit
   geometrically.
   - the second tile can be added only if the first tile is not
   adjacent to the current assembly
   - Double tiles are removed as follows:
   - Only the first of the two tiles has a non-zero off rate.
   - If it is chosen, the binding energy is the sum of the
   energies attaching both tiles to the assembly.
   If it is removed, both tiles are removed.
(RS)
   3/19/05 Adding blast damage for demonstrating self-repair.  alpha, beta, gamma
command-line args. Wander allows the original seed tile to be lost. It's a bit strange
when no_fission is used, because fission is tested sequentially as tiles in the square are
being removed, rather than after the whole square has been removed.

   5/14/05 Added an anneal feature, at which Gse starts at some value and
   declines by a time constant until it reaches its final value.  It shouldn't be hard to
implement a melt feature, with the exact reverse behavior, but right now we require that
Gse initial be higher than Gse final.  Note that it would be horribly inefficient to
implement an exact annealing feature of this type, because all rates would have to be
updated after every event.  As a compromise we update Gse (and therefore the rates) some
fixed number of (for now 100) times per time constant. -- RS 7/30/05 Added testing to
xgrow -- new option "testing" will test whether a tile set obeys detailed balance once it
reaches equilibrium.  Takes a while to do, and effectively forever if you choose a tile
set that takes a looooong time to reach equilibrium, like the Sierpinski set. -- RS
   11/20/05 Fixed a bug in count_flakes().  It's still very brittle.  -EW
   12/12/05 Fixed up the handling of Gse and tp->Gse, which was leading to the GUI display
being occasionally off. -EW 12/15/05 Added flake loading option to offset the seed
assembly position -EW

   .... various undocumented improvements by Erik Winfree (EW), Rebecca Schulman (RS), and
Constantine Evans (CE) ....

2008 Changed tile and bond data type to short int, for > 256 types.  For smaller memory
footprint, xgrow-small is avalailable. (CE)

   3/1/09 Fixed aTAM bug where it hangs if there's no possible move.  (EW & CE)

   TO DO List:

   * If the tile set specifies a stoichiometry of 0 (e.g. for the seed), the simulation
can freak out.
   * export ALL may actually be ALL-BUT-ONE.  Test and fix, please?
   * add a DILUTE/CONCENTRATE button.
   * files should be flushed after each use.  arrayfile should *append*, not overwrite,
properly adding to flake{} array at the appropriate index?
   * option for chunk_fission, but disallowing larger flake_fission events, would be
desirable, so as to have a simulation in detailed balance.
   * for the same reason, we'd need a chunk_fusion option -- calculates equilibrium
concentrations for all dimers & 2x2 blocks, and adds these as on-event options.  total
on-rate would include these; in circumstances where they can't occur, result would just be
a non-event.  There are so many 2x2 blocks, we'd need a binary tree to efficiently select
one. (Not a bad idea for large tile sets, even w/o chunks.)  To avoid non-events, we could
have 4 "empty" arrays, one for each size.  This would be a significant change to the code.
* Mism() and errortile()
   ought to be redefined to be sensible for g() type glues implementing complementary
   (rather than self-complementary) sticky ends.  E.g., precalculate bestmatch[bond_type]
   * tile files should allow sticky ends to be specified by string labels rather
   than numbers.  num_tiles and num_bindings should be ignored, and inferred from
the file.  (compatibility with existing xgrow tile files must be maintained.)
   * Would like auto-rotate for tile files, to simulate DAO-E flip, DAE-O symmetry, etc.
   * Is it possible to reverse the random-number generator, so as to "reverse time"?
   That way, you see something interesting, and you can ask "How did that happen again?"
   In that case, it would be useful to have GUI to change update_rate.
   * Would be nice to change other parameters in addition to Gse, Gmc.  E.g. inc/dec by
   left/right mouse click on the numbers in the display.
   * Something like "multiflakes=100@27" argument does "the right thing"
   by adding 100 flakes for each tile type, at Gfc=27+stoich,
   with each seed centered in the field,
   defaults to "wander"
   * Event counter should simply have more bits!  Currently, it wraps around...
   * In no-fission mode, one can get caught up on very fast -- but disallowed --
   dissociations, which must be rejected.  This is not good.
(not sure if this is still true.  EW 11/10/03  I think it is fixed. EW 1/9/04)

   * rubberbanding for "puncture"
   * green dot (current value), red dot (selection cursor) for Gse/Gmc mouse choice
   * fp->G nan and other discrepencies should be tracked down.
   * Very Much want to implement a Surface mode, with Gsf for the surface, and
   allowing disconnected flakes.

   * importfile has problems.
   * for blast, need option that ignores seed
   * for fission, need option that chooses largest component, not component connected to
the seed
   * program should die when user closes the window, rather than only when the user
presses "quit"
   * when requested size is larger than the screen, the program should quit rather than
   crash.  (See the code near "BadMatch", search for it below.)

   *  The (undocumented) hydrolysis mode is now broken; it crashes.  This model is a 2D
generalization of microtubule / actin hydrolysis that leads to treadmilling and dynamic
instability; in 1999 I was curious if it could be utilized to decrease error rates.
   *  Xgrow still sometimes hangs (or doesn't allow you to quit) when there is no
available move, but it's not obvious when this happens.  Example: "xgrow barseed T=2
block=2 seed=180,180,1"  Here, tile 5 should be the seed.
   *  With multiple flakes, if the simulation is running, only the last flake is shown. To
see other flakes, you have to pause and then use next/big/prev to find one of interest.
But then when you run the simulation again, it pops to the last flake.  I think one would
like to keep looking at the same flake you were looking at before pressing RUN.  Example:
"xgrow barcode block=2 size=128 Gmc=15 Gse=10 addflakes=130,130,1:100@24"


   Compiling:  run "make"

   */

#include <X11/Xatom.h>
#include <X11/Xlib.h>
#include <X11/Xos.h>
#include <X11/Xutil.h>
#include <assert.h>
#include <libgen.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "grow.h"
#ifdef TESTING_OK
#include "xgrow-tests.h"
#endif

/* lattice dimensions (plus two for boundaries): */
/* NCOLS should be a multiple of bytes per long word (BPW) */
/* THIS is for X bitmap stuff, so add total 4 for boundary */
#define NBDY 2

/* The max number of info lines printed at the top of the window,
   for spacing purposes (will be multiplied by font_height). */
#define NUM_INFO_LINES 6

/* Aesthetic margin of image area. */
#define PLAY_MARGIN 10

#define IS_ARG_MATCH(arg, s) (strncmp(arg, s, strlen(s)) == 0)

long int translate[MAXTILETYPES]; /* for converting colors */
int paused = 0, errorc = 0, errors = 0, sampling = 0;
int export_mode = 0, export_flake_n = 1, export_movie_n = 1, export_movie = 0;
FILE *export_fp = NULL;
static char *progname;

#define BUFFER_LEN 256

char stringbuffer[BUFFER_LEN];
char tileset_name[BUFFER_LEN];
int testing = 0;
int initial_rc = 1;

/* various window stuff */
Display *display;
int screen;
Window window, quitbutton, pausebutton, playground, restartbutton, colorbutton,
    sidebutton, flakebutton, seedbutton, tempbutton, fissionbutton, samplebutton,
    exportbutton, cleanbutton;
GC gc, gcr, gccolor;
XImage *spinimage = NULL;
XFontStruct *font = NULL;
int font_height;
XSizeHints size_hints;
long event_mask;
int depth;
int darkcolor, lightcolor, black, white, errorcolor, goodcolor, hydrocolor, hydroerrcolor,
    strongcolor, weakcolor, nullcolor;
char *tile_colors[MAXTILETYPES] = {"black",     "blue",     "red",        "green",
                                   "yellow",    "gold",     "purple",     "white",
                                   "dark blue", "dark red", "dark green", "wheat",
                                   "orange",    "cyan",     "light grey"};

/* simulation parameters */
// tube *tp; flake *fp;

int first_tile, second_tile;

int NROWS, NCOLS, VOLUME, WINDOWWIDTH, WINDOWHEIGHT;
int block = 1; /* default to small blocks; calling with argument changes this */
int linear;
FILE *tracefp, *datafp, *arrayfp, *tilefp, *largeflakefp;
FILE *untiltilescountfp;
int tthresh = 0;
int N, num_bindings, tileb_length;
int clean_cycles = 0;
double clean_X = 1.0;
int fill_cycles = 0;
double fill_X = 1.0;
double error_radius = 0.0;
double repair_unique_T = 2.0;
int repair_unique = 0;
double tmax;
int emax, smax, fsmax, smin, mmax;
char *stripe_args = NULL;
int XXX = 1; /* If 1, draw the simulation, otherwise don't.  (as inferred from the effect
                of -nw)*/
int import = 0; /* Are we importing flakes? */
int import_flake_size = 0;
int import_offset_i = 0, import_offset_j = 0;
char font_selection_string[256] = "9x13";

const char NEWLINE[2] = "\n\0";

struct flake_param {
   int seed_i, seed_j, seed_n, N;
   double Gfc;
   FILE *import_from;
   struct flake_param *next_param;
};
struct flake_param *fparam = NULL, *fprm;

int count_flakes(FILE *flake_file);

typedef struct tube_params_struct {
   // formerly shared with grow
   double blast_rate_alpha; //  =0;
   double blast_rate_beta;  //  =4;  // k>3 required for finite rate of blasting a given
                            //  tile in infinite size flakes (gamma=0)
   double blast_rate_gamma; //  =0;
   double blast_rate;       //  =0;
   double min_strength;     //  =1;
   int wander;              //  ,
   int periodic, linear, fission_allowed, zero_bonds_allowed;
   int *present_list;    //  =NULL;
   int present_list_len; //  =0;
   int untiltiles;       //  =0,
   int untiltilescount;  //=0;
   int **tileb;
   double *strength;
   double **glue;
   double *stoic;
   double anneal_g;
   double anneal_t;
   int updates_per_RC;
   double anneal_h;
   double anneal_s;
   double startC;
   double endC;
   double seconds_per_C;
   int double_tiles;
   int *dt_right;
   int *dt_left;
   int *dt_down;
   int *dt_up;
   int hydro;
   double k;
   double Gmc;
   double Gse;
   double Gmch;
   double Gseh;
   double Ghyd;
   double Gas;
   double Gam;
   double Gae;
   double Gah;
   double Gao;
   double T;
   double tinybox;
   int seed_i;
   int seed_j;
   int seed_n;
   double Gfc;
   char *tile_names[MAXTILETYPES];
   int size;
   int size_P;
   int update_rate;
} tube_params;

void set_default_params(tube_params *params) {
   params->blast_rate_alpha = 0;
   params->blast_rate_beta = 4; // k>3 required for finite rate of blasting a given tile
                                // in infinite size flakes (gamma=0)
   params->blast_rate_gamma = 0;
   params->blast_rate = 0;
   params->min_strength = 1;
   params->wander = 0;
   params->periodic = 0;
   params->linear = 0;
   params->fission_allowed = 0;
   params->zero_bonds_allowed = 0;
   params->present_list = NULL;
   params->present_list_len = 0;
   params->untiltiles = 0;
   params->untiltilescount = 0;
   params->double_tiles = 0;
   params->Gmc = 17;
   params->Gse = 8.6;
   params->k = 1000000.0;
   params->T = 0;
   params->Gmch = 30;
   params->Gseh = 0;
   params->Ghyd = 30;
   params->Gas = 30;
   params->Gam = 15;
   params->Gae = 30;
   params->Gah = 30;
   params->Gao = 10;
   params->seed_i = 250;
   params->seed_j = 250;
   params->seed_n = 1;
   params->hydro = 0;
   params->Gfc = 0;
   params->size = 256;
   params->size_P = 8;
   params->update_rate = 10000;
}

int parse_arg_line(char *arg, tube_params *params) {
   if (IS_ARG_MATCH(arg, "block="))
      block = MAX(1, MIN(30, atoi(&arg[6])));
   else if (IS_ARG_MATCH(arg, "size="))
      params->size = MAX(32, MIN(4096, atoi(&arg[5])));
   else if (IS_ARG_MATCH(arg, "rand=")) {
      srand48(atoi(&arg[5]));
      srandom(atoi(&arg[5]));
   } else if (IS_ARG_MATCH(arg, "k="))
      params->k = atof(&arg[2]);
   else if (IS_ARG_MATCH(arg, "Gmc="))
      params->Gmc = atof(&arg[4]);
   else if (IS_ARG_MATCH(arg, "Gse="))
      params->Gse = atof(&arg[4]);
   // FIXME FIXME FIXME: HYDRO OPTIONS ARE DISABLED
   // else if (IS_ARG_MATCH(arg,"Gmch=")) {hydro=1; Gmch=atof(&arg[5]);}
   // else if (IS_ARG_MATCH(arg,"Gseh=")) {hydro=1; Gseh=atof(&arg[5]);}
   // else if (IS_ARG_MATCH(arg,"Ghyd=")) {hydro=1; Ghyd=atof(&arg[5]);}
   else if (IS_ARG_MATCH(arg, "Gas=")) {
      params->hydro = 1;
      params->Gas = atof(&arg[4]);
   } else if (IS_ARG_MATCH(arg, "Gam=")) {
      params->hydro = 1;
      params->Gam = atof(&arg[4]);
   } else if (IS_ARG_MATCH(arg, "Gae=")) {
      params->hydro = 1;
      params->Gae = atof(&arg[4]);
   }
   // else if (IS_ARG_MATCH(arg,"Gah=")) {hydro=1; Gah=atof(&arg[4]);}
   // else if (IS_ARG_MATCH(arg,"Gao=")) {hydro=1; Gao=atof(&arg[4]);}
   else if (IS_ARG_MATCH(arg, "Gfc=")) {
      params->Gfc = atof(&arg[4]);
   } else if (IS_ARG_MATCH(arg, "stoic=")) {
      double sts;
      int stn;
      char *sp;
      sp = strchr(&arg[6], '[');
      if (sp != NULL) {
         stn = atoi(&arg[6]);
         sts = atof(sp + 1);
         if (stn <= N && stn > 0 && sts > 0) {
            params->stoic[stn] = sts;
            fprintf(stderr, "stoic of tile %d = [%f]\n", stn, sts);
         } else {
            if (sts == 0) {
               fprintf(stderr,
                       "Stoic == 0.0 is not allowed; use a small but positive value.\n");
               return -1;
            }
            fprintf(stderr, "Requested stoic of tile %d = [%f] is impossible!\n", stn,
                    sts);
            return -1;
         }
      }
   } else if (IS_ARG_MATCH(arg, "T="))
      params->T = atof(&arg[2]);
   else if (IS_ARG_MATCH(arg, "pause=True"))
      paused = 1;
   else if (IS_ARG_MATCH(arg, "pause=False"))
      paused = 0;
   else if (IS_ARG_MATCH(arg, "pause"))
      paused = 1;
   else if (IS_ARG_MATCH(arg, "periodic=True"))
      params->periodic = 1;
   else if (IS_ARG_MATCH(arg, "periodic=False"))
      params->periodic = 0;
   else if (IS_ARG_MATCH(arg, "periodic"))
      params->periodic = !params->periodic;
   else if (IS_ARG_MATCH(arg, "wander=True"))
      params->wander = 1;
   else if (IS_ARG_MATCH(arg, "wander=False"))
      params->wander = 0;
   else if (IS_ARG_MATCH(arg, "wander"))
      params->wander = !params->wander;
   else if (IS_ARG_MATCH(arg, "fission=off"))
      params->fission_allowed = 0;
   else if (IS_ARG_MATCH(arg, "fission=on"))
      params->fission_allowed = 1;
   else if (IS_ARG_MATCH(arg, "fission=chunk"))
      params->fission_allowed = 2;
   else if (IS_ARG_MATCH(arg, "no_fission"))
      params->fission_allowed = 0;
   else if (IS_ARG_MATCH(arg, "fission"))
      params->fission_allowed = 1;
   else if (IS_ARG_MATCH(arg, "chunk_fission"))
      params->fission_allowed = 2;
   else if (IS_ARG_MATCH(arg, "blast_rate_alpha="))
      params->blast_rate_alpha = atof(&arg[17]);
   else if (IS_ARG_MATCH(arg, "blast_rate_beta="))
      params->blast_rate_beta = atof(&arg[16]);
   else if (IS_ARG_MATCH(arg, "blast_rate_gamma="))
      params->blast_rate_gamma = atof(&arg[17]);
   else if (IS_ARG_MATCH(arg, "zero_bonds"))
      params->zero_bonds_allowed = 1;
   else if (IS_ARG_MATCH(arg, "no_zero_bonds"))
      params->zero_bonds_allowed = 0;
   else if (IS_ARG_MATCH(arg, "seed=")) {
      char *p = (&arg[5]);
      params->seed_i = atoi(p);
      if ((p = strchr(p, ',')) != NULL) {
         params->seed_j = atoi(p + 1);
         if ((p = strchr(p + 1, ',')) != NULL)
            params->seed_n = atoi(p + 1);
      }
   } else if (IS_ARG_MATCH(arg, "addflakes=")) {
      char *p = (&arg[10]);
      fprm = (struct flake_param *)malloc(sizeof(struct flake_param));
      fprm->seed_i = fprm->seed_j = 130;
      fprm->seed_n = 1;
      fprm->N = 1;
      fprm->Gfc = 0;
      fprm->import_from = NULL;
      fprm->next_param = fparam;
      fparam = fprm;
      fparam->seed_i = atoi(p);
      if ((p = strchr(p, ',')) != NULL) {
         fparam->seed_j = atoi(p + 1);
         if ((p = strchr(p + 1, ',')) != NULL) {
            fparam->seed_n = atoi(p + 1);
            if ((p = strchr(p + 1, ':')) != NULL) {
               fparam->N = atoi(p + 1);
               if ((p = strchr(p + 1, '@')) != NULL) {
                  fparam->Gfc = atof(p + 1);
               }
            }
         }
      }
   } else if (IS_ARG_MATCH(arg, "tinybox=")) {
      char *p = (&arg[8]);
      params->tinybox = atof(p);
   } else if (IS_ARG_MATCH(arg, "stripe=")) {
      stripe_args = (&arg[7]);
      params->periodic = 1;
      params->wander = 1;
   } else if (strcmp(arg, "-nw") == 0)
      XXX = 0;
   else if (IS_ARG_MATCH(arg, "window=True"))
      XXX = 1;
   else if (IS_ARG_MATCH(arg, "window=False"))
      XXX = 0;
   else if (strcmp(arg, "-linear") == 0)
      linear = 1;
   else if (strncmp(arg, "linan=", 6) == 0) {
      char *c;
      c = &arg[6];
      params->anneal_h = atof(c);
      c = strchr(c, ',') + 1;
      if (c == NULL) {
         fprintf(stderr, "Usage : linan=h,s,C0,Cfin,dt ");
         return -1;
      }
      params->anneal_s = atof(c);
      c = strchr(c, ',') + 1;
      if (c == NULL) {
         fprintf(stderr, "Usage : linan=h,s,C0,Cfin,dt ");
         return -1;
      }
      params->startC = atof(c);
      c = strchr(c, ',') + 1;
      if (c == NULL) {
         fprintf(stderr, "Usage : linan=h,s,C0,Cfin,dt ");
         return -1;
      }
      params->endC = atof(c);
      c = strchr(c, ',') + 1;
      if (c == NULL) {
         fprintf(stderr, "Usage : linan=h,s,C0,Cfin,dt ");
         return -1;
      }
      params->seconds_per_C = atof(c);
   } else if (strncmp(arg, "anneal=", 7) == 0) {
      char *c;
      c = &arg[7];
      params->anneal_g = atof(c);
      c = strchr(c, ',');
      if (c == NULL) {
         fprintf(stderr, "Anneal option requires both an initial Gse and a time constant "
                         "of annealing.\n");
         return -1;
      }
      params->anneal_t = atof(c + 1);
      if (params->anneal_t <= 0) {
         fprintf(stderr, "Time constant for annealing must be positive.\n");
         return -1;
      }
   } else if (IS_ARG_MATCH(arg, "doubletile=")) {
      char *p;
      first_tile = atoi(&arg[11]);
      p = strchr(&arg[12], ',');
      if (p == NULL) {
         fprintf(stderr, "Two tiles with a comma separating them are required for the "
                         "double tile option.\n");
         return -1;
      }
      second_tile = atoi(p + 1);
      if (params->tileb[first_tile][1] != params->tileb[second_tile][3]) {
         fprintf(stderr,
                 "Pieces of a double tile must have a matching bond in between them.\n");
         return -1;
      }
      params->dt_right[first_tile] = second_tile;
      params->dt_left[second_tile] = first_tile;
      params->double_tiles = 1;
   } else if (IS_ARG_MATCH(arg, "vdoubletile=")) {
      char *p;
      first_tile = atoi(&arg[12]);
      p = strchr(&arg[13], ',');
      if (p == NULL) {
         fprintf(stderr, "Two tiles with a comma separating them are required for the "
                         "vdouble tile option.\n");
         return -1;
      }
      second_tile = atoi(p + 1);
      if (params->tileb[first_tile][2] != params->tileb[second_tile][0]) {
         fprintf(stderr,
                 "Pieces of a vdouble tile must have a matching bond in between them.\n");
         return -1;
      }
      params->dt_down[first_tile] = second_tile;
      params->dt_up[second_tile] = first_tile;
      params->double_tiles = 1;
   } else if (IS_ARG_MATCH(arg, "update_rate="))
      params->update_rate = MAX(1, MIN(atol(&arg[12]), 10000000));
   else if (IS_ARG_MATCH(arg, "tracefile="))
      tracefp = fopen(strtok(&arg[10], NEWLINE), "a");
   else if (IS_ARG_MATCH(arg, "movie")) {
      export_mode = 2;
      export_movie = 1;
   } else if (IS_ARG_MATCH(arg, "tmax="))
      tmax = atof(&arg[5]);
   else if (IS_ARG_MATCH(arg, "emax="))
      emax = atoi(&arg[5]);
   else if (IS_ARG_MATCH(arg, "smax="))
      smax = atoi(&arg[5]);
   else if (IS_ARG_MATCH(arg, "mmax="))
      mmax = atoi(&arg[6]);
   else if (IS_ARG_MATCH(arg, "smin="))
      smin = atoi(&arg[5]);
   else if (IS_ARG_MATCH(arg, "fsmax="))
      fsmax = atoi(&arg[6]);
   else if (IS_ARG_MATCH(arg, "untiltiles=")) {
      int i = 0;
      char *pos;
      pos = arg;
      while (pos != NULL) {
         pos = strchr(pos + 1, ',');
         params->present_list_len++;
      }
      params->present_list = (int *)malloc(params->present_list_len * sizeof(int));
      pos = &arg[11];
      while ((pos - 1) != NULL) {
         params->present_list[i++] = atoi(pos);
         pos = strchr(pos, ',') + 1;
      }
      params->untiltiles = 1;
   } else if (IS_ARG_MATCH(arg, "untiltilescount=")) {
      int i = 0;
      char *pos;
      params->untiltilescount = 1;
      pos = arg;
      while (pos != NULL) {
         pos = strchr(pos + 1, ',');
         params->present_list_len++;
      }
      params->present_list = (int *)malloc(params->present_list_len * sizeof(int));
      pos = &arg[16];
      while ((pos - 1) != NULL) {
         params->present_list[i++] = atoi(pos);
         pos = strchr(pos, ',') + 1;
      }
      params->untiltilescount = 1;
   } else if (IS_ARG_MATCH(arg, "clean_cycles="))
      clean_cycles = atoi(&arg[13]);
   else if (IS_ARG_MATCH(arg, "clean_X="))
      clean_X = atof(&arg[8]);
   else if (IS_ARG_MATCH(arg, "fill_cycles="))
      fill_cycles = atoi(&arg[12]);
   else if (IS_ARG_MATCH(arg, "fill_X="))
      fill_X = atof(&arg[7]);
   else if (IS_ARG_MATCH(arg, "error_radius="))
      error_radius = atof(&arg[13]);
   else if (IS_ARG_MATCH(arg, "repair_unique_T=")) {
      repair_unique_T = atof(&arg[15]);
      repair_unique = 1;
   } else if (IS_ARG_MATCH(arg, "datafile="))
      datafp = fopen(strtok(&arg[9], NEWLINE), "a");
   else if (IS_ARG_MATCH(arg, "largeflakedatafile=")) {
      char *c;
      c = strchr(arg, ',') + 1;
      tthresh = atoi(&arg[19]);
      largeflakefp = fopen(c, "a");
   } else if (IS_ARG_MATCH(arg, "untiltilescountfile=")) {
      untiltilescountfp = fopen(&arg[20], "a");
   } else if (IS_ARG_MATCH(arg, "arrayfile="))
      arrayfp = fopen(strtok(&arg[10], NEWLINE), "w");
   else if (IS_ARG_MATCH(arg, "exportfile="))
      export_fp = fopen(strtok(&arg[11], NEWLINE), "w");
   else if (strncmp(arg, "testing", 7) == 0) {
      testing = 1;
   } else if (IS_ARG_MATCH(arg, "import_offset=")) {
      char *p = (&arg[14]);
      import_offset_i = atoi(p);
      if ((p = strchr(p, ',')) != NULL) {
         import_offset_j = atoi(p + 1);
      }
   } else if (IS_ARG_MATCH(arg, "importfile")) {
      char *p = (&arg[11]);
      FILE *import_fp;
      import = 1;
      fprm = (struct flake_param *)malloc(sizeof(struct flake_param));
      fprm->seed_i = fprm->seed_j = 130;
      fprm->seed_n = 1;
      fprm->N = 1;
      fprm->Gfc = 0;
      fprm->import_from = NULL;
      fprm->next_param = fparam;
      fparam = fprm;
      /*
       * Need to figure out seed_i and seed_j after loading in the cells.
       * Randomly choose spots from the smallest flake until we find a tile
       * that is present in each flake.
       */
      fparam->seed_i = fparam->seed_j = params->size / 2;
      fparam->seed_n = 1;

      if (IS_ARG_MATCH(arg, "importfile=")) {
         char imp_fn[512], *arg_fn;
         import_fp = fopen(arg_fn = strtok(p, NEWLINE), "r");
         if (import_fp == NULL) {
            snprintf(&imp_fn[0], 512, "%s", arg_fn);
            import_fp = fopen(&imp_fn[0], "r");
         }
         if (import_fp == NULL) {
            snprintf(&imp_fn[0], 512, "%s.seed", arg_fn);
            import_fp = fopen(&imp_fn[0], "r");
         }
         char *tileset_dir;
         tileset_dir = dirname(tileset_name);
         if (import_fp == NULL) {
            snprintf(&imp_fn[0], 512, "%s/%s", tileset_dir, arg_fn);
            import_fp = fopen(&imp_fn[0], "r");
         }
         if (import_fp == NULL) {
            snprintf(&imp_fn[0], 512, "%s/%s.seed", tileset_dir, arg_fn);
            import_fp = fopen(&imp_fn[0], "r");
         }
      } else {
         import_fp = fopen("xgrow_export_output", "r");
      }
      /* If the file does not exist. */
      if (import_fp == NULL) {
         fprintf(stderr, "Error: Import file not found.\n");
         return -1;
      }
      fparam->import_from = import_fp;
      if ((p = strtok(NULL, NEWLINE)) != NULL)
         fparam->Gfc = atof(p);
      else
         fparam->Gfc = 0;
      fparam->N = count_flakes(import_fp);
   } else if (IS_ARG_MATCH(arg, "min_strength=")) {
      params->min_strength = atof(&arg[13]);
   } else if (IS_ARG_MATCH(arg, "font=")) {
      strcpy(font_selection_string, strtok(&arg[5], "\n\0"));
   } else {
      fprintf(stderr, "Could not parse argument '%s'\n", arg);
      return -1;
   }
   return 0;
}

#define rsc read_skip_comment(tilefp)
void read_skip_comment(FILE *fp)
// anything after a "%" gets ignored
{
   int c;
   fscanf(fp, " ");
   c = fgetc(fp);
   while (c == '%') {
      fgets(&stringbuffer[0], 256, fp);
      // fprintf(stderr,"%%%s",stringbuffer);
      fscanf(fp, " ");
      c = fgetc(fp);
   }
   ungetc(c, fp);
}

void read_tilefile(FILE *tilefp, tube_params *params) {
   float strength_float, glue_float, stoic_float;
   int i, j, k;
   int temp_char;
   char s[255], **btnames;
   int n, m, r;
   int return_code;

   rewind(tilefp);
   // needs to be read twice,
   //  once to get sim specs, field size, etc
   //  once to get colors after X window is open
   rsc;
   fscanf(tilefp, "tile edges matches {{N E S W}*}\n");
   rsc;

   if (1 != fscanf(tilefp, "num tile types=%d\n", &N)) {
      fprintf(stderr, "Reading tile file: expected num tile types.\n");
      exit(-1);
   }
   rsc;
   if (1 != fscanf(tilefp, "num binding types=%d\n", &num_bindings)) {
      fprintf(stderr, "Reading tile file: expected num binding types.\n");
      exit(-1);
   }
   rsc;

   if (N > MAXTILETYPES) {
      fprintf(stderr, "Reading tile file: too many tile types.\n");
      exit(-1);
   }
   rsc;
   if (num_bindings > MAXTILETYPES) {
      fprintf(stderr, "Reading tile file: too many binding types.\n");
      exit(-1);
   }
   rsc;

   btnames = (char **)malloc((num_bindings + 1) * sizeof(char *));
   for (k = 0; k <= num_bindings; k++)
      btnames[k] = "null"; // until overwritten by tile file specification;

   temp_char = getc(tilefp);
   ungetc(temp_char, tilefp);
   if (temp_char == 'b') {
      r = 0;
      fscanf(tilefp, "binding type names=%n", &r);
      rsc;
      if (r != 19) {
         fprintf(stderr, "Reading tile file: expected binding type names.\n");
         exit(-1);
      }
      r = 0;
      fscanf(tilefp, "{%n", &r);
      rsc;
      if (r != 1) {
         fprintf(stderr, "Reading tile file: expected binding type names {.\n");
         exit(-1);
      }

      for (i = 1; i <= num_bindings; i++) {
         if (1 != fscanf(tilefp, "%100s", &s[0])) {
            fprintf(stderr, "Reading tile file: expected binding type %d name.\n", i);
            exit(-1);
         }
         if (s[strlen(s) - 1] == '}') {
            ungetc('}', tilefp);
            s[strlen(s) - 1] = 0;
         }
         if (strlen(s) == 0) {
            fprintf(stderr, "Reading tile file: expected binding type %d name.\n", i);
            exit(-1);
         }
         btnames[i] = strdup(s);
         rsc;
         if (index("0123456789", btnames[i][0]) != NULL) {
            fprintf(
                stderr,
                "Reading tile file: binding type %d name ('%s') cannot begin with 0-9.\n",
                i, s);
            exit(-1);
         }
      }
      rsc;
      r = 0;
      fscanf(tilefp, "}%n", &r);
      rsc;
      if (r != 1) {
         fprintf(stderr, "Reading tile file: expected binding type names }.\n");
         exit(-1);
      }
   }

   r = 0;
   fscanf(tilefp, "tile edges=%n", &r);
   rsc; // fprintf(stderr, "r=%d\n",r);
   if (r != 11) {
      fprintf(stderr, "Reading tile file: expected `tile edges=' declaration.\n");
      exit(-1);
   }
   tileb_length = N + 1;
   params->tileb = (int **)calloc(sizeof(int *), tileb_length);
   params->stoic = (double *)calloc(sizeof(double), tileb_length);
   r = 0;
   fscanf(tilefp, "{%n", &r);
   rsc;
   if (r != 1) {
      fprintf(stderr, "Reading tile file: expected tile set start {.\n");
      exit(-1);
   }
   params->tileb[0] = (int *)calloc(sizeof(int), 4);
   for (j = 0; j < 4; j++) {
      params->tileb[0][j] = 0;
   }
   params->stoic[0] = 0;
   for (i = 1; i < tileb_length; i++) {
      params->tileb[i] = (int *)calloc(sizeof(int), 4);
      r = 0;
      fscanf(tilefp, "{%n", &r);
      rsc;
      if (r != 1) {
         fprintf(stderr, "Reading tile file: expected tile %d start def {. \n", i);
         exit(-1);
      }
      /* read in the four binding types {N E S W} */
      for (j = 0; j < 4; j++) {
         rsc;
         temp_char = getc(tilefp);
         ungetc(temp_char, tilefp);
         if (index("0123456789", temp_char) != NULL) {
            if (1 != fscanf(tilefp, "%d", &params->tileb[i][j])) {
               fprintf(stderr, "Reading tile file: expected tile %d bond %d value.\n", i,
                       j);
               exit(-1);
            }
         } else {
            if (1 != fscanf(tilefp, "%[^} \t]", &s[0])) {
               fprintf(stderr, "Reading tile file: expected tile %d bond %d's name.\n", i,
                       j);
               exit(-1);
            }
            for (k = 0; k <= num_bindings; k++)
               if (strcmp(s, btnames[k]) == 0)
                  break;
            if (k <= num_bindings) {
               params->tileb[i][j] = k;
            } else {
               fprintf(
                   stderr,
                   "Reading tile file: expected tile %d bond %d's name '%s' unknown.\n",
                   i, j, s);
               exit(-1);
            }
         }
      }
      rsc;
      r = 0;
      fscanf(tilefp, "}%n", &r);
      rsc;
      if (r != 1) {
         fprintf(stderr, "Reading tile file: expected tile %d end def }. \n", i);
         exit(-1);
      }
      if (fscanf(tilefp, "[%g]", &stoic_float))
         params->stoic[i] = stoic_float;
      else
         params->stoic[i] = 1.0;
      rsc;
      if (fscanf(tilefp, "(%200[^)])", &stringbuffer[0])) {
         tile_colors[i] = (char *)malloc(strlen(stringbuffer) + 2);
         strcpy(tile_colors[i], stringbuffer);
      }
      if (fscanf(tilefp, "<%200[^>]>", &stringbuffer[0])) {
         params->tile_names[i] = (char *)malloc(strlen(stringbuffer) + 2);
         strcpy(params->tile_names[i], stringbuffer);
      }
      fscanf(tilefp, "\n");
      rsc;
   }
   r = 0;
   fscanf(tilefp, "}%n", &r);
   rsc;
   if (r != 1) {
      fprintf(stderr, "Reading tile file: expected tile set end }.\n");
      exit(-1);
   }
   // fprintf(stderr, "Tile set loaded (%d tiles)\n",N);
   // for (i=1;i<tileb_length;i++) {
   //   for (j=0;j<4;j++) fprintf(stderr, "%d ",tileb[i][j]); fprintf(stderr, "\n");
   // }

   params->glue = (double **)calloc(sizeof(double *), num_bindings + 1);
   for (i = 0; i <= num_bindings; i++) {
      int j;
      params->glue[i] = (double *)calloc(sizeof(double), num_bindings + 1);
      for (j = 0; j <= num_bindings; j++) {
         params->glue[i][j] = 0;
      } // necessary??  -- EW
   }
   params->strength = (double *)calloc(sizeof(double), num_bindings + 1);
   params->dt_left = (int *)calloc(sizeof(int), N + 1);
   params->dt_right = (int *)calloc(sizeof(int), N + 1);
   params->dt_up = (int *)calloc(sizeof(int), N + 1);
   params->dt_down = (int *)calloc(sizeof(int), N + 1);

   temp_char = getc(tilefp);
   ungetc(temp_char, tilefp);
   if (temp_char == 'b') {
      r = 0;
      fscanf(tilefp, "binding strengths=%n", &r);
      rsc;
      if (r != 18) {
         fprintf(stderr, "Reading tile file: expected binding strength defs.\n");
         exit(-1);
      }
      r = 0;
      fscanf(tilefp, "{%n", &r);
      rsc;
      if (r != 1) {
         fprintf(stderr, "Reading tile file: expected binding strength defs {.\n");
         exit(-1);
      }
      params->strength[0] = 0; /* bond type 0 ("null") always has strength 0 */
      for (i = 1; i <= num_bindings; i++) {
         if (1 != fscanf(tilefp, "%g", &strength_float)) {
            fprintf(stderr, "Reading tile file: expected binding strength %d value.\n",
                    i);
            exit(-1);
         }
         params->strength[i] = (double)strength_float;
         // fprintf(stderr, "strength for se #%d = %g\n",i,strength[i]);
      }
      rsc;
      r = 0;
      fscanf(tilefp, "}%n", &r);
      rsc;
      if (r != 1) {
         fprintf(stderr, "Reading tile file: expected binding strength defs }.\n");
         exit(-1);
      }
   }
   while ((temp_char = fgetc(tilefp)) == 'g') {
      if (3 != fscanf(tilefp, "(%d,%d)=%g\n", &n, &m, &glue_float)) {
         fprintf(stderr, "Reading tile file: expected glue def.\n");
         exit(-1);
      }
      rsc;
      // printf ("Glue float is g(%d,%d)=%g\n",n,m,glue_float);
      params->glue[n][m] = (double)glue_float;
      params->glue[m][n] = (double)glue_float;
   }
   ungetc(temp_char, tilefp);
   rsc;
   // fprintf(stderr, "Bond strengths loaded (%d bond types)\n",num_bindings);
   // for (i=1;i<=num_bindings;i++) fprintf(stderr, "%f ",strength[i]); fprintf(stderr,
   // "\n");

   rsc;
   while (fgets(&stringbuffer[0], 256, tilefp) != NULL) {
      return_code = parse_arg_line(&stringbuffer[0], params);
      if (return_code) {
         exit(-1);
      }
      rsc;
   }

   fclose(tilefp);

   // better not free this: some strings are constants, some are alloc'ed.
   // for (i=1; i<=num_bindings;i++) free(btnames[i]); free(btnames);
}

void getargs(int argc, char **argv, tube_params *params) {
   int i;
   struct flake_param *fprm;
   struct timeval tv;
   gettimeofday(&tv, NULL);
   srand48(tv.tv_usec);
   srandom(tv.tv_usec);
   /* NOTE: Disabled to allow compilation on 64-bit systems; doesn't seem to cause
      problems. (cge, 091028) if (sizeof(long) != 4) { fprintf(stderr, "Error: sizeof long
      (%d) should be 4\n", (int)sizeof(long int)); exit(-1);
      }
      */
   if (argc == 2 && strcmp(argv[1], "--") == 0) {
      printf("usage: xgrow tilefile [option=#]... \n");
      printf(" tilefile is an input file that specifies tiles\n");
      printf(" options:\n");
      printf("  block=  display block size, 1...10\n");
      printf("  size=   field side length (power-of-two) [default 256]\n");
      printf("  rand=   random number seed\n");
      printf("  T=      threshold T (relative to Gse) for irreversible Tile Assembly "
             "Model\n");
      printf("  k=      hybridization rate constant (/sec)\n");
      printf("  Gmc=    initiation free energy  (units kT)\n");
      printf("  Gse=    interaction free energy per binding\n");
      printf("  Gmch=   initiation free energy  for hydrolyzed units\n");
      printf("  Gseh=   interaction free energy for hydrolyzed units\n");
      printf("  Ghyd=   free energy of hydrolysis\n");
      printf("  Gas=    activation energy for spontaneous hydrolysis\n");
      printf("  Gam=    activation energy for mismatched sticky ends\n");
      printf("  Gae=    activation energy for unmatched sticky ends\n");
      printf("  Gah=    activation energy for hydrolyzed neighbors\n");
      printf("  Gao=    delta a. e. for output vs input-triggers hydrolysis\n");
      printf("  Gfc=    log concentration of flakes (otherwise no depletion)\n");
      printf("  stoic=n[s]            set stoichiometry of tile n to relative value s\n");
      printf("  anneal=g,t            anneal Gse from g to given Gse with time constant "
             "t\n");
      printf("  linan=h,s,C0,Cfin,dt  do a linear anneal, where delta G at each "
             "temperature is\n");
      printf("                        governed by (Delta) h and (Delta) s in kcals per "
             "mol.  \n");
      printf("                        Go from C0 > Cfin changing temps 1 degree every dt "
             "secs \n");
      printf("                        (in incremenets of 0.1C).  Ignores Gse value.\n");
      printf("  seed=i,j,n            seed tile type n at position i,j\n");
      printf("  tinybox=V             use dynamic flakes in a box of volume V (in "
             "liters).\n");
      printf("  addflakes=i,j,n:N@Gfc simulate N separate flakes\n");
      printf("  stripe=o[:p,w]*       width w stripe with p errors, offset o\n");
      printf("  wander                wandering `seed' designation\n");
      printf("  fission               can tile be removed if two flakes result?\n");
      printf("  no_fission             the answer is no [default]\n");
      printf("  chunk_fission         allow pairs & 2x2 blocks to dissociate as one "
             "(implies fission)\n");
      printf("  blast_rate_alpha      square kxk holes are blasted with this per-tile "
             "rate for 1x1 [default=0]\n"
             "                         (rate relative to tile addition, i.e. scaled by "
             "total concentration)\n");
      printf("  blast_rate_beta        rate scales as 1/k^beta [default=4]\n");
      printf(
          "  blast_rate_gamma       rate also scales as exp(-gamma*(k-1)) [default=0]\n");
      printf(
          "  zero_bonds            can tiles be added if they bond with 0 strength?\n");
      printf("  no_zero_bonds          the answer is no [default]\n");
      printf("  min_strength=         set minimum bond energy below which attachments "
             "are considered incorrect\n");
      printf("  periodic              periodic boundary conditions\n");
      printf("  -linear               simulate linear A B tiles, write errs > stdout \n");
      printf("  -nw                   no X window (only if ?max set)\n");
      printf("  update_rate=          update display every so-many events\n");
      printf("  tracefile=            append datafile info (see below) EVERY so-many "
             "events\n");
      printf("  movie                 export MATLAB-format flake array information EVERY "
             "so-many events\n");
      printf("  tmax=                 quit after time t has passed\n");
      printf("  emax=                 quit after e events have occurred\n");
      printf("  smax=                 quit when the fragment or total size of fragments "
             "is size s\n");
      printf("  mmax=                quit when there are at least mm mismatches by "
             "xgrow's interpretation\n");
      printf("  smin=                 quit when the fragment or total size of fragments "
             "goes to or below size s\n");
      printf("  fsmax=                 quit when a single fragment reaches size s\n");
      printf("  untiltiles=           quit when all (numbered) tiles in the "
             "comma-delineated list are in the assembly\n");
      printf("  clean_cycles=         at end, remove how many layers of weakly attached "
             "tiles?\n"
             "                        [default=0]\n");
      printf("  clean_X=              for cleaning, minimal ratio of off-rate to on-rate "
             "[default=1.0]\n");
      printf("  fill_cycles=          at end, add how many layers of strongly attached "
             "tiles?\n"
             "                        [default=0]\n");
      printf("  fill_X=               for filling, minimal ratio of off-rate to on-rate "
             "[default=1.0]\n");
      printf("  error_radius=         when writing to file, #mismatches counts only "
             "those for which \n"
             "                        all surounding tiles are present (after "
             "clean/fill) [default=0]\n");
      printf("  repair_unique_T=      alternative clean/fill, called Rx: remove "
             "mismatches, fill in interior sites \n"
             "                        if there is a unique strength-T tile, then fill in "
             "by strongest tile\n");
      printf("  datafile=             append Gmc, Gse, ratek, time, size, #mismatched "
             "se, events, perimeter, dG, dG_bonds for each flake\n");
      printf("  arrayfile=            output MATLAB-format flake array information on "
             "exit (after cleaning)\n");
      printf("  exportfile=           on-request output of MATLAB-format flake array "
             "information\n");
      printf("                        [defaults to 'xgrow_export_output']\n");
      printf("  importfile=FILENAME   import all flakes from FILENAME.\n");
      printf("  importfile            import all flakes from xgrow_export_output.\n");
      printf("  pause                 start in paused state; wait for user to request "
             "simulation to start.\n");
      printf("  testing               run automated tests instead of a simulation.\n");
      exit(0);
   }
   if (argc == 1) {
      fprintf(stderr,
              "* First argument must be a tile file!\nTry 'xgrow --' for help.\n");
      exit(1);
   }

   tmax = 0;
   emax = 0;
   smax = 0;
   fsmax = 0;
   smin = -1;
   mmax = 0;
   linear = 0;
   datafp = NULL;
   arrayfp = NULL;
   largeflakefp = NULL;
   untiltilescountfp = NULL;

   // reset tile colors and names that haven't been assigned
   for (i = 15; i < MAXTILETYPES; i++)
      tile_colors[i] = NULL;
   for (i = 0; i < MAXTILETYPES; i++)
      params->tile_names[i] = NULL;

   if ((snprintf(tileset_name, 256, "%s", argv[1]), tilefp = fopen(tileset_name, "r")) !=
       NULL)
      read_tilefile(tilefp, params);
   else if ((snprintf(tileset_name, 256, "%s.tiles", argv[1]),
             tilefp = fopen(tileset_name, "r")) != NULL)
      read_tilefile(tilefp, params);
   else {
      fprintf(stderr,
              "* First argument must be a tile file!\nTry 'xgrow --' for help.\n");
      exit(1);
   }

   for (i = 2; i < argc; i++) {
      parse_arg_line(argv[i], params);
   }
   if (tmax == 0 && emax == 0 && smax == 0 && mmax == 0 && fsmax == 0 && smin == -1) {
      fprintf(stderr, "No max setting: forcing UI mode.\n");
      XXX = 1;
   }
   if (params->hydro && params->fission_allowed == 2) {
      fprintf(stderr, "* Current implementation does not allow chunk_fission and "
                      "hydrolysis simultaneously.\n");
      exit(1);
   }
   if (params->double_tiles) {
      // if (fission_allowed != 0) {
      // fprintf(stderr, "Double tiles cannot be used with fission or chunk_fission
      // currently.\n"); exit(0);
      //}
      if (params->hydro) {
         fprintf(stderr, "Double tiles cannot be used with hydrolysis currently.\n");
         exit(1);
      }
      if (fill_cycles || repair_unique) {
         fprintf(stderr,
                 "Double tiles cannot be used with repair or fill cycles currently.\n");
         exit(1);
      }
      if (params->blast_rate > 0) {
         fprintf(stderr, "Blasting cannot be used with double tiles currently.\n");
         exit(1);
      }
   }

   for (params->size_P = 5; (1 << params->size_P) < params->size; params->size_P++)
      ;
   params->size = (1 << params->size_P);
   // if (XXX) {
   //   if (size*block > 800) block=800/size;
   //   if (block==0) { size=512; block=1; size_P=9; }
   //}

   if (params->blast_rate_alpha > 0) {
      int kb;
      fprintf(stderr, "blast_rate: alpha = %f, beta = %f, gamma = %f\n",
              params->blast_rate_alpha, params->blast_rate_beta,
              params->blast_rate_gamma);
      for (kb = 1; kb < params->size; kb++)
         params->blast_rate += params->blast_rate_alpha *
                               exp(-params->blast_rate_gamma * (kb - 1)) /
                               pow(kb * 1.0, params->blast_rate_beta);
      fprintf(stderr, "total blast_rate per site = %f\n", params->blast_rate);
   }

   NROWS = (params->size + NBDY * 2);
   NCOLS = (params->size + NBDY * 2);
   VOLUME = (NROWS * NCOLS);

   params->T = params->T * params->Gse;

   /* for consistency's sake, insert seed_i, seed_j, seed_n, Gfc into fparam */
   /* if no addflake commands were issued */
   if (fparam == NULL) {
      fprm = (struct flake_param *)malloc(sizeof(struct flake_param));
      fprm->seed_i = params->seed_i;
      fprm->seed_j = params->seed_j;
      fprm->seed_n = params->seed_n;
      fprm->import_from = NULL;
      fprm->N = (params->tinybox == 0);
      fprm->Gfc = params->Gfc;
      fprm->next_param = fparam;
      fparam = fprm;
   }

   // make sure all seed tiles are on the board
   for (fprm = fparam; fprm != NULL; fprm = fprm->next_param) {
      while (fprm->seed_i >= params->size)
         fprm->seed_i /= 2;
      while (fprm->seed_j >= params->size)
         fprm->seed_j /= 2;
   }
   if (fparam != NULL) {
      params->seed_i = fparam->seed_i;
      params->seed_j = fparam->seed_j;
   }

   if (!XXX) {
      fprintf(stderr, " Starting simulation (1st seed=%d,%d,%d) on %d x %d board.\n",
              params->seed_i, params->seed_j, params->seed_n, params->size, params->size);
   }
}

void write_datalines(tube_params *params, tube *tp, flake *fp, FILE *out, char *text) {
   flake *fpp;
   int perimeter;
   double dG_bonds;

   for (fpp = tp->flake_list; fpp != NULL; fpp = fpp->next_flake) {
      if (strcmp(text, "") == 0)
         fpp = fp;
      perimeter = calc_perimeter(fpp);
      dG_bonds = calc_dG_bonds(fpp);
      if (tp->hydro)
         fprintf(out, " %f %f %f %f %f %f %f %f %f ", params->Gseh, params->Gmch,
                 params->Ghyd, params->Gas, params->Gam, params->Gae, params->Gah,
                 params->Gao, params->Gfc);
      fprintf(out, " %f %f %f %f %d %d %lld %d %f %f%s", tp->Gmc, tp->Gse, tp->k, tp->t,
              fpp->tiles, fpp->mismatches, tp->events, perimeter, fpp->G, dG_bonds, text);
      if (strcmp(text, "") == 0)
         break;
   }
   fflush(out);
}

void write_largeflakedata(tube *tp, FILE *filep) {
   flake *fp;
   int large_flakes = 0;
   for (fp = tp->flake_list; fp != NULL; fp = fp->next_flake) {
      if (fp->tiles > tthresh) {
         large_flakes++;
      }
   }
   fprintf(filep, "%d\n", large_flakes);
}

void write_untiltilescountdata(tube *tp, FILE *filep) {
   fprintf(filep, "%d\n", tp->untiltilescount);
}

void write_flake(tube_params *params, tube *tp, FILE *filep, char *mode, flake *fp) {
   int n, row, col, tl;
   int size = (1 << fp->P);
   if (filep != NULL) {
      if (strcmp(mode, "flake") == 0)
         n = export_flake_n++;
      else if (strcmp(mode, "movie") == 0)
         n = export_movie_n++;
      else
         n = 1;
      fprintf(filep, "\n%s{%d}={ ...\n", mode, n);
      fprintf(filep, "[ ");
      write_datalines(params, tp, fp, filep, "");
      fprintf(filep, " ],...\n  [");
      for (tl = 1; tl <= tp->N; tl++)
         fprintf(filep, " %8.5f", (tp->conc[tl] > 0) ? -log(tp->conc[tl]) : 0);
      fprintf(filep, " ],...\n  [");
      for (row = 0; row < params->size; row++) {
         for (col = 0; col < params->size; col++)
            fprintf(filep, " %d", fp->Cell(row, col));
         fprintf(filep, "; ...\n");
      }
      fprintf(filep, " ] };\n\n");
      fflush(filep);
   }
}

void export_flake(tube_params *params, tube *tp, char *mode, flake *fp) {
   if (export_fp == NULL)
      export_fp = fopen("xgrow_export_output", "a+");
   fprintf(stderr, "Writing flake #%d...\n", export_flake_n);
   write_flake(params, tp, export_fp, mode, fp);
   fflush(export_fp);
}

/*    READ in by Shaun Lee
 * Count the number of flakes in an exported file.
 * Doesn't actually count the number of flakes, it just watches
 * the flake number and returns the last one it finds.
 */
int count_flakes(FILE *flake_file) {
   int n, i, flake_size, row, garbage, end;
   /*
    * This assumes the parameters (which we want to ignore)
    * are less than 20 charaters long (Overkill.  14 should be enough)
    */
   char line[20];
   int lnum = 0;
   flake_size = 0;

   /* Run through once to make sure the flakes are the right size. */

   fscanf(flake_file, "\nflake{%d}={ ...\n[", &n);
   lnum += 2; // the line number counting is very approximate.

   /* For debugging */
   fprintf(stderr, "Reading flake number: %d\n", n);

   /* Run through the parameters, waiting for ],... to appear twice. */
   for (i = 0; i < 2; i++) {
      while (!feof(flake_file)) {
         fscanf(flake_file, "%s", line);
         if (strcmp(line, "],...") == 0)
            break;
      }
   }
   if (feof(flake_file)) {
      fprintf(stderr, "Reached EOF on flake number %d without reading params.\n", n);
      exit(-1);
   }

   /* Read in the '[' */
   fscanf(flake_file, "%s", line);
   lnum += 3; // for params  (if written by xgrow)

   while (1) {
      fscanf(flake_file, "%s", line);
      // fprintf(stderr, "Just read in loop near line %d: %s\n", lnum, line); //****
      if (strcmp(line, "...") == 0) /* then we've reached the end of the line */
      {
         break;
         lnum++;
      } else
         flake_size++;
   }

   if (flake_size > import_flake_size)
      import_flake_size = flake_size;

   /*
    * Now finish processing the first flake.
    */

   for (row = 2; row <= flake_size; row++) {
      //        fprintf(stderr, "Now on row: %d\n", row);
      for (i = 0; i <= flake_size; i++) {
         fscanf(flake_file, "%d", &garbage);
         /*	  fprintf(stderr, "%d:%d ", i, garbage); */
      }
      fscanf(flake_file, "%s", line);
      /* Make sure that we've actually read the end of line. */
      if (strcmp(";", line) != 0) {
         fprintf(stderr,
                 "Error parsing input file near line %d. Expected ';' but read '%s'\n",
                 lnum, line);
         exit(-1);
      }
      fscanf(flake_file, "%s", line);
      lnum++;
      if (strcmp("...", line) != 0) {
         fprintf(stderr,
                 "Error parsing input file near line %d. Expected '...' but read '%s'\n",
                 lnum, line);
         exit(-1);
      }
   }

   row--; /* stupid off by one things. */

   if (row != flake_size) {
      fprintf(stderr,
              "Error near line %d: Flake dimensions do not match.\n%d rows and "
              "%d columns\n",
              lnum, row, flake_size);
      exit(-1);
   }

   fscanf(flake_file, "%s", line);
   assert(strcmp(line, "]") == 0);
   fscanf(flake_file, "%s", line);
   lnum++;
   assert(strcmp(line, "};") == 0);

   /* Watch for end of file as well as find the next flake number. */
   end = fscanf(flake_file, "\nflake{%d}={ ...\n[", &n);
   lnum += 2;

   /* Now go through the rest of the file. */
   while (end == 1) {
      /* For debugging */
      fprintf(stderr, "Reading flake number: %d\n", n); //****

      /* Run through the parameters, waiting for ],... to appear twice. */
      for (i = 0; i < 2; i++) {
         while (!feof(flake_file)) {
            fscanf(flake_file, "%s", line);
            if (strcmp(line, "],...") == 0)
               break;
         }
      }
      if (feof(flake_file)) {
         fprintf(stderr, "Reached EOF on flake number %d without reading params.\n", n);
         exit(-1);
      }

      fprintf(stderr, " skipped param lines\n");
      fscanf(flake_file, "%s", line);
      lnum += 3;
      assert(strcmp(line, "[") == 0);
      for (row = 1; row <= flake_size; row++) {
         for (i = 0; i < flake_size; i++) {
            fscanf(flake_file, "%s", line);
         }
         fscanf(flake_file, "%s", line);
         lnum++;
         assert(strcmp("...", line) == 0);
      }

      fscanf(flake_file, "%s", line);
      assert(strcmp("]", line) == 0);
      fscanf(flake_file, "%s", line);
      lnum++;
      assert(strcmp("};", line) == 0);

      /* Watch for end of file as well as find the next flake number. */
      end = fscanf(flake_file, "\nflake{%d}={ ...\n[", &n);
      lnum += 2;
   } /* end while */

   fprintf(stderr, "Found %d flakes\n", n);

   return n;
}

/* Import flakes by Shaun Lee*/
/*
 * Import flake data for flake n from the file and store it into
 * flake fp, then recalc_g at the end.  If there is no corresponding flake
 * in the file, it doesn't modify the flake.
 */
void import_flake(tube *tp, flake *current_flake, FILE *flake_file, int flake_number) {
   fpos_t flake_start;
   int tile_type, read_flake_number, translate_i, translate_j, i, j, end, flake_size;
   char line[20];
   int lnum = 0;
   flake_size = 0;

   /* Just to be sure. */
   rewind(flake_file);

   /* Find the flake we want. */
   fscanf(flake_file, "\nflake{%d}={ ...\n[", &read_flake_number);
   lnum++; // the line number counting is very approximate.

   while (read_flake_number != flake_number) {
      /* Skip through the parameters, waiting for ],... to appear twice. */
      for (i = 0; i < 2; i++) {
         while (1) {
            /* Assumes the parameters are terminated by '],...' */
            if (fscanf(flake_file, "%s", line) == 0) {
               fprintf(stderr,
                       "Error near line %d: Expected ],... at end of parameters.\n",
                       lnum + 1);
               exit(-1);
            }
            lnum++;
            if (strcmp(line, "],...") == 0)
               break;
         }
      }

      /* Skip through the data for this flake */
      while (1) {
         if (fscanf(flake_file, "%s", line) == 0) {
            fprintf(stderr, "Error near line %d: Expected }; at end of data.\n",
                    lnum + 1);
            exit(-1);
         }
         lnum++;
         if (strcmp(line, "};") == 0)
            break;
      }

      end = fscanf(flake_file, "\n\nflake{%d}={", &read_flake_number);

      if (!end) {
         fprintf(stderr, "Flake %d not found in file (near line %d).\n", flake_number,
                 lnum);
         return;
      }
      lnum += 2;
   }

   /* Skip through the parameters, waiting for ],... to appear twice. */
   for (i = 0; i < 2; i++) {
      while (1) {
         fscanf(flake_file, "%s", line);
         lnum++;
         if (strcmp(line, "],...") == 0)
            break;
      }
   }
   /* Read in the starting '[' */
   fscanf(flake_file, "%s", line);
   lnum++;
   assert(strcmp("[", line) == 0);

   /*
    * Save the current position in the file stream so that we can
    * import data values after calculating the size of the import
    * flake.
    */
   fgetpos(flake_file, &flake_start);

   /* Find the width of the flake. */
   while (1) {
      fscanf(flake_file, "%s", line);
      if (strcmp(line, "...") == 0)
         break;
      else
         flake_size++;
   }

   /* Go back to the starting point. */
   fsetpos(flake_file, &flake_start);

   /*
    * We want to import the flake into the middle of the empty flake, so
    * calculate how much we need to translate every cell by.
    */
   int size = (1 << current_flake->P);
   translate_i = (size - flake_size) / 2 + import_offset_i;
   translate_j = (size - flake_size) / 2 + import_offset_j;

   /* we already plopped down a seed tile.  erase it if the loaded assemble won't. */
   change_cell(current_flake, current_flake->seed_i, current_flake->seed_j, 0);

   for (i = 0; i < flake_size; i++) {
      for (j = 0; j < (flake_size - 1); j++) {
         fscanf(flake_file, "%d", &tile_type);
         change_cell(current_flake, translate_i + i, translate_j + j, tile_type);
      }
      fscanf(flake_file, "%d;", &tile_type);
      change_cell(current_flake, translate_i + i, translate_j + j, tile_type);

      fscanf(flake_file, "%s", line);
      lnum++;
      assert(strcmp(line, "...") == 0);
   }

   /* Now choose a random active cell and set it as the flake's seed. */
   srand(time(0));
   i = flake_size * (((double)rand()) / ((double)RAND_MAX)) + translate_i;
   j = flake_size * (((double)rand()) / ((double)RAND_MAX)) + translate_j;
   while ((current_flake->Cell(i, j)) == 0 ||
          tp->dt_left[current_flake->Cell(i, j)]) // FIXME: vdouble
   {
      i = flake_size * (((double)rand()) / ((double)RAND_MAX)) + translate_i;
      j = flake_size * (((double)rand()) / ((double)RAND_MAX)) + translate_j;
   }
   current_flake->seed_i = i;
   current_flake->seed_j = j;
   current_flake->seed_n = current_flake->Cell(i, j);

   /*    fprintf(stderr, "Just set seed for flake %d to (%d,%d), tile type %d\n",
    * flake_number, i, j, current_flake->Cell(i,j)); */
   current_flake->events = 0;

   recalc_G(current_flake);
}

void closeargs(tube *tp, tube_params *params) {
   int i;
   flake *fpp;

   // cleans all flakes  (removes "temporary" tiles on growth edge)
   for (fpp = tp->flake_list; fpp != NULL; fpp = fpp->next_flake)
      clean_flake(fpp, clean_X, clean_cycles);

   // fills in all flakes  (adds tiles in holes or on growth edge that "ought" to be
   // there)
   for (fpp = tp->flake_list; fpp != NULL; fpp = fpp->next_flake)
      fill_flake(fpp, fill_X, fill_cycles);

   // alternative fix-up routine; redundant with clean/fill, but supposedly more reliable
   if (repair_unique)
      for (fpp = tp->flake_list; fpp != NULL; fpp = fpp->next_flake)
         repair_flake(fpp, repair_unique_T, tp->Gse);

   // recalc's all mismatch #'s if error_radius is given -- keep only those not near an
   // empty cell.
   if (error_radius > .5) {
      for (fpp = tp->flake_list; fpp != NULL; fpp = fpp->next_flake)
         error_radius_flake(fpp, error_radius);
   }

   /* output information for *all* flakes */
   if (datafp != NULL) {
      write_datalines(params, tp, tp->flake_list, datafp, "\n");
      fclose(datafp);
   }

   if (largeflakefp != NULL) {
      write_largeflakedata(tp, largeflakefp);
      fclose(largeflakefp);
   }

   if (untiltilescountfp != NULL) {
      write_untiltilescountdata(tp, untiltilescountfp);
      fclose(untiltilescountfp);
   }
   if (arrayfp != NULL) {
      export_flake_n = 1;
      for (fpp = tp->flake_list; fpp != NULL; fpp = fpp->next_flake) {
         write_flake(params, tp, arrayfp, "flake", fpp);
      }
      fclose(arrayfp);
   }
   if (export_fp != NULL)
      fclose(export_fp);

   // free memory
   free(params->stoic);

   while (fparam != NULL) {
      fprm = fparam->next_param;
      free(fparam);
      fparam = fprm;
   }

   free_tube(tp); // this frees all the flakes too
}

// (i,j) should not be empty. returns 0 if OK, 1 if mismatches with E or S or N or W, or
// unbound s.e.
#define errortile(i, j)                                                                  \
   (((tp->tileb[fp->Cell(i, j)][1] != tp->tileb[fp->Cell(i, (j) + 1)][3] &&              \
      tp->tileb[fp->Cell(i, j)][1] * tp->tileb[fp->Cell(i, (j) + 1)][3] > 0 &&           \
      (fp->tube->glue)[tp->tileb[fp->Cell(i, j)][1]]                                     \
                      [tp->tileb[fp->Cell(i, (j) + 1)][3]] < min_strength) ||            \
     (tp->tileb[fp->Cell(i, j)][3] != tp->tileb[fp->Cell(i, (j)-1)][1] &&                \
      tp->tileb[fp->Cell(i, j)][3] * tp->tileb[fp->Cell(i, (j)-1)][1] > 0 &&             \
      (fp->tube->glue)[tp->tileb[fp->Cell(i, j)][3]][tp->tileb[fp->Cell(i, (j)-1)][1]] < \
          min_strength) ||                                                               \
     (tp->tileb[fp->Cell(i, j)][2] != tp->tileb[fp->Cell((i) + 1, j)][0] &&              \
      tp->tileb[fp->Cell(i, j)][2] * tp->tileb[fp->Cell((i) + 1, (j))][0] >              \
          0 && /* This one works!? */                                                    \
      (fp->tube->glue)[tp->tileb[fp->Cell(i, j)][2]]                                     \
                      [tp->tileb[fp->Cell((i) + 1, j)][0]] < min_strength) ||            \
     (tp->tileb[fp->Cell(i, j)][0] != tp->tileb[fp->Cell((i)-1, j)][2] &&                \
      tp->tileb[fp->Cell(i, j)][0] * tp->tileb[fp->Cell((i)-1, (j))][2] > 0 &&           \
      (fp->tube->glue)[tp->tileb[fp->Cell(i, j)][0]][tp->tileb[fp->Cell((i)-1, j)][2]] < \
          min_strength))                                                                 \
        ? 1                                                                              \
        : 0)

int getcolor(tube *tp, flake *fp, int i, int j, int err) {
   int min_strength = tp->min_strength;
   int size = (1 << fp->P);
   return (
       (fp->Cell(i, j) == 0)
           ? translate[0]
           : ((err == 1)
                  ? (errortile(i, j) ? errorcolor : goodcolor)
                  : ((err == 2) ? ((fp->Cell(i, j) > fp->N / 2)
                                       ? (errortile(i, j) ? hydroerrcolor : hydrocolor)
                                       : (errortile(i, j) ? errorcolor : goodcolor))
                                : translate[fp->Cell(i, j)])));
}

#define getcolordij(tp, fp, row, col, di, dj, err)                                       \
   ((((row) + (di)) >= 0 && ((row) + (di)) < size && ((col) + (dj)) >= 0 &&              \
     ((col) + (dj)) < size)                                                              \
        ? getcolor(tp, fp, (row) + (di), (col) + (dj), err)                              \
        : 0)

/* NOTE: requires 2^P < NCOLS+2*NBDY */
void showpic(tube *tp, flake *fp,
             int err) /* display the field */ // err param is ignored!
{
   int row, col, i1, i2, color, j, j1, j2, blocktop = block;
   char *picture = (*spinimage).data;
   static int last_display_type = 0;
   int new_display = (last_display_type != (10 * errors + errorc));
   last_display_type = (10 * errors + errorc); // re-draw everything when colormap changes
   int size = (1 << fp->P);
   if (block > 4)
      blocktop = block - 1;
   if (8 == (*spinimage).depth) {
      if (block > 1) /* I wish I knew how to do this faster */
         for (row = 0; row < size; row++)
            for (col = 0; col < size; col++) {
               color = getcolor(tp, fp, row, col, err);
               j = block * ((col + NBDY) + block * NCOLS * (row + NBDY));
               if (color != picture[j] || new_display) {
                  for (i1 = 0; i1 < blocktop; i1++) {
                     j1 = i1 * block * NCOLS + j;
                     for (i2 = 0; i2 < blocktop; i2++)
                        picture[j1 + i2] = color;
                  }
                  if (block > 4) {
                     int n = fp->Cell(row, col);
                     int Ccolm = lightcolor, Ccolp = lightcolor, Crowp = lightcolor,
                         Crowm = lightcolor;
                     if (errors == 1) {
                        int ncolm = fp->Cell(row, col - 1),
                            ncolp = fp->Cell(row, col + 1);
                        int nrowm = fp->Cell(row - 1, col),
                            nrowp = fp->Cell(row + 1, col);
                        Ccolm = weakcolor;
                        Ccolp = weakcolor;
                        Crowp = weakcolor;
                        Crowm = weakcolor;
                        if (tp->Gse_EW[n][ncolm] > 1.5 * tp->Gse)
                           Ccolm = strongcolor;
                        if (tp->Gse_EW[ncolp][n] > 1.5 * tp->Gse)
                           Ccolp = strongcolor;
                        if (tp->Gse_NS[n][nrowp] > 1.5 * tp->Gse)
                           Crowp = strongcolor;
                        if (tp->Gse_NS[nrowm][n] > 1.5 * tp->Gse)
                           Crowm = strongcolor;
                        if (tp->Gse_EW[n][ncolm] < 0.5 * tp->Gse)
                           Ccolm = nullcolor;
                        if (tp->Gse_EW[ncolp][n] < 0.5 * tp->Gse)
                           Ccolp = nullcolor;
                        if (tp->Gse_NS[n][nrowp] < 0.5 * tp->Gse)
                           Crowp = nullcolor;
                        if (tp->Gse_NS[nrowm][n] < 0.5 * tp->Gse)
                           Crowm = nullcolor;
                     }
                     if (n == 0 && fp->Cell(row, col - 1) == 0)
                        Ccolm = translate[0];
                     if (n == 0 && fp->Cell(row, col + 1) == 0)
                        Ccolp = translate[0];
                     if (n == 0 && fp->Cell(row - 1, col) == 0)
                        Crowm = translate[0];
                     if (n == 0 && fp->Cell(row + 1, col) == 0)
                        Crowp = translate[0];
                     for (i1 = 0; i1 < blocktop; i1++) // col-1 side
                        picture[i1 * block * NCOLS + j - 1] = Ccolm;
                     for (i1 = 0; i1 < blocktop; i1++) // col+1 side
                        picture[i1 * block * NCOLS + j + blocktop] = Ccolp;
                     for (i1 = 0; i1 < blocktop; i1++) // row+1 side
                        picture[blocktop * block * NCOLS + j + i1] = Crowp;
                     for (i1 = 0; i1 < blocktop; i1++) // row-1 side
                        picture[-block * NCOLS + j + i1] = Crowm;
                  }
               }
            }
      else {
         for (row = 0, j = NCOLS * NBDY + NBDY; row < size; row++, j += 2 * NBDY)
            for (col = 0; col < size; col++, j++)
               picture[j] = getcolor(tp, fp, row, col, err);
      }
   } else {          /* depth is not == 8, use xputpixel (this is really ugly) */
      if (block > 1) /* I wish I knew how to do this faster */
         for (row = 0; row < size; row++)
            for (col = 0; col < size; col++) {
               color = getcolor(tp, fp, row, col, err);
               if (color != XGetPixel(spinimage, j1 = block * (col + NBDY),
                                      j2 = block * (row + NBDY)) ||
                   new_display) {
                  for (i2 = 0; i2 < blocktop; i2++)
                     for (i1 = 0; i1 < blocktop; i1++)
                        XPutPixel(spinimage, j1 + i1, j2 + i2, color);
                  if (block > 4) {
                     for (i1 = 0; i1 < blocktop; i1++) {
                        int n = fp->Cell(row, col);
                        int Ccolm = lightcolor, Ccolp = lightcolor, Crowp = lightcolor,
                            Crowm = lightcolor;
                        if (errors == 1) {
                           int ncolm = fp->Cell(row, col - 1),
                               ncolp = fp->Cell(row, col + 1);
                           int nrowm = fp->Cell(row - 1, col),
                               nrowp = fp->Cell(row + 1, col);
                           Ccolm = weakcolor;
                           Ccolp = weakcolor;
                           Crowp = weakcolor;
                           Crowm = weakcolor;
                           if (tp->Gse_EW[n][ncolm] > 1.5 * tp->Gse)
                              Ccolm = strongcolor;
                           if (tp->Gse_EW[ncolp][n] > 1.5 * tp->Gse)
                              Ccolp = strongcolor;
                           if (tp->Gse_NS[n][nrowp] > 1.5 * tp->Gse)
                              Crowp = strongcolor;
                           if (tp->Gse_NS[nrowm][n] > 1.5 * tp->Gse)
                              Crowm = strongcolor;
                           if (tp->Gse_EW[n][ncolm] < 0.5 * tp->Gse)
                              Ccolm = nullcolor;
                           if (tp->Gse_EW[ncolp][n] < 0.5 * tp->Gse)
                              Ccolp = nullcolor;
                           if (tp->Gse_NS[n][nrowp] < 0.5 * tp->Gse)
                              Crowp = nullcolor;
                           if (tp->Gse_NS[nrowm][n] < 0.5 * tp->Gse)
                              Crowm = nullcolor;
                        }
                        if (n == 0 && fp->Cell(row, col - 1) == 0)
                           Ccolm = translate[0];
                        if (n == 0 && fp->Cell(row, col + 1) == 0)
                           Ccolp = translate[0];
                        if (n == 0 && fp->Cell(row - 1, col) == 0)
                           Crowm = translate[0];
                        if (n == 0 && fp->Cell(row + 1, col) == 0)
                           Crowp = translate[0];
                        XPutPixel(spinimage, j1 - 1, j2 + i1, Ccolm);        // col-1 side
                        XPutPixel(spinimage, j1 + blocktop, j2 + i1, Ccolp); // col+1 side
                        XPutPixel(spinimage, j1 + i1, j2 + blocktop, Crowp); // row+1 side
                        XPutPixel(spinimage, j1 + i1, j2 - 1, Crowm);        // row-1 side
                     }
                  }
               }
            }
      else
         for (row = 0; row < size; row++)
            for (col = 0; col < size; col++) {
               color = getcolor(tp, fp, row, col, err);
               XPutPixel(spinimage, col + NBDY, row + NBDY, color);
            }
   }
   XPutImage(display, playground, gc, spinimage, 0, 0, 0, 0, block * NCOLS,
             block * NROWS);
   return;
}

/* NOTE: requires 2^P < NCOLS+2*NBDY */
void add_sample_pic(tube *tp, flake *fp, int err) /* add sample the field */
{
   int row, col, i1, i2, di = 0, ddi, j1, j2, dj = 0, ddj, blocktop = block;
   int color, oldcolor, n_tries = 10000, n, collision = 0, anything = 0;
   int size = (1 << fp->P);
   // we will try n_tries random offsets before giving up on
   // avoiding collision.
   // this gets a di, dj
   XDrawImageString(display, samplebutton, gcr, 0, font_height, "  SAMPLING   ", 13);
   XPutImage(display, playground, gc, spinimage, 0, 0, 0, 0, block * NCOLS,
             block * NROWS);

   for (n = 1; n < n_tries && (collision == 1 || anything == 0); n++) {
      di = random() % (2 * size - 1) - size;
      dj = random() % (2 * size - 1) - size;
      //   fprintf(stderr, "trying to place flake %d at %d %d\n", fp, di, dj);
      collision = 0;
      anything = 0;
      for (row = 0; row < size; row++)
         for (col = 0; col < size; col++) {
            color = 0;
            for (ddi = -1; ddi <= 1; ddi++)
               for (ddj = -1; ddj <= 1; ddj++)
                  color =
                      MAX(getcolordij(tp, fp, row, col, di + ddi, dj + ddj, err), color);
            oldcolor = XGetPixel(spinimage, j1 = block * (col + NBDY),
                                 j2 = block * (row + NBDY));
            if (oldcolor > 0 && color > 0)
               collision = 1;
            if (getcolordij(tp, fp, row, col, di, dj, err) > 0)
               anything++;
         }
      anything = anything && (anything == fp->tiles);
   }
   // if (n==n_tries)  fprintf(stderr, "---GAVE UP---\n"); else fprintf(stderr,
   // "!!!SUCCESS!!!\n");

   if (collision == 0) {
      if (block > 4)
         blocktop = block - 1;
      // always do this the slow way!
      for (row = 0; row < size; row++)
         for (col = 0; col < size; col++) {
            oldcolor = XGetPixel(spinimage, j1 = block * (col + NBDY),
                                 j2 = block * (row + NBDY));
            color = MAX(getcolordij(tp, fp, row, col, di, dj, err), oldcolor);
            if (color != oldcolor) {
               for (i2 = 0; i2 < blocktop; i2++)
                  for (i1 = 0; i1 < blocktop; i1++)
                     XPutPixel(spinimage, j1 + i1, j2 + i2, color);
               if (block > 4) {
                  for (i1 = 0; i1 < blocktop; i1++) {
                     XPutPixel(spinimage, j1 - 1, j2 + i1, lightcolor); // col-1 side
                     XPutPixel(spinimage, j1 + blocktop, j2 + i1,
                               lightcolor); // col+1 side
                     XPutPixel(spinimage, j1 + i1, j2 + blocktop,
                               lightcolor);                             // row+1 side
                     XPutPixel(spinimage, j1 + i1, j2 - 1, lightcolor); // row-1 side
                  }
               }
            }
         }
   }
   XPutImage(display, playground, gc, spinimage, 0, 0, 0, 0, block * NCOLS,
             block * NROWS);
   XDrawImageString(display, samplebutton, gcr, 0, font_height, "   sample    ", 13);
   return;
}

/* NOTE: requires 2^P < NCOLS+2*NBDY */
void showphase(tube_params *params) /* replace tiles by phase diagram T=1 & T=2 lines */
{
   int row, col, color, i1, i2, j1, j2, blocktop = block;
   if (block > 4)
      blocktop = block - 1;
   for (row = 0; row < params->size; row++)
      for (col = 0; col < params->size; col++) {
         color = ((params->size - row) / 2 == col || (params->size - row) == col)
                     ? translate[7]
                     : translate[0];
         j1 = block * (col + NBDY);
         j2 = block * (row + NBDY);
         for (i2 = 0; i2 < blocktop; i2++)
            for (i1 = 0; i1 < blocktop; i1++)
               XPutPixel(spinimage, j1 + i1, j2 + i2, color);
      }
   XPutImage(display, playground, gc, spinimage, 0, 0, 0, 0, block * NCOLS,
             block * NROWS);
   return;
}

/* fix up the pause button */
void setpause(int value) {
   paused = value;
   if (paused == 0)
      sampling = 0;
   if (paused)
      XDrawImageString(display, pausebutton, gcr, 0, font_height, "  run/PAUSE  ", 13);
   else
      XDrawImageString(display, pausebutton, gcr, 0, font_height, "  RUN/pause  ", 13);
}

/* fix up the colors button */
void settilecolor(tube *tp, int value) {
   errorc = value;
   if (tp->hydro) {
      if (errorc == 2)
         XDrawImageString(display, colorbutton, gcr, 0, font_height, "tile/err/HYD", 12);
      else if (errorc == 1)
         XDrawImageString(display, colorbutton, gcr, 0, font_height, "tile/ERR/hyd", 12);
      else if (errorc == 0)
         XDrawImageString(display, colorbutton, gcr, 0, font_height, "TILE/err/hyd", 12);
   } else {
      if (errorc == 1)
         XDrawImageString(display, colorbutton, gcr, 0, font_height, "  tile/ERR  ", 12);
      else if (errorc == 0)
         XDrawImageString(display, colorbutton, gcr, 0, font_height, "  TILE/err  ", 12);
   }
}

void setsidecolor(int value) {
   errors = value;
   if (errors == 1)
      XDrawImageString(display, sidebutton, gcr, 0, font_height, "  box/BONDS ", 12);
   else if (errors == 0)
      XDrawImageString(display, sidebutton, gcr, 0, font_height, "  BOX/bonds ", 12);
}

/* fix up the export button */
void setexport(int value) {
   export_mode = value;
   if (export_mode == 2)
      if (export_movie == 0)
         XDrawImageString(display, exportbutton, gcr, 0, font_height, "export[MOVIE]",
                          13);
      else
         switch (random() % 3) {
         case 0:
            XDrawImageString(display, exportbutton, gcr, 0, font_height, "ExPoRt[MOVIE]",
                             13);
            break;
         case 1:
            XDrawImageString(display, exportbutton, gcr, 0, font_height, "eXpOrT[MOVIE]",
                             13);
            break;
         case 2:
            XDrawImageString(display, exportbutton, gcr, 0, font_height, "EXPORT[MOVIE]",
                             13);
            break;
         }
   else if (export_mode == 1)
      XDrawImageString(display, exportbutton, gcr, 0, font_height, "export [ALL] ", 13);
   else if (export_mode == 0)
      XDrawImageString(display, exportbutton, gcr, 0, font_height, "export [ONE] ", 13);
}

/* fix up the seed mode button */
void setwander(tube *tp, int value) {
   tp->wander = value;
   if (tp->wander)
      XDrawImageString(display, seedbutton, gcr, 0, font_height, "fixed/WANDER", 12);
   else
      XDrawImageString(display, seedbutton, gcr, 0, font_height, "FIXED/wander", 12);
}

/* fix up the fission mode button */
void setfission(tube *tp, int value) {
   tp->fission_allowed = value;
   if (tp->hydro && value == 2)
      tp->fission_allowed = 0;
   if (tp->fission_allowed > 0)
      if (tp->fission_allowed == 1)
         XDrawImageString(display, fissionbutton, gcr, 0, font_height, " fission  OK ",
                          13);
      else
         XDrawImageString(display, fissionbutton, gcr, 0, font_height, "chunk fission",
                          13);
   else
      XDrawImageString(display, fissionbutton, gcr, 0, font_height, " NO  fission ", 13);
}

/* this fixes the window up whenever it is uncovered */
void repaint(tube_params *params, tube *tp, flake *fp) {
   int i = 0;
   XDrawString(display, quitbutton, gcr, 0, font_height, "    quit     ", 13);
   XDrawString(display, restartbutton, gcr, 0, font_height, "   restart   ", 13);
   XDrawString(display, cleanbutton, gcr, 0, font_height, "clean/fill/Rx", 13);
   XDrawString(display, samplebutton, gcr, 0, font_height, "   sample    ", 13);
   XDrawString(display, flakebutton, gcr, 0, font_height, "next/big/prev", 13);
   XDrawString(display, tempbutton, gcr, 0, font_height, " cool   heat ", 13);
   setexport(export_mode);
   setpause(paused);
   settilecolor(tp, errorc);
   setsidecolor(errors);
   setwander(tp, tp->wander);
   setfission(tp, tp->fission_allowed);

   /* write various strings */
   if (fp) {
      snprintf(stringbuffer, BUFFER_LEN, "flake %d (%d by %d%s, seed %d @ (%d,%d))",
               fp->flake_ID, (1 << fp->P), (1 << fp->P), tp->periodic ? ", periodic" : "",
               fp->seed_n, fp->seed_i, fp->seed_j);
      XDrawImageString(display, window, gc, 5, (++i) * font_height, stringbuffer,
                       strlen(stringbuffer));
      snprintf(stringbuffer, BUFFER_LEN, "%lld events, %d tiles, %d mismatches         ",
               fp->events, fp->tiles, fp->mismatches);
      XDrawImageString(display, window, gc, 5, (++i) * font_height, stringbuffer,
                       strlen(stringbuffer));
   }

   snprintf(stringbuffer, BUFFER_LEN, "([DX] = %g uM, T = %5.3f C, 5-mer s.e.)    ",
            1000000.0 * 20.0 * exp(-tp->Gmc), 4000 / (tp->Gse / 5 + 11) - 273.15);
   XDrawImageString(display, window, gc, 5, (++i) * font_height, stringbuffer,
                    strlen(stringbuffer));

   if (tp->T > 0)
      snprintf(stringbuffer, BUFFER_LEN, "Gmc=%4.1f  Gse=%4.1f  k=%6.0f   T=%4.1f",
               tp->Gmc, tp->Gse, tp->k, tp->T / tp->Gse);
   else
      snprintf(stringbuffer, BUFFER_LEN, "Gmc=%4.1f  Gse=%4.1f  k=%6.0f   ", tp->Gmc,
               tp->Gse, tp->k);
   XDrawImageString(display, window, gc, 5, (++i) * font_height, stringbuffer,
                    strlen(stringbuffer));

   if (tp->hydro) {
      snprintf(stringbuffer, BUFFER_LEN, "Gmch=%4.1f  Gseh=%4.1f  Ghyd=%4.1f",
               params->Gmch, params->Gseh, params->Ghyd);
      XDrawString(display, window, gc, 5, (++i) * font_height, stringbuffer,
                  strlen(stringbuffer));
      snprintf(stringbuffer, BUFFER_LEN,
               "Gas=%4.1f Gam=%4.1f Gae=%4.1f Gah=%4.1f Gao=%4.1f", params->Gas,
               params->Gam, params->Gae, params->Gah, params->Gao);
      XDrawString(display, window, gc, 5, (++i) * font_height, stringbuffer,
                  strlen(stringbuffer));
   }
   if (fp) {
      snprintf(stringbuffer, BUFFER_LEN, "t = %12.3f sec; G = %12.3f      ", tp->t,
               (tp->wander || tp->conc[fp->seed_n] == 0)
                   ? fp->G
                   : (fp->G + log(tp->conc[fp->seed_n])));
   }
   XDrawImageString(display, window, gc, 5, (++i) * font_height, stringbuffer,
                    strlen(stringbuffer));
   snprintf(stringbuffer, BUFFER_LEN,
            "%lld events (%llda,%lldd,%lldh,%lldf), %lld tiles total %s      ",
            tp->events, tp->stat_a, tp->stat_d, tp->stat_h, tp->stat_f,
            tp->stat_a - tp->stat_d + tp->num_flakes, tp->ewrapped ? "[wrapped]" : "");
   XDrawImageString(display, window, gc, 5, (++i) * font_height, stringbuffer,
                    strlen(stringbuffer));

   if (fp && fp->flake_conc > 0) {
      int tl;
      snprintf(stringbuffer, BUFFER_LEN, "Gfc=%4.1f; Gmc=[", -log(fp->flake_conc));
      for (tl = 1; tl <= tp->N && strlen(stringbuffer) < 230; tl++)
         sprintf(stringbuffer + strlen(stringbuffer), " %4.1f",
                 (tp->conc[tl] > 0) ? -log(tp->conc[tl]) : 0);
      sprintf(stringbuffer + strlen(stringbuffer), " ]");
      XDrawImageString(display, window, gc, 5, (++i) * font_height, stringbuffer,
                       strlen(stringbuffer));
   }

   XDrawString(display, window, gc, PLAY_MARGIN, WINDOWHEIGHT - 2 - font->descent,
               "left: identify, middle: puncture, right: Gmc/Gse | EW, RS, CGE '98-'21",
               70);

   if (!sampling && fp)
      showpic(tp, fp, errorc);
   else
      XPutImage(display, playground, gc, spinimage, 0, 0, 0, 0, block * NCOLS,
                block * NROWS);
}

/* a lot of this is taken from the basicwin program in the
   Xlib Programming Manual */
void openwindow(int argc, char **argv) {
   char *window_name;
   char *icon_name = "xgrow";
   Pixmap icon_pixmap;
   char *display_name = NULL;
   XEvent report;
   XColor xcolor, colorcell;
   Colormap cmap;
   char *buffer;
   int i, j;
#define icon_bitmap_width 16
#define icon_bitmap_height 16
   static char icon_bitmap_bits[] = {0x1f, 0xf8, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0x88,
                                     0x1f, 0x88, 0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8,
                                     0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0xff, 0xff,
                                     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};

   window_name = (char *)malloc(256 + 10);
   snprintf(window_name, 256 + 10, "xgrow: %s", tileset_name);

   /* open up the display */
   if ((display = XOpenDisplay(display_name)) == NULL) {
      fprintf(stderr, "%s: cannot connect to X server %s\n", progname,
              XDisplayName(display_name));
      exit(-1);
   }

   /* In order, we try here:
      1. Any user-specified font (set by font=), which is done via font_selection_string.
         If the user didn't specify anything, this has a default values as well (9x15).
      2. 6x13.

      This might seem old-fashioned, but on modern systems, bitmap font availability is
      actually getting smaller and less standardized. */
   if ((font = XLoadQueryFont(display, font_selection_string)) == NULL) {
      fprintf(stderr, "%s: Cannot open %s font\n", progname, font_selection_string);
      if ((font = XLoadQueryFont(display, "6x13")) == NULL) {
         fprintf(stderr, "%s: Cannot open fallback 6x13 font\n", progname);
         exit(-1);
      }
   }

   int throwaway;

   XCharStruct maxchar;
   XTextExtents(font, "EXPORT[MOVIE]", 13, &throwaway, &throwaway, &throwaway, &maxchar);

   font_height = font->ascent + font->descent;

// 4 is 2*BUTTON_BORDER
#define BUTTON_SEP (4 + (font_height / 3))
#define BUTTON_WIDTH (maxchar.width + 2)
#define BUTTON_HEIGHT (font_height + 4)
#define BOTTOM_MARGIN (font_height + 4)
#define BUTTON_BORDER 2

   WINDOWWIDTH = (MAX(block * NCOLS, 256) + 3 * PLAY_MARGIN + BUTTON_WIDTH);
   WINDOWHEIGHT = (3 * PLAY_MARGIN + font_height * NUM_INFO_LINES + BOTTOM_MARGIN +
                   MAX(MAX(block * NROWS, 256), 12 * BUTTON_HEIGHT + 11 * BUTTON_SEP));

   screen = DefaultScreen(display);
   depth = DefaultDepth(display, screen);
   cmap = DefaultColormap(display, screen);
   /* color? This is not the right way to do it, but .... */
   if (1 == depth) {
      fprintf(stderr, "Sorry but this program needs a color monitor.\n");
      exit(-1);
   }
   black = BlackPixel(display, screen);
   white = WhitePixel(display, screen);
   if (XAllocNamedColor(display, cmap, "firebrick", &colorcell, &xcolor))
      darkcolor = colorcell.pixel;
   if (XAllocNamedColor(display, cmap, "wheat", &colorcell, &xcolor))
      lightcolor = colorcell.pixel;
   if (XAllocNamedColor(display, cmap, "red", &colorcell, &xcolor))
      nullcolor = colorcell.pixel;
   if (XAllocNamedColor(display, cmap, "wheat", &colorcell, &xcolor))
      weakcolor = colorcell.pixel;
   if (XAllocNamedColor(display, cmap, "green", &colorcell, &xcolor))
      strongcolor = colorcell.pixel;
   if (XAllocNamedColor(display, cmap, "red", &colorcell, &xcolor))
      errorcolor = colorcell.pixel;
   if (XAllocNamedColor(display, cmap, "green", &colorcell, &xcolor))
      goodcolor = colorcell.pixel;
   if (XAllocNamedColor(display, cmap, "gold", &colorcell, &xcolor))
      hydrocolor = colorcell.pixel;
   if (XAllocNamedColor(display, cmap, "pink", &colorcell, &xcolor))
      hydroerrcolor = colorcell.pixel;
   if (XAllocNamedColor(display, cmap, "black", &colorcell, &xcolor))
      translate[0] = colorcell.pixel;

   for (i = 1; i < MAXTILETYPES; i++)
      if (tile_colors[i] != NULL &&
          XAllocNamedColor(display, cmap, tile_colors[i], &colorcell, &xcolor))
         translate[i] = colorcell.pixel;
      else
         translate[i] = translate[(i > 14) ? ((i - 1) % 14 + 1) : 0];

   /* make the main window */
   window = XCreateSimpleWindow(display, RootWindow(display, screen), 100, 0, WINDOWWIDTH,
                                WINDOWHEIGHT, 4, black, lightcolor);

   /* make the icon */
   icon_pixmap = XCreateBitmapFromData(display, window, icon_bitmap_bits,
                                       icon_bitmap_width, icon_bitmap_height);

   size_hints.flags = PPosition | PSize | PMinSize;
   size_hints.min_width = WINDOWWIDTH;
   size_hints.min_height = WINDOWHEIGHT;
#ifdef X11R3
   size_hints.x = x;
   size_hints.y = y;
   size_hints.width = WINDOWWIDTH;
   size_hints.height = WINDOWHEIGHT;
   XSetStandardProperties(display, window, window_name, icon_name, icon_pixmap, argv,
                          argc, &size_hints);
#else
   {
      XWMHints wm_hints;
      XClassHint class_hints;
      XTextProperty windowName, iconName;
      if (XStringListToTextProperty(&window_name, 1, &windowName) == 0) {
         fprintf(stderr, "%s: structure allocation for windowName failed.\n", progname);
         exit(-1);
      }
      if (XStringListToTextProperty(&icon_name, 1, &iconName) == 0) {
         fprintf(stderr, "%s: structure allocation for iconName failed.\n", progname);
         exit(-1);
      }
      wm_hints.initial_state = NormalState;
      wm_hints.input = True;
      wm_hints.icon_pixmap = icon_pixmap;
      wm_hints.flags = StateHint | IconPixmapHint | InputHint;
      class_hints.res_name = progname;
      class_hints.res_class = "Basicwin";
      XSetWMProperties(display, window, &windowName, &iconName, argv, argc, &size_hints,
                       &wm_hints, &class_hints);
   }
#endif

   /* make the buttons */
   i = 1;
   quitbutton =
       XCreateSimpleWindow(display, window, WINDOWWIDTH - BUTTON_WIDTH - PLAY_MARGIN,
                           WINDOWHEIGHT - PLAY_MARGIN - BOTTOM_MARGIN - i * BUTTON_HEIGHT,
                           BUTTON_WIDTH, BUTTON_HEIGHT, BUTTON_BORDER, black, darkcolor);
   i++;
   restartbutton =
       XCreateSimpleWindow(display, window, WINDOWWIDTH - BUTTON_WIDTH - PLAY_MARGIN,
                           WINDOWHEIGHT - PLAY_MARGIN - BOTTOM_MARGIN -
                               i * BUTTON_HEIGHT - (i - 1) * BUTTON_SEP,
                           BUTTON_WIDTH, BUTTON_HEIGHT, BUTTON_BORDER, black, darkcolor);
   i++;
   cleanbutton =
       XCreateSimpleWindow(display, window, WINDOWWIDTH - BUTTON_WIDTH - PLAY_MARGIN,
                           WINDOWHEIGHT - PLAY_MARGIN - BOTTOM_MARGIN -
                               i * BUTTON_HEIGHT - (i - 1) * BUTTON_SEP,
                           BUTTON_WIDTH, BUTTON_HEIGHT, BUTTON_BORDER, black, darkcolor);
   i++;
   seedbutton =
       XCreateSimpleWindow(display, window, WINDOWWIDTH - BUTTON_WIDTH - PLAY_MARGIN,
                           WINDOWHEIGHT - PLAY_MARGIN - BOTTOM_MARGIN -
                               i * BUTTON_HEIGHT - (i - 1) * BUTTON_SEP,
                           BUTTON_WIDTH, BUTTON_HEIGHT, BUTTON_BORDER, black, darkcolor);
   i++;
   fissionbutton =
       XCreateSimpleWindow(display, window, WINDOWWIDTH - BUTTON_WIDTH - PLAY_MARGIN,
                           WINDOWHEIGHT - PLAY_MARGIN - BOTTOM_MARGIN -
                               i * BUTTON_HEIGHT - (i - 1) * BUTTON_SEP,
                           BUTTON_WIDTH, BUTTON_HEIGHT, BUTTON_BORDER, black, darkcolor);
   i++;
   tempbutton =
       XCreateSimpleWindow(display, window, WINDOWWIDTH - BUTTON_WIDTH - PLAY_MARGIN,
                           WINDOWHEIGHT - PLAY_MARGIN - BOTTOM_MARGIN -
                               i * BUTTON_HEIGHT - (i - 1) * BUTTON_SEP,
                           BUTTON_WIDTH, BUTTON_HEIGHT, BUTTON_BORDER, black, darkcolor);
   i++;
   flakebutton =
       XCreateSimpleWindow(display, window, WINDOWWIDTH - BUTTON_WIDTH - PLAY_MARGIN,
                           WINDOWHEIGHT - PLAY_MARGIN - BOTTOM_MARGIN -
                               i * BUTTON_HEIGHT - (i - 1) * BUTTON_SEP,
                           BUTTON_WIDTH, BUTTON_HEIGHT, BUTTON_BORDER, black, darkcolor);
   i++;
   samplebutton =
       XCreateSimpleWindow(display, window, WINDOWWIDTH - BUTTON_WIDTH - PLAY_MARGIN,
                           WINDOWHEIGHT - PLAY_MARGIN - BOTTOM_MARGIN -
                               i * BUTTON_HEIGHT - (i - 1) * BUTTON_SEP,
                           BUTTON_WIDTH, BUTTON_HEIGHT, BUTTON_BORDER, black, darkcolor);
   i++;
   exportbutton =
       XCreateSimpleWindow(display, window, WINDOWWIDTH - BUTTON_WIDTH - PLAY_MARGIN,
                           WINDOWHEIGHT - PLAY_MARGIN - BOTTOM_MARGIN -
                               i * BUTTON_HEIGHT - (i - 1) * BUTTON_SEP,
                           BUTTON_WIDTH, BUTTON_HEIGHT, BUTTON_BORDER, black, darkcolor);
   i++;
   pausebutton =
       XCreateSimpleWindow(display, window, WINDOWWIDTH - BUTTON_WIDTH - PLAY_MARGIN,
                           WINDOWHEIGHT - PLAY_MARGIN - BOTTOM_MARGIN -
                               i * BUTTON_HEIGHT - (i - 1) * BUTTON_SEP,
                           BUTTON_WIDTH, BUTTON_HEIGHT, BUTTON_BORDER, black, darkcolor);
   i++;
   colorbutton =
       XCreateSimpleWindow(display, window, WINDOWWIDTH - BUTTON_WIDTH - PLAY_MARGIN,
                           WINDOWHEIGHT - PLAY_MARGIN - BOTTOM_MARGIN -
                               i * BUTTON_HEIGHT - (i - 1) * BUTTON_SEP,
                           BUTTON_WIDTH, BUTTON_HEIGHT, BUTTON_BORDER, black, darkcolor);
   i++;
   sidebutton =
       XCreateSimpleWindow(display, window, WINDOWWIDTH - BUTTON_WIDTH - PLAY_MARGIN,
                           WINDOWHEIGHT - PLAY_MARGIN - BOTTOM_MARGIN -
                               i * BUTTON_HEIGHT - (i - 1) * BUTTON_SEP,
                           BUTTON_WIDTH, BUTTON_HEIGHT, BUTTON_BORDER, black, darkcolor);
   i++;
   playground = XCreateSimpleWindow(display, window, PLAY_MARGIN,
                                    2 * PLAY_MARGIN + NUM_INFO_LINES * font_height,
                                    block * NCOLS, block * NROWS, 2, translate[4], white);

   /* pick the events to look for */
   event_mask = ExposureMask | ButtonPressMask | StructureNotifyMask;
   XSelectInput(display, window, event_mask);
   event_mask = ButtonPressMask;
   /* note that with this simple mask if one just covers a button
      it will not get redrawn.  I wonder if anyone will notice?  If I put
      the exposuremask in here, things flash irritatingly on being uncovered. */
   XSelectInput(display, quitbutton, event_mask);
   XSelectInput(display, pausebutton, event_mask);
   XSelectInput(display, restartbutton, event_mask);
   XSelectInput(display, cleanbutton, event_mask);
   XSelectInput(display, colorbutton, event_mask);
   XSelectInput(display, sidebutton, event_mask);
   XSelectInput(display, flakebutton, event_mask);
   XSelectInput(display, seedbutton, event_mask);
   XSelectInput(display, fissionbutton, event_mask);
   XSelectInput(display, tempbutton, event_mask);
   XSelectInput(display, samplebutton, event_mask);
   XSelectInput(display, exportbutton, event_mask);
   event_mask =
       ButtonReleaseMask | ButtonPressMask | PointerMotionHintMask | ButtonMotionMask;
   XSelectInput(display, playground, event_mask);

   /* make graphics contexts:
      gc for black on white (actually, lightcolor)
      gccolor for background and buttons
      gcr for reverse video  */

   gc = XCreateGC(display, window, 0, NULL);
   XSetFont(display, gc, font->fid);
   XSetForeground(display, gc, black);
   XSetBackground(display, gc, lightcolor);
   /* speed up? */
   XSetPlaneMask(display, gc,
                 black | white | translate[0] | translate[1] | translate[2] |
                     translate[3] | translate[16] | translate[4]);

   gcr = XCreateGC(display, window, 0, NULL);
   XSetFont(display, gcr, font->fid);
   XSetForeground(display, gcr, lightcolor);
   XSetBackground(display, gcr, darkcolor);

   gccolor = XCreateGC(display, window, 0, NULL);
   XSetFont(display, gccolor, font->fid);
   XSetForeground(display, gccolor, darkcolor);
   XSetBackground(display, gccolor, lightcolor);

   /* show the window and buttons */
   XMapWindow(display, window);
   XMapWindow(display, quitbutton);
   XMapWindow(display, pausebutton);
   XMapWindow(display, restartbutton);
   XMapWindow(display, cleanbutton);
   XMapWindow(display, colorbutton);
   if (block > 4)
      XMapWindow(display, sidebutton);
   if (!(fparam->N == 1 && fparam->next_param == NULL))
      XMapWindow(display, flakebutton);
   XMapWindow(display, seedbutton);
   XMapWindow(display, fissionbutton);
   XMapWindow(display, tempbutton);
   XMapWindow(display, samplebutton);
   XMapWindow(display, exportbutton);
   XMapWindow(display, playground);

   /* make image structure */
   /* wait for playground to be displayed before proceeding */
   i = 1; /* a flag */
   while (i) {
      XNextEvent(display, &report);
      switch (report.type) {
      case Expose:
         if (report.xexpose.window != playground)
            i = 0;
      default:
         break;
      }
   }

   /* FIXME: I Don't know why 6 is the right value here!? 8 doesn't work */
   buffer = malloc(block * block * NROWS * NCOLS * depth / 6);
   spinimage = XCreateImage(display, CopyFromParent, depth, ZPixmap, 0, buffer,
                            block * NCOLS, block * NROWS, 32, 0);
   XInitImage(spinimage);

   if (NULL == spinimage) {
      fprintf(stderr, "trouble creating image structure\n");
      exit(-1);
   }

   /* make sure everything get written first time */
   for (i = 0; i < block * NROWS; i++)
      for (j = 0; j < block * NCOLS; j++)
         XPutPixel(spinimage, i, j, translate[0]);
}

void cleanup() {
   XUnloadFont(display, font->fid);
   XFreeGC(display, gc);
   XFreeGC(display, gcr);
   XFreeGC(display, gccolor);
   XCloseDisplay(display);
   XDestroyImage(spinimage);
   exit(1);
}

void identify(tube_params *params, tube *tp, flake *fp, int x, int y) {
   int i, j, t;
   int size = (1 << fp->P);
   i = MIN(MAX(0, y - 2), size - 1);
   j = MIN(MAX(0, x - 2), size - 1);
   t = fp->Cell(i, j);
   snprintf(stringbuffer, BUFFER_LEN,
            "([DX] = %g uM, T = %5.3f C, 5-mer s.e.)    "
            "%s #%d = {%d %d %d %d} at (%d,%d)              ",
            1000000.0 * 20.0 * exp(-tp->Gmc), 4000 / (tp->Gse / 5 + 11) - 273.15,
            (params->tile_names[t] != NULL) ? params->tile_names[t] : "tile", t,
            tp->tileb[t][0], tp->tileb[t][1], tp->tileb[t][2], tp->tileb[t][3], i, j);
   XDrawImageString(display, window, gc, 5, 3 * font_height, stringbuffer,
                    strlen(stringbuffer));
}

void set_params(tube *tp, tube_params *params) {
   int i, j, n;

   /* make our own copy of tile set and strengths, so the calling program
   can do as it pleases with it's tileb & strength information */
   for (i = 0; i <= tp->N; i++)
      for (j = 0; j < 4; j++)
         tp->tileb[i][j] = params->tileb[i][j];
   for (i = 0; i <= tp->num_bindings; i++)
      tp->strength[i] = params->strength[i];
   for (i = 0; i <= tp->num_bindings; i++) {
      for (j = 0; j <= tp->num_bindings; j++) {
         tp->glue[i][j] = params->glue[i][j];
      }
   }

   tp->blast_rate_alpha = params->blast_rate_alpha;
   tp->blast_rate_beta = params->blast_rate_beta;
   tp->blast_rate_gamma = params->blast_rate_gamma;
   tp->blast_rate = params->blast_rate;
   tp->min_strength = params->min_strength;
   tp->wander = params->wander;
   tp->periodic = params->periodic;
   tp->fission_allowed = params->fission_allowed;
   tp->zero_bonds_allowed = params->zero_bonds_allowed;
   tp->present_list = params->present_list;
   tp->present_list_len = params->present_list_len;
   tp->untiltiles = params->untiltiles;
   tp->untiltilescount = params->untiltilescount;
   tp->hydro = params->hydro;
   tp->T = params->T;
   tp->k = params->k;
   tp->kas = exp(-params->Gas);
   tp->kao = exp(-params->Gao);
   tp->kam = exp(-params->Gam);
   tp->kae = exp(-params->Gae);
   tp->kah = exp(-params->Gah);
   tp->Gmc = params->Gmc;
   tp->anneal_g = params->anneal_g;
   tp->anneal_t = params->anneal_t;
   tp->anneal_h = params->anneal_h;
   tp->anneal_s = params->anneal_s;
   tp->startC = params->startC;
   tp->endC = params->endC;
   tp->seconds_per_C = params->seconds_per_C;
   tp->tinybox = params->tinybox;
   tp->default_seed_i = params->seed_i;
   tp->default_seed_j = params->seed_j;
   tp->initial_Gfc = params->Gfc;
   tp->Gmch = params->Gmch;

   tp->currentC = tp->startC;
   tp->Gse = params->Gse;
   tp->update_freq = params->updates_per_RC;

   tp->dt_right = params->dt_right;
   tp->dt_left = params->dt_left;

   tp->dt_up = params->dt_up;
   tp->dt_down = params->dt_down;

   /* set tp->Gcb from Ghyd */
   for (n = 0; n <= tp->N; n++)
      tp->Gcb[n] = 0;
   if (tp->hydro)
      for (n = tp->N / 2 + 1; n <= tp->N; n++)
         tp->Gcb[n] = params->Ghyd;

   /* set tp->conc from Gmc... and from Gmch for hydrolysis rules */
   tp->conc[0] = 0;
   for (n = 1; n <= tp->N; n++)
      tp->conc[0] +=
          (tp->conc[n] =
               exp(-((n > tp->N / 2 && tp->hydro) ? params->Gmch : params->Gmc)) *
               params->stoic[n]);
}

int main(int argc, char **argv) {
   int x, y, b, i, j;
   int clear_x = 0, clear_y = 0;
   int mousing = 0;
   int stat = 0;
   double new_Gse, new_Gmc;
   XEvent report;
   progname = argv[0];

   tube_params *params = (tube_params *)malloc(sizeof(tube_params));
   set_default_params(params);

   for (i = 0; i < MAXTILETYPES; i++) {
      translate[i] = 0; /* translate[] from tiles to colors */
   }

   getargs(argc, argv, params);

   if (import_flake_size > params->size) {
      fprintf(stderr, "Error: The input flake size is larger than the "
                      "simulation flake size.\n");
      exit(-1);
   }

   if (params->hydro) { /* automatically double the number of tiles */
      params->tileb = (int **)realloc(params->tileb, sizeof(int *) * (2 * N + 1));
      for (i = 1; i <= N; i++) {
         params->tileb[i + N] = (int *)calloc(sizeof(int), 4);
         for (j = 0; j < 4; j++)
            params->tileb[i + N][j] = params->tileb[i][j];
      }
      params->stoic = (double *)realloc(params->stoic, sizeof(double) * (2 * N + 1));
      for (i = 1; i <= N; i++)
         params->stoic[i + N] = params->stoic[i];
      N = N * 2;
   }

   if (linear) {
      linear_simulate(params->k, params->Gmc, params->Gse, tmax, emax, smax, mmax);
      return 0;
   }

   tube *tp;
   flake *fp;

   if (testing) {
      tp = init_tube(params->size_P, N, num_bindings);
      set_params(tp, params);
      setup_tube(tp);
#ifdef TESTING_OK
      run_xgrow_tests(tp, Gmc, Gse, seed_i, seed_j, seed_n, size);
#endif
      return 0;
   }

   if (XXX)
      openwindow(argc, argv);

   /* fprintf(stderr, "xgrow: tile set read, beginning simulation\n"); */

   /* set initial state */
   tp = init_tube(params->size_P, N, num_bindings);
   set_params(tp, params);
   setup_tube(tp);

   fprm = fparam;
   int size = (1 << tp->P);
   /* initialize flakes */
   while (fprm != NULL) {
      int fn;
      for (fn = 1; fn <= fprm->N; fn++) {
         if (tp->dt_left[fprm->seed_n]) {
            fprm->seed_n = tp->dt_left[fprm->seed_n]; // FIXME: vdouble
            fprm->seed_j--;
         }

         fp = init_flake(params->size_P, N, fprm->seed_i, fprm->seed_j, fprm->seed_n,
                         fprm->Gfc, tp->present_list_len, tp->periodic);

         insert_flake(fp, tp);

         if (tp->dt_right[fprm->seed_n]) {
            change_cell(fp, params->seed_i, params->seed_j + 1,
                        tp->dt_right[fprm->seed_n]);
            fp->seed_is_double_tile = 1;
         }
         assert(!tp->dt_left[fprm->seed_n]);
         if (tp->dt_down[fprm->seed_n]) {
            change_cell(fp, params->seed_i + 1, params->seed_j,
                        tp->dt_down[fprm->seed_n]);
            fp->seed_is_vdouble_tile = 1;
         }
         assert(!tp->dt_up[fprm->seed_n]);

         if (fprm->import_from != NULL) {
            fprintf(
                stderr,
                "WARNING: In imported flakes, the seed position is chosen randomly.\n");
            import_flake(tp, fp, fprm->import_from, fn);
         }
      }
      fprm = fprm->next_param;
   }

   //   print_tree(tp->flake_tree,0,'*');

   if (stripe_args != NULL) {
      /* STRIPE OPTION HAS HARDCODED TILESET NONSENSE -- A BUG */
      /* SHOULD REPLACE THIS OPTION BY READING INITIAL FLAKE FROM FILE */
      int i, j, k, w;
      double p;
      char *s = stripe_args;
      char XOR[2][2] = {{4, 7}, {6, 5}}; /* XOR[S][E] */
      char c, cc;
      i = params->size - 1;
      j = atoi(s) % params->size;
      for (k = 0; k < params->size; k++)
         change_cell(fp, (i - k + params->size) % params->size, (j + k) % params->size,
                     4 + random() % 4);
      fp->seed_i = i;
      fp->seed_j = j;
      fp->seed_n = fp->Cell(i, j);
      s = strchr(s, ':');
      while (s != NULL) {
         p = atof(s + 1);
         s = strchr(s, ',');
         if (s != NULL) {
            w = atoi(s + 1);
            s = strchr(s, ':');
            for (; w > 0; w--) {
               i = (i - 1 + params->size) % params->size;
               for (k = 0; k < params->size; k++) {
                  cc = c = XOR[(fp->Cell((i - k + 1 + params->size) % params->size,
                                         (j + k) % params->size) -
                                4) /
                               2][(fp->Cell((i - k + params->size) % params->size,
                                            (j + k + 1) % params->size) -
                                   4) /
                                  2];
                  if (drand48() < p)
                     do
                        cc = 4 + random() % 4;
                     while (cc == c);
                  change_cell(fp, (i - k + params->size) % params->size,
                              (j + k) % params->size, cc);
                  tp->events--; /* don't count these as events */
                                /* ERROR: stats are also modified!!! */
               }
            }
         }
      }
      /* no corner or boundary tiles for stripe simulations */
      tp->conc[0] -= tp->conc[1];
      tp->conc[0] += (tp->conc[1] = exp(-35));
      tp->conc[0] -= tp->conc[2];
      tp->conc[0] += (tp->conc[2] = exp(-35));
      tp->conc[0] -= tp->conc[3];
      tp->conc[0] += (tp->conc[3] = exp(-35));
   }

   // fprintf(stderr, "flake initialized, size_P=%d, size=%d\n",size_P,size);

   new_Gse = tp->Gse;
   new_Gmc = tp->Gmc;
   if (tracefp != NULL)
      write_datalines(params, tp, fp, tracefp, "\n");

   // fprintf(stderr, "tmax=%f  emax=%d  smax=%d\n",tmax,emax,smax);

   if (XXX)
      repaint(params, tp, fp);

   /* loop forever, looking for events */
   while ((tmax == 0 || tp->t < tmax) && (emax == 0 || tp->events < emax) &&
          (smax == 0 || tp->stat_a - tp->stat_d < smax) &&
          (mmax == 0 || tp->stat_m < mmax) &&
          (smin == -1 || tp->stat_a - tp->stat_d > smin) &&
          (fsmax == 0 || tp->largest_flake_size < fsmax) &&
          (tp->seconds_per_C == 0 || tp->currentC > tp->endC) &&
          !(tp->untiltiles && tp->all_present)) {

      if (!XXX) {
         simulate(tp, params->update_rate, tmax, emax, smax, fsmax, smin, mmax);
         if (tracefp != NULL)
            write_datalines(params, tp, fp, tracefp, "\n");
         if (export_mode == 2 && export_movie == 1)
            export_flake(params, tp, "movie", fp);
      } else {
         if (0 == paused && 0 == mousing && !XPending(display)) {
            simulate(tp, params->update_rate, tmax, emax, smax, fsmax, smin, mmax);
            fp = tp->flake_list;
            assert(!fp || !tp->tinybox ||
                   ((!fp->seed_is_double_tile && fp->tiles > 1) || fp->tiles > 2));
            if (tracefp != NULL)
               write_datalines(params, tp, fp, tracefp, "\n");
            if (export_mode == 2 && export_movie == 1)
               export_flake(params, tp, "movie", fp);
            if (fp && fp->flake_conc > 0)
               recalc_G(fp);
            // make sure displayed G is accurate for conc's
            // hopefully this won't slow things down too much.
            stat++;
            if (stat == 1) {
               stat = 0;
               repaint(params, tp, fp);
            }
         }
         if (paused | mousing | XPending(display)) {
            XNextEvent(display, &report);
            switch (report.type) {
            case Expose:
               if (report.xexpose.count != 0)
                  break; /* more in queue, wait for them */
               repaint(params, tp, fp);
               break;
            case ConfigureNotify:
               break;
            case MotionNotify:
               if (report.xbutton.window == playground) {
                  int newx, newy;
                  Window root, child;
                  unsigned int keys_buttons;
                  int window_x, window_y;
                  if (tp->hydro)
                     break; /* don't know how to reset params */
                  x = report.xbutton.x / block;
                  y = report.xbutton.y / block;
                  b = report.xbutton.button;
                  if (mousing == 0) { /* identify while moving around */
                     identify(params, tp, fp, x, y);
                  }
                  if (mousing == 3) {
                     /* was sketch(x,y,b); now change Gse & Gmc */
                     new_Gse = (30.0 * x) / params->size;
                     new_Gmc = 30 - (30.0 * y) / params->size;
                     /* draw current Gse, Gmc values */
                     snprintf(stringbuffer, BUFFER_LEN,
                              "Gmc=%4.1f->%4.1f  Gse=%4.1f->%4.1f", tp->Gmc, new_Gmc,
                              tp->Gse, new_Gse);
                     XDrawImageString(display, window, gc, 5, 3 * font_height,
                                      stringbuffer, strlen(stringbuffer));
                     /* this is necessary in order to get more than one MotionNotify
                   event; I don't know why. */
                  }
                  if (!XQueryPointer(display, playground, &root, &child, &window_x,
                                     &window_y, &newx, &newy, &keys_buttons)) {
                     mousing = 0;
                     fprintf(stderr, "Weird X feature\n");
                  }
               }
               break;
            case ButtonRelease:
               if (mousing == 3) { // change Gmc, Gse w/ visual GUI
                  fprintf(stderr, "Changing Gmc -> %f, Gse -> %f \n", new_Gmc, new_Gse);
                  reset_params(tp, tp->Gmc, tp->Gse, new_Gmc, new_Gse, params->Gseh);
                  if (params->Gfc > 0)
                     params->Gfc += (new_Gmc - tp->Gmc);
                  fprm = fparam;
                  while (fprm != NULL) { /* fix up all flake info, for restarting */
                     if (fprm->Gfc > 0)
                        fprm->Gfc += (new_Gmc - tp->Gmc);
                     fprm = fprm->next_param;
                  }
                  tp->Gse = new_Gse;
                  tp->Gmc = new_Gmc;
                  showpic(tp, fp, errorc);
               } else if (mousing == 1) {
                  /* clear a region, i.e., "puncture" */
                  int i, j, mi, Mi, mj, Mj;
                  x = report.xbutton.x / block;
                  y = report.xbutton.y / block;
                  b = report.xbutton.button;
                  if (clear_y > y) {
                     mi = MIN(MAX(0, y - 2), params->size - 1);
                     Mi = MIN(MAX(0, clear_y - 2), params->size - 1);
                  } else {
                     mi = MIN(MAX(0, clear_y - 2), params->size - 1);
                     Mi = MIN(MAX(0, y - 2), params->size - 1);
                  }
                  if (clear_x > x) {
                     mj = MIN(MAX(0, x - 2), params->size - 1);
                     Mj = MIN(MAX(0, clear_x - 2), params->size - 1);
                  } else {
                     mj = MIN(MAX(0, clear_x - 2), params->size - 1);
                     Mj = MIN(MAX(0, x - 2), params->size - 1);
                  }
                  if (mi != Mi && mj != Mj) {
                     int fa = tp->fission_allowed;
                     if (mi <= fp->seed_i && fp->seed_i <= Mi && mj <= fp->seed_j &&
                         fp->seed_j <= Mj) {
                        int si = fp->seed_i, sj = fp->seed_j;
                        // better move the seed out of the way, if possible
                        for (i = mi; i <= Mi; i++)
                           if (fp->Cell(i, mj - 1) > 0) {
                              si = i;
                              sj = (mj - 1 + params->size) % params->size;
                           }
                        for (i = mi; i <= Mi; i++)
                           if (fp->Cell(i, Mj + 1) > 0) {
                              si = i;
                              sj = (Mj + 1) % params->size;
                           }
                        for (j = mj; j <= Mj; j++)
                           if (fp->Cell(mi - 1, j) > 0) {
                              si = (mi - 1 + params->size) % params->size;
                              sj = j;
                           }
                        for (j = mj; j <= Mj; j++)
                           if (fp->Cell(Mi + 1, j) > 0) {
                              si = (Mi + 1) % params->size;
                              sj = j;
                           }
                        change_seed(fp, si, sj);
                     }
                     tp->fission_allowed = 1;
                     for (i = mi; i <= Mi; i++)
                        for (j = mj; j <= Mj; j++) {
                           if (fp->Cell(i, j) > 0) {
                              change_cell(fp, i, j, 0);
                              flake_fission(fp, i, j);
                           }
                        }
                     tp->fission_allowed = fa;
                  }
                  showpic(tp, fp, errorc);
               }
               mousing = 0;
               break;
            case ButtonPress:
               if (report.xbutton.window == quitbutton) {
                  closeargs(tp,
                            params); // This cleans up and exports things if instructed.
                  exit(0);
               } else if (report.xbutton.window == pausebutton) {
                  setpause(1 - paused);
                  repaint(params, tp, fp);
               } else if (report.xbutton.window == restartbutton) {
                  free_tube(tp);
                  tp = init_tube(params->size_P, N, num_bindings);
                  set_params(tp, params);
                  fprm = fparam;
                  while (fprm != NULL) {
                     int fn;
                     for (fn = 1; fn <= fprm->N; fn++) {
                        if (tp->dt_left[fprm->seed_n]) {
                           fprm->seed_n = tp->dt_left[fprm->seed_n]; // FIXME: vdouble
                           fprm->seed_j--;
                        }
                        insert_flake(fp =
                                         init_flake(params->size_P, N, fprm->seed_i,
                                                    fprm->seed_j, fprm->seed_n, fprm->Gfc,
                                                    tp->present_list_len, tp->periodic),
                                     tp);
                        if (tp->dt_right[fprm->seed_n]) {
                           change_cell(fp, params->seed_i, params->seed_j + 1,
                                       tp->dt_right[fprm->seed_n]);
                           fp->seed_is_double_tile = 1;
                        }
                        assert(!tp->dt_left[fprm->seed_n]);

                        if (fprm->import_from != NULL) {
                           fprintf(stderr, "WARNING: In imported flakes, the seed "
                                           "position is chosen randomly.\n");
                           import_flake(tp, fp, fprm->import_from, fn);
                        }
                     }
                     fprm = fprm->next_param;
                  }
                  repaint(params, tp, fp);
               } else if (report.xbutton.window ==
                          colorbutton) { // show tiles or error or hyd
                  settilecolor(tp, tp->hydro ? (errorc + 1) % 3 : (errorc + 1) % 2);
                  repaint(params, tp, fp);
               } else if (report.xbutton.window ==
                          sidebutton) { // show box or null/weak/strong
                  setsidecolor((errors + 1) % 2);
                  repaint(params, tp, fp);
               } else if (report.xbutton.window == seedbutton) {
                  setwander(tp, 1 - tp->wander);
                  repaint(params, tp, fp);
               } else if (report.xbutton.window == fissionbutton) {
                  tp->fission_allowed = (tp->fission_allowed + 1) % 3;
                  setfission(tp, tp->fission_allowed);
                  update_all_rates(tp);
                  repaint(params, tp, fp);
               } else if (report.xbutton.window == cleanbutton) {
                  if (x < 50) { // do a clean_flake cycle
                     fprintf(stderr, "Cleaning 1 cycle, clean_X=%f\n", clean_X);
                     clean_flake(fp, clean_X, 1);
                  } else if (x < 90) { // do a fill_flake cycle
                     fprintf(stderr, "Filling 1 cycle, fill_X=%f\n", fill_X);
                     fill_flake(fp, fill_X, 1);
                  } else {
                     fprintf(stderr, "Repairing, uniqueness for T=%f\n", repair_unique_T);
                     repair_flake(fp, repair_unique_T, tp->Gse);
                  }
                  repaint(params, tp, fp);
               } else if (report.xbutton.window == exportbutton) {
                  if (x > 70) { // change from ONE to ALL to MOVIE
                     setexport((export_mode + 1) % 3);
                  } else if (export_mode == 0) {
                     // output to file (unless MOVIE mode already)
                     export_flake(params, tp, "flake", fp);
                  } else if (export_mode == 1) {
                     flake *tfp = tp->flake_list;
                     while (tfp->next_flake != NULL) {
                        export_flake(params, tp, "flake", tfp);
                        tfp = tfp->next_flake;
                     }
                  } else if (export_mode == 2) {
                     // turn movie mode on & off
                     export_movie = !export_movie;
                  }
                  repaint(params, tp, fp);
               } else if (report.xbutton.window == samplebutton) {
                  flake *tfp;
                  int n, num_big = 0;
                  // stop simulation if you're sampling.
                  setpause(1);
                  sampling = 1;
                  // pick a (if possible non-monomer) sample and add to the field
                  for (tfp = tp->flake_list; tfp != NULL; tfp = tfp->next_flake) {
                     if (tfp->tiles > 1)
                        num_big++;
                  }
                  n = random() % (num_big > 0 ? num_big : tp->num_flakes);
                  tfp = tp->flake_list;
                  while (num_big > 0 && tfp->tiles == 1) {
                     tfp = tfp->next_flake;
                  }
                  for (i = 0; i < n; i++) {
                     tfp = tfp->next_flake;
                     while (num_big > 0 && tfp->tiles == 1) {
                        tfp = tfp->next_flake;
                     }
                  }
                  add_sample_pic(tp, tfp, errorc);
                  repaint(params, tp, fp);
               } else if (report.xbutton.window == flakebutton) {
                  flake *tfp = tp->flake_list;
                  // cycle through flakes.
                  if (x > 80) {
                     fp = fp->next_flake;
                     if (fp == NULL)
                        fp = tp->flake_list;
                  } else if (x < 40) {
                     while (tfp && tfp->next_flake != NULL && tfp->next_flake != fp)
                        tfp = tfp->next_flake;
                     if (tfp == NULL)
                        fp = tp->flake_list;
                     else
                        fp = tfp;
                  } else { // find the biggest flake.
                     fp = tfp;
                     if (fp) {
                        while (tfp->next_flake != NULL) {
                           tfp = tfp->next_flake;
                           if (tfp->tiles > fp->tiles)
                              fp = tfp;
                        }
                     }
                  }
                  //             print_tree(tp->flake_tree,0,'*');
                  sampling = 0;
                  repaint(params, tp, fp);
               } else if (report.xbutton.window == tempbutton) { // change Gse w/ button
                  if (tp->hydro)
                     break; /* don't know how to reset params */
                  if (x > 60)
                     new_Gse = tp->Gse - 0.1;
                  else
                     new_Gse = tp->Gse + 0.1;
                  reset_params(tp, tp->Gmc, tp->Gse, new_Gmc, new_Gse, params->Gseh);
                  tp->Gse = new_Gse;
                  repaint(params, tp, fp);
               } else if (report.xbutton.window == playground) // we're in ButtonPress
               {
                  if (b == 3) {
                     if (tp->hydro)
                        break; /* don't know how to reset params */
                     new_Gse = (30.0 * x) / params->size;
                     new_Gmc = 30 - (30.0 * y) / params->size;
                     /* draw current Gse, Gmc values */
                     snprintf(stringbuffer, BUFFER_LEN,
                              "Gmc=%4.1f->%4.1f  Gse=%4.1f->%4.1f", tp->Gmc, new_Gmc,
                              tp->Gse, new_Gse);
                     XDrawImageString(display, window, gc, 5, 3 * font_height,
                                      stringbuffer, strlen(stringbuffer));
                     /* later: if down button, draw T=1/T=2 diagram */
                     mousing = 3;
                     showphase(params);
                  } else if (b == 2) { // "puncture"
                     mousing = 1;
                     clear_x = x;
                     clear_y = y;      /* prepare to clear a region */
                  } else if (b == 1) { // "identify"
                     identify(params, tp, fp, x, y);
                  }
               } else {
                  /*  what here ?  */
               }
               break;
            default:
               break;
            } /* end of switch */
         }    /* end of if XPending */
      }
   } /* end of while(...) if...else */

   closeargs(tp, params);
   return 0;
} /* end of main */
