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
   But in cases of big erroneous holes, it can (a) result in many mismatches, since non-unique
   filling proceeds somewhat arbitrarily to completion, and (b) result in a non-connected
   flake, which is weird if the simulation is continued. (So far, non-connected flakes always 
	 get cleaned up well, but I'd still call this a BUG.)   Also, since only interior holes
   get filled in, any error that stalls growth or leaves an inlet will not get counted.
   Also added command-line "pause", primarily for examining saved files w/o starting simulator.
   2/25/04 Fixed bug in locally_fission_proof() call when chunk_fission is set. 
(Sometimes bogus value resulted in non-fission when fission should have occurred.)
   6/4/04  Bond types can now be given as names.  See spiral.tiles for example.
   6/20/04 Added command-line option for setting tile stoichiometry.
   Made sure file is flushed immediately after export.
   3/8/05  Fixed bug in chunk_fission that caused bogus "impossible chunk chosen" errors (RS)

   3/8/05 Implemented a "doubletile" function in configuration files
   that allows two tile types to be considered a connected,
   horizontal rectangle for purposes of tile additions and removals.  This works as follows:
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
   3/19/05 Adding blast damage for demonstrating self-repair.  alpha, beta, gamma command-line args.
   Wander allows the original seed tile to be lost.
   It's a bit strange when no_fission is used, because fission is tested
   sequentially as tiles in the square are being removed, rather than
   after the whole square has been removed. 

   5/14/05 Added an anneal feature, at which Gse starts at some value and
   declines by a time constant until it reaches its final value.  It shouldn't be hard to implement a melt feature, with the exact reverse behavior, but right now we require that Gse initial be higher than Gse final.  Note that it would be horribly inefficient to implement an exact annealing feature of this type, because all rates would have to be updated after every event.  As a compromise we update Gse (and therefore the rates) some fixed number of (for now 100) times per time constant. -- RS
   7/30/05 Added testing to xgrow -- new option "testing" will test whether a tile set obeys detailed balance once it reaches equilibrium.  Takes a while to do, and effectively forever if you choose a tile set that takes a looooong time to reach equilibrium, like the Sierpinski set. -- RS
   11/20/05 Fixed a bug in count_flakes().  It's still very brittle.  -EW
   12/12/05 Fixed up the handling of Gse and tp->Gse, which was leading to the GUI display being occasionally off. -EW
   12/15/05 Added flake loading option to offset the seed assembly position -EW

   .... various undocumented improvements by Erik Winfree (EW), Rebecca Schulman (RS), and Constantine Evans (CE) ....

2008 Changed tile and bond data type to short int, for > 256 types.  For smaller memory footprint, xgrow-small is avalailable. (CE)

   3/1/09 Fixed aTAM bug where it hangs if there's no possible move.  (EW & CE)

   TO DO List:

   * If the tile set specifies a stoichiometry of 0 (e.g. for the seed), the simulation can freak out.
   * export ALL may actually be ALL-BUT-ONE.  Test and fix, please?
   * add a DILUTE/CONCENTRATE button.
   * files should be flushed after each use.  arrayfile should *append*, not overwrite, properly
   adding to flake{} array at the appropriate index?
   * option for chunk_fission, but disallowing larger flake_fission events, would be desirable, so
   as to have a simulation in detailed balance.
   * for the same reason, we'd need a chunk_fusion option -- calculates equilibrium concentrations
   for all dimers & 2x2 blocks, and adds these as on-event options.  total on-rate would include
   these; in circumstances where they can't occur, result would just be a non-event.  There are 
   so many 2x2 blocks, we'd need a binary tree to efficiently select one. (Not a bad idea for
	 large tile sets, even w/o chunks.)  To avoid non-events, we could have 4 "empty" arrays, one
   for each size.  This would be a significant change to the code.
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
   * for fission, need option that chooses largest component, not component connected to the seed
   * program should die when user closes the window, rather than only when the user presses "quit"
   * when requested size is larger than the screen, the program should quit rather than 
   crash.  (See the code near "BadMatch", search for it below.)

   *  The (undocumented) hydrolysis mode is now broken; it crashes.  This model is a 2D generalization of microtubule / actin hydrolysis that leads to treadmilling and dynamic instability; in 1999 I was curious if it could be utilized to decrease error rates.
   *  Xgrow still sometimes hangs (or doesn't allow you to quit) when there is no available move, but it's not obvious when this happens.  Example: "xgrow barseed T=2 block=2 seed=180,180,1"  Here, tile 5 should be the seed.
   *  With multiple flakes, if the simulation is running, only the last flake is shown.  To see other flakes, you have to pause and then use next/big/prev to find one of interest.  But then when you run the simulation again, it pops to the last flake.  I think one would like to keep looking at the same flake you were looking at before pressing RUN.  Example: "xgrow barcode block=2 size=128 Gmc=15 Gse=10 addflakes=130,130,1:100@24"


   Compiling:  run "make"

   */

# include <X11/Xlib.h>
# include <X11/Xutil.h>
# include <X11/Xos.h>
# include <X11/Xatom.h> 
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <sys/time.h>
# include <unistd.h>
# include <math.h>
# include <assert.h>
# include <limits.h>

# include "grow.h"
#ifdef TESTING_OK
# include "xgrow-tests.h"
#endif

   /* lattice dimensions (plus two for boundaries): */
   /* NCOLS should be a multiple of bytes per long word (BPW) */
   /* THIS is for X bitmap stuff, so add total 4 for boundary */
# define NBDY 2
   /* THESE BECOME VARIABLES, NOW THAT FIELD SIZE IS VARIABLE 
# define NROWS (256+NBDY*2)
# define NCOLS (256+NBDY*2)
# define VOLUME (NROWS*NCOLS)

# define WINDOWWIDTH (block*NCOLS+PLAYLEFT+BOUNDWIDTH+20)
# define WINDOWHEIGHT (PLAYTOP+block*NROWS+10)
*/
# define PLAYTOP 120
# define PLAYLEFT 10
# define BOUNDHEIGHT 80
# define LATTICEHEIGHT 60
# define BOUNDWIDTH  140
# define BARHEIGHT 45
# define BARWIDTH  135
# define BARTOP 55
# define BARLEFT 4
# define XHEIGHT 40
# define NSIZE 42
#ifdef SMALL
# define MAXTILETYPES 256
#else
# define MAXTILETYPES USHRT_MAX
#endif

#define IS_ARG_MATCH(arg,s) (strncmp(arg,s,strlen(s))==0)
   
   long int translate[MAXTILETYPES]; /* for converting colors */
   int paused=0, errorc=0, errors=0, sampling=0;
   int export_mode=0, export_flake_n=1, export_movie_n=1, export_movie=0; 
   FILE *export_fp=NULL;
   int update_rate=10000;
   static char *progname;
   char stringbuffer[256];
   char tileset_name[256];
   int testing = 0;
   int initial_rc = 1;
   int *present_list=NULL;
   int present_list_len=0;
   int untiltiles=0,untiltilescount=0;


   /* various window stuff */
   Display *display;
   int screen;
   Window window,quitbutton,pausebutton,playground,
   restartbutton,colorbutton,sidebutton,flakebutton,seedbutton,tempbutton, 
   fissionbutton, samplebutton, exportbutton,cleanbutton;
GC gc, gcr, gccolor;
XImage *spinimage=NULL;
XFontStruct *font=NULL;
int font_height;
XSizeHints size_hints;
long event_mask;
int depth;
int darkcolor,lightcolor,black,white,
    errorcolor,goodcolor,hydrocolor,hydroerrcolor,strongcolor,weakcolor,nullcolor;
char *tile_colors[MAXTILETYPES]={ "black",
   "blue",      "red",      "green",      "yellow", "gold",   "purple", "white", 
   "dark blue", "dark red", "dark green", "wheat",  "orange", "cyan",   "light grey"};


/* simulation parameters */
tube *tp; flake *fp; 
double Gse, Gmc, ratek, T;
double Gseh, Gmch, Ghyd, Gas, Gam, Gae, Gah, Gao, Gfc;

double blast_rate_alpha=0;
double blast_rate_beta=4;  // k>3 required for finite rate of blasting a given tile in infinite size flakes (gamma=0)
double blast_rate_gamma=0;
double blast_rate=0;
int first_tile, second_tile;
int double_tiles=0;
int vdouble_tiles=0;
double min_strength=1;

int NROWS,NCOLS,VOLUME,WINDOWWIDTH,WINDOWHEIGHT;
int size=256, size_P=8; 
int block=1; /* default to small blocks; calling with argument changes this */
int wander, periodic, linear, fission_allowed, zero_bonds_allowed; 
FILE *tracefp, *datafp, *arrayfp, *tilefp, *largeflakefp;
FILE *untiltilescountfp;
int tthresh=0;
double *strength;
double **glue;
int *dt_right, *dt_left;
int *dt_down, *dt_up;
int N, num_bindings, tileb_length;
int **tileb; double *stoic;
int hydro; 
int clean_cycles=0; double clean_X=1.0; int fill_cycles=0; double fill_X=1.0; 
double error_radius=0.0; double repair_unique_T=2.0; int repair_unique=0;
double tmax; int emax, smax, fsmax, smin, mmax;
int seed_i,seed_j,seed_n;
double tinybox = 0;
double anneal_h, anneal_s, startC, endC, seconds_per_C = 0;
double anneal_g, anneal_t = 0;
char *stripe_args=NULL;
int XXX=1;  /* If 1, draw the simulation, otherwise don't.  (as inferred from the effect of -nw)*/
/* This has been adjusted so that updating the temperature and thus
   the required rates takes about 1% of the simulation time */
int updates_per_RC = 1000;

int import=0; /* Are we importing flakes? */
int import_flake_size = 0;
int import_offset_i=0, import_offset_j=0;

char newline[2];

struct flake_param {
   int seed_i,seed_j,seed_n,N; double Gfc; FILE *import_from;  struct flake_param *next_param;
};
struct flake_param *fparam=NULL, *fprm; 

int count_flakes(FILE *flake_file);

int parse_arg_line(char *arg)
{
   if (IS_ARG_MATCH(arg,"block=")) 
      block=MAX(1,MIN(30,atoi(&arg[6])));
   else if (IS_ARG_MATCH(arg,"size=")) 
      size=MAX(32,MIN(4096,atoi(&arg[5])));
   else if (IS_ARG_MATCH(arg,"rand=")) 
   { srand48(atoi(&arg[5])); srandom(atoi(&arg[5])); }
   else if (IS_ARG_MATCH(arg,"k=")) ratek=atof(&arg[2]);
   else if (IS_ARG_MATCH(arg,"Gmc=")) Gmc=atof(&arg[4]);
   else if (IS_ARG_MATCH(arg,"Gse=")) Gse=atof(&arg[4]);
   // FIXME FIXME FIXME: HYDRO OPTIONS ARE DISABLED
   //else if (IS_ARG_MATCH(arg,"Gmch=")) {hydro=1; Gmch=atof(&arg[5]);}
   //else if (IS_ARG_MATCH(arg,"Gseh=")) {hydro=1; Gseh=atof(&arg[5]);}
   //else if (IS_ARG_MATCH(arg,"Ghyd=")) {hydro=1; Ghyd=atof(&arg[5]);}
   else if (IS_ARG_MATCH(arg,"Gas=")) {hydro=1; Gas=atof(&arg[4]);}
   else if (IS_ARG_MATCH(arg,"Gam=")) {hydro=1; Gam=atof(&arg[4]);}
   else if (IS_ARG_MATCH(arg,"Gae=")) {hydro=1; Gae=atof(&arg[4]);}
   //else if (IS_ARG_MATCH(arg,"Gah=")) {hydro=1; Gah=atof(&arg[4]);}
   //else if (IS_ARG_MATCH(arg,"Gao=")) {hydro=1; Gao=atof(&arg[4]);}
   else if (IS_ARG_MATCH(arg,"Gfc=")) {Gfc=atof(&arg[4]);}
   else if (IS_ARG_MATCH(arg,"stoic=")) {
      double sts; int stn; char *sp;
      sp=strchr(&arg[6],'[');
      if (sp!=NULL) {
	 stn=atoi(&arg[6]); sts=atof(sp+1);
	 if (stn<=N && stn>0 && sts>0) {
	    stoic[stn]=sts;
	    fprintf(stderr,"stoic of tile %d = [%f]\n",stn,sts);
	 } else { 
	    if (sts==0) {
	       fprintf(stderr,"Stoic == 0.0 is not allowed; use a small but positive value.\n"); 
	       return -1;
	    }
	    fprintf(stderr,"Requested stoic of tile %d = [%f] is impossible!\n",stn,sts); 
	    return -1;
	 }
      }
   }
   else if (IS_ARG_MATCH(arg,"T=")) T=atof(&arg[2]);
   else if (IS_ARG_MATCH(arg,"pause=True")) paused=1;
   else if (IS_ARG_MATCH(arg,"pause=False")) paused=0;
   else if (IS_ARG_MATCH(arg,"pause")) paused=1;
   else if (IS_ARG_MATCH(arg,"periodic=True")) periodic=1;
   else if (IS_ARG_MATCH(arg,"periodic=False")) periodic=0;
   else if (IS_ARG_MATCH(arg,"periodic")) periodic=!periodic;
   else if (IS_ARG_MATCH(arg,"wander=True")) wander=1;
   else if (IS_ARG_MATCH(arg,"wander=False")) wander=0;
   else if (IS_ARG_MATCH(arg,"wander")) wander=!wander;
   else if (IS_ARG_MATCH(arg,"fission=off")) fission_allowed=0;
   else if (IS_ARG_MATCH(arg,"fission=on")) fission_allowed=1;
   else if (IS_ARG_MATCH(arg,"fission=chunk")) fission_allowed=2; 
   else if (IS_ARG_MATCH(arg,"no_fission")) fission_allowed=0;
   else if (IS_ARG_MATCH(arg,"fission")) fission_allowed=1;
   else if (IS_ARG_MATCH(arg,"chunk_fission")) fission_allowed=2;
   else if (IS_ARG_MATCH(arg,"blast_rate_alpha=")) blast_rate_alpha=atof(&arg[17]);
   else if (IS_ARG_MATCH(arg,"blast_rate_beta="))  blast_rate_beta =atof(&arg[16]);
   else if (IS_ARG_MATCH(arg,"blast_rate_gamma=")) blast_rate_gamma=atof(&arg[17]);
   else if (IS_ARG_MATCH(arg,"zero_bonds")) zero_bonds_allowed=1;
   else if (IS_ARG_MATCH(arg,"no_zero_bonds")) zero_bonds_allowed=0;
   else if (IS_ARG_MATCH(arg,"seed=")) {
      char *p=(&arg[5]);
      seed_i=atoi(p);
      if ((p=strchr(p,','))!=NULL) {
	 seed_j=atoi(p+1);
	 if ((p=strchr(p+1,','))!=NULL) seed_n=atoi(p+1);
      }
   }
   else if (IS_ARG_MATCH(arg,"addflakes=")) {
      char *p=(&arg[10]);
      fprm = (struct flake_param *)malloc(sizeof(struct flake_param)); 
      fprm->seed_i=fprm->seed_j=130; fprm->seed_n=1; fprm->N=1; fprm->Gfc=0;
      fprm->import_from=NULL;
      fprm->next_param=fparam;
      fparam=fprm;
      fparam->seed_i=atoi(p);
      if ((p=strchr(p,','))!=NULL) {
	 fparam->seed_j=atoi(p+1);
	 if ((p=strchr(p+1,','))!=NULL) {
	    fparam->seed_n=atoi(p+1);
	    if ((p=strchr(p+1,':'))!=NULL) {
	       fparam->N=atoi(p+1);
	       if ((p=strchr(p+1,'@'))!=NULL) {
		  fparam->Gfc=atof(p+1);
	       }
	    }
	 }
      }
   }
   else if (IS_ARG_MATCH(arg,"tinybox=")) {
      char *p=(&arg[8]);
      tinybox = atof(p);
   }
   else if (IS_ARG_MATCH(arg,"stripe=")) 
   { stripe_args=(&arg[7]); periodic=1; wander=1; }
   else if (strcmp(arg,"-nw")==0) XXX=0;
   else if (IS_ARG_MATCH(arg,"window=True")) XXX=1;
   else if (IS_ARG_MATCH(arg,"window=False")) XXX=0;
   else if (strcmp(arg,"-linear")==0) linear=1;
   else if (strncmp(arg,"linan=",6) == 0) {
      char *c;
      c = &arg[6];
      anneal_h = atof(c);
      c = strchr(c,',') + 1;
      if (c == NULL) {
	 fprintf(stderr,"Usage : linan=h,s,C0,Cfin,dt ");
	 return -1;
      }
      anneal_s = atof(c);
      c = strchr(c,',') + 1;
      if (c == NULL) {
	 fprintf(stderr,"Usage : linan=h,s,C0,Cfin,dt ");
	 return -1;
      }
      startC = atof(c);
      c = strchr(c,',') + 1;
      if (c == NULL) {
	 fprintf(stderr,"Usage : linan=h,s,C0,Cfin,dt ");
	 return -1;
      }
      endC = atof(c);
      c = strchr(c,',') + 1;
      if (c == NULL) {
	 fprintf(stderr,"Usage : linan=h,s,C0,Cfin,dt ");
	 return -1;
      }
      seconds_per_C = atof(c);
   }
   else if (strncmp(arg,"anneal=",7) == 0) {
      char *c;
      c = &arg[7];
      anneal_g = atof(c);
      c = strchr(c,',');
      if (c == NULL) {
	 fprintf(stderr,"Anneal option requires both an initial Gse and a time constant of annealing.\n");
	 return -1;
      }
      anneal_t = atof(c+1);
      if (anneal_t <= 0) {
	 fprintf(stderr,"Time constant for annealing must be positive.\n");
	 return -1;
      }
   }
   else if (IS_ARG_MATCH(arg,"doubletile=")) {
      char *p;
      first_tile = atoi(&arg[11]);
      p = strchr(&arg[12],',');
      if (p == NULL) {
	 fprintf(stderr,"Two tiles with a comma separating them are required for the double tile option.\n");
	 return -1;
      }
      second_tile = atoi(p + 1);
      if (tileb[first_tile][1] != tileb[second_tile][3]) {
	 fprintf(stderr,"Pieces of a double tile must have a matching bond in between them.\n");
	 return -1;
      }
      dt_right[first_tile] = second_tile;
      dt_left[second_tile] = first_tile;
      double_tiles = 1;
   }
   else if (IS_ARG_MATCH(arg,"vdoubletile=")) {
      char *p;
      first_tile = atoi(&arg[12]);
      p = strchr(&arg[13],',');
      if (p == NULL) {
	 fprintf(stderr,"Two tiles with a comma separating them are required for the vdouble tile option.\n");
	 return -1;
      }
      second_tile = atoi(p + 1);
      if (tileb[first_tile][2] != tileb[second_tile][0]) {
	 fprintf(stderr,"Pieces of a vdouble tile must have a matching bond in between them.\n");
	 return -1;
      }
      dt_down[first_tile] = second_tile;
      dt_up[second_tile] = first_tile;
      vdouble_tiles = 1;
   }
   else if (IS_ARG_MATCH(arg,"update_rate=")) 
      update_rate=MAX(1,MIN(atol(&arg[12]),10000000));
   else if (IS_ARG_MATCH(arg,"tracefile=")) tracefp=fopen(strtok(&arg[10],newline), "a");
   else if (IS_ARG_MATCH(arg,"movie")) { export_mode=2; export_movie=1; }
   else if (IS_ARG_MATCH(arg,"tmax=")) tmax=atof(&arg[5]);
   else if (IS_ARG_MATCH(arg,"emax=")) emax=atoi(&arg[5]);
   else if (IS_ARG_MATCH(arg,"smax=")) smax=atoi(&arg[5]);
   else if (IS_ARG_MATCH(arg,"mmax=")) mmax=atoi(&arg[6]);
   else if (IS_ARG_MATCH(arg,"smin=")) smin=atoi(&arg[5]);
   else if (IS_ARG_MATCH(arg,"fsmax=")) fsmax=atoi(&arg[6]);
   else if (IS_ARG_MATCH(arg,"untiltiles=")) {
      int i=0;
      char *pos;
      pos = arg;
      while (pos != NULL) {
	 pos = strchr (pos+1,',');
	 present_list_len++;
      }
      present_list = (int *) malloc(present_list_len*sizeof(int));
      pos = &arg[11];
      while ((pos-1) != NULL) {
	 present_list[i++] = atoi(pos);
	 pos = strchr (pos,',') + 1;
      }
      untiltiles = 1;
   }
   else if (IS_ARG_MATCH(arg,"untiltilescount=")) {
      int i=0;
      char *pos;
      untiltilescount=1;
      pos = arg;
      while (pos != NULL) {
	 pos = strchr (pos+1,',');
	 present_list_len++;
      }
      present_list = (int *) malloc(present_list_len*sizeof(int));
      pos = &arg[16];
      while ((pos-1) != NULL) {
	 present_list[i++] = atoi(pos);
	 pos = strchr (pos,',') + 1;
      }
      untiltilescount = 1;
   }
   else if (IS_ARG_MATCH(arg,"clean_cycles=")) clean_cycles=atoi(&arg[13]);
   else if (IS_ARG_MATCH(arg,"clean_X=")) clean_X=atof(&arg[8]);
   else if (IS_ARG_MATCH(arg,"fill_cycles=")) fill_cycles=atoi(&arg[12]);
   else if (IS_ARG_MATCH(arg,"fill_X=")) fill_X=atof(&arg[7]);
   else if (IS_ARG_MATCH(arg,"error_radius=")) error_radius=atof(&arg[13]);
   else if (IS_ARG_MATCH(arg,"repair_unique_T=")) { repair_unique_T=atof(&arg[15]); repair_unique=1; }
   else if (IS_ARG_MATCH(arg,"datafile=")) datafp=fopen(strtok(&arg[9],newline), "a");
   else if (IS_ARG_MATCH(arg,"largeflakedatafile=")) {
      char *c;
      c =strchr(arg,',') + 1;
      tthresh = atoi(&arg[19]);
      largeflakefp=fopen(c, "a");
   }
   else if (IS_ARG_MATCH(arg,"untiltilescountfile=")) {
      untiltilescountfp=fopen(&arg[20], "a");
   }
   else if (IS_ARG_MATCH(arg,"arrayfile=")) arrayfp=fopen(strtok(&arg[10],newline), "w");
   else if (IS_ARG_MATCH(arg,"exportfile=")) export_fp=fopen(strtok(&arg[11],newline), "w");
   else if (strncmp(arg,"testing",7) == 0) {
      testing = 1;
   }
   else if (IS_ARG_MATCH(arg,"import_offset=")) {
      char *p=(&arg[14]);
      import_offset_i=atoi(p);
      if ((p=strchr(p,','))!=NULL) {
	 import_offset_j=atoi(p+1);
      }
   }
   else if (IS_ARG_MATCH(arg,"importfile"))
   {
      char *p=(&arg[11]);
      FILE *import_fp;
      import = 1;
      fprm = (struct flake_param *)malloc(sizeof(struct flake_param));
      fprm->seed_i = fprm->seed_j = 130;
      fprm->seed_n = 1;
      fprm->N = 1;
      fprm->Gfc=0;
      fprm->import_from=NULL;
      fprm->next_param = fparam;
      fparam = fprm;
      /*
       * Need to figure out seed_i and seed_j after loading in the cells.
       * Randomly choose spots from the smallest flake until we find a tile
       * that is present in each flake.
       */
      fparam->seed_i = fparam->seed_j = size / 2;
      fparam->seed_n = 1;
      /* It's easier to use characters here than fight cvs's varying
	 habits with non-alphabetic characters */
      newline[0] = 10;
      newline[1] = 0;
      if (IS_ARG_MATCH(arg,"importfile=")) {
	 char imp_fn[256], *arg_fn;
	 import_fp = fopen(arg_fn=strtok(p,newline), "r");
	 if (import_fp == NULL) { sprintf(&imp_fn[0],"%s",arg_fn); import_fp = fopen(&imp_fn[0],"r"); }
	 if (import_fp == NULL) { sprintf(&imp_fn[0],"%s.seed",arg_fn); import_fp = fopen(&imp_fn[0],"r"); }
	 if (import_fp == NULL) { sprintf(&imp_fn[0],"tilesets/%s",arg_fn); import_fp = fopen(&imp_fn[0],"r"); }
	 if (import_fp == NULL) { sprintf(&imp_fn[0],"tilesets/%s.seed",arg_fn); import_fp = fopen(&imp_fn[0],"r"); }
      } else {
	 import_fp = fopen("xgrow_export_output", "r");
      }
      /* If the file does not exist. */
      if (import_fp == NULL) 
      {
	 fprintf(stderr, "Error: Import file not found.\n");
	 return -1;
      }
      fparam->import_from = import_fp;
      if ((p = strtok(NULL, newline)) != NULL)
	 fparam->Gfc=atof(p);
      else
	 fparam->Gfc=0;
      fparam->N = count_flakes(import_fp);
   }
   else if (IS_ARG_MATCH(arg,"min_strength=")) {min_strength=atof(&arg[13]);}
   else {
      fprintf(stderr,"Could not parse argument '%s'\n",arg); 
      return -1;
   }
   return 0;
}

#define rsc read_skip_comment(tilefp)
void read_skip_comment(FILE *fp)
   // anything after a "%" gets ignored
{ int c;
   fscanf(fp," "); c=fgetc(fp); 
   while (c=='%') {
      fgets(&stringbuffer[0],256,fp); 
      // fprintf(stderr,"%%%s",stringbuffer);
      fscanf(fp," "); c=fgetc(fp); 
   }
   ungetc(c,fp); 
}

void read_tilefile(FILE *tilefp) 
{ 
   float strength_float, glue_float, stoic_float; int i,j,k;
   int temp_char; char s[255], **btnames;
   int n,m,r;
   int return_code;

   rewind(tilefp); 
   // needs to be read twice, 
   //  once to get sim specs, field size, etc
   //  once to get colors after X window is open
   rsc;
   fscanf(tilefp,"tile edges matches {{N E S W}*}\n"); rsc;

   if (1!=fscanf(tilefp,"num tile types=%d\n",&N)) 
   { fprintf(stderr,"Reading tile file: expected num tile types.\n"); exit(-1); } 
   rsc;
   if (1!=fscanf(tilefp,"num binding types=%d\n", &num_bindings))
   { fprintf(stderr,"Reading tile file: expected num binding types.\n"); exit(-1); } 
   rsc;

   if (N > MAXTILETYPES) 
   { fprintf(stderr,"Reading tile file: too many tile types.\n"); exit(-1); }
   rsc;
   if (num_bindings > MAXTILETYPES) 
   { fprintf(stderr,"Reading tile file: too many binding types.\n"); exit(-1); } 
   rsc;

   btnames=(char **)malloc((num_bindings+1)*sizeof(char *));  
   for (k=0;k<=num_bindings;k++) btnames[k]="null";   // until overwritten by tile file specification;

   temp_char = getc (tilefp); ungetc(temp_char, tilefp);
   if (temp_char == 'b') {
      r=0; fscanf(tilefp,"binding type names=%n",&r); rsc;
      if (r!=19) { fprintf(stderr,"Reading tile file: expected binding type names.\n"); exit(-1); }
      r=0; fscanf(tilefp,"{%n",&r); rsc;
      if (r!=1) { fprintf(stderr,"Reading tile file: expected binding type names {.\n"); exit(-1); }

      for (i=1;i<=num_bindings;i++) {
	 if (1!=fscanf(tilefp,"%100s",&s[0]))
	 { fprintf(stderr,"Reading tile file: expected binding type %d name.\n",i); exit(-1); }
	 if (s[strlen(s)-1]=='}') { ungetc('}',tilefp); s[strlen(s)-1]=0; }
	 if (strlen(s)==0)
	 { fprintf(stderr,"Reading tile file: expected binding type %d name.\n",i); exit(-1); }
	 btnames[i]=malloc(strlen(s)); strcpy(btnames[i],s); rsc;
	 if (index("0123456789",btnames[i][0])!=NULL)
	 { fprintf(stderr,"Reading tile file: binding type %d name ('%s') cannot begin with 0-9.\n",i,s); exit(-1); }
	 // printf("read bond type %d name = '%s'\n",i,btnames[i]);
      } rsc;
      r=0; fscanf(tilefp,"}%n",&r); rsc;
      if (r!=1) { fprintf(stderr,"Reading tile file: expected binding type names }.\n"); exit(-1); }
   }

   r=0; fscanf(tilefp,"tile edges=%n",&r); rsc; //printf("r=%d\n",r);
   if (r!=11) { fprintf(stderr,"Reading tile file: expected `tile edges=' declaration.\n"); exit(-1); }
   tileb_length = N+1;
   tileb = (int**) calloc(sizeof(int*),tileb_length);
   stoic = (double*) calloc(sizeof(double),tileb_length);
   r=0; fscanf(tilefp,"{%n",&r); rsc;
   if (r!=1) { fprintf(stderr,"Reading tile file: expected tile set start {.\n"); exit(-1); }
   tileb[0] = (int*) calloc(sizeof(int),4);
   for (j=0;j<4;j++) {
      tileb[0][j] = 0;
   }
   stoic[0]=0;
   for (i=1;i<tileb_length;i++) {
      tileb[i] = (int*) calloc(sizeof(int),4);
      r=0; fscanf(tilefp,"{%n",&r); rsc; 
      if (r!=1) { fprintf(stderr,"Reading tile file: expected tile %d start def {. \n",i); exit(-1); }
      /* read in the four binding types {N E S W} */
      for (j=0;j<4;j++) { rsc;
	 temp_char = getc (tilefp); ungetc(temp_char, tilefp);
	 if (index("0123456789",temp_char)!=NULL) {
	    if (1!=fscanf(tilefp,"%d",&tileb[i][j]))
	    { fprintf(stderr,"Reading tile file: expected tile %d bond %d value.\n",i,j); exit(-1); } 
	 } else {
	    if (1!=fscanf(tilefp,"%[^} \t]",&s[0]))
	    { fprintf(stderr,"Reading tile file: expected tile %d bond %d's name.\n",i,j); exit(-1); } 
	    for (k=0;k<=num_bindings;k++) if (strcmp(s,btnames[k])==0) break;
	    if (k<=num_bindings)
	    { tileb[i][j]=k; }
	    else
	    { fprintf(stderr,"Reading tile file: expected tile %d bond %d's name '%s' unknown.\n",i,j,s); exit(-1); } 
	 }
      } rsc;
      r=0; fscanf(tilefp,"}%n",&r); rsc; 
      if (r!=1) { fprintf(stderr,"Reading tile file: expected tile %d end def }. \n",i); exit(-1); }
      if (fscanf(tilefp,"[%g]",&stoic_float)) 
	 stoic[i]=stoic_float; else stoic[i]=1.0; rsc;
      if (fscanf(tilefp,"(%200[^)])",&stringbuffer[0])) {
	 tile_colors[i]=(char *)malloc(strlen(stringbuffer)+2);
	 strcpy(tile_colors[i],stringbuffer);
      }
      fscanf(tilefp,"\n"); rsc;
   }
   r=0; fscanf(tilefp,"}%n",&r); rsc; 
   if (r!=1) { fprintf(stderr,"Reading tile file: expected tile set end }.\n"); exit(-1); }
   // printf("Tile set loaded (%d tiles)\n",N);
   // for (i=1;i<tileb_length;i++) {
   //   for (j=0;j<4;j++) printf("%d ",tileb[i][j]); printf("\n");
   // }


   glue = (double **) calloc(sizeof (double *),num_bindings + 1);
   for (i=0;i<=num_bindings;i++) { int j;
      glue[i] = (double *) calloc(sizeof (double),num_bindings + 1);
      for (j=0;j<=num_bindings;j++) { glue[i][j]=0; }  // necessary??  -- EW
   }
   strength = (double*) calloc(sizeof(double),num_bindings+1); 
   dt_left = (int *) calloc(sizeof(int),N+1);
   dt_right = (int *) calloc(sizeof(int),N+1);
   dt_up = (int *) calloc(sizeof(int),N+1);
   dt_down = (int *) calloc(sizeof(int),N+1);

   temp_char = getc (tilefp);  ungetc(temp_char, tilefp);
   if (temp_char == 'b') {
      r=0; fscanf(tilefp,"binding strengths=%n",&r); rsc;
      if (r!=18) { fprintf(stderr,"Reading tile file: expected binding strength defs.\n"); exit(-1); }
      r=0; fscanf(tilefp,"{%n",&r); rsc;
      if (r!=1) { fprintf(stderr,"Reading tile file: expected binding strength defs {.\n"); exit(-1); }
      strength[0]=0; /* bond type 0 ("null") always has strength 0 */
      for (i=1;i<=num_bindings;i++) {
	 if (1!=fscanf(tilefp,"%g",&strength_float))
	 { fprintf(stderr,"Reading tile file: expected binding strength %d value.\n",i); exit(-1); }
	 strength[i]=(double)strength_float;
	 // printf("strength for se #%d = %g\n",i,strength[i]);
      } rsc;
      r=0; fscanf(tilefp,"}%n",&r); rsc;
      if (r!=1) { fprintf(stderr,"Reading tile file: expected binding strength defs }.\n"); exit(-1); }
   }
   while ((temp_char=fgetc(tilefp)) == 'g') {
      if (3!=fscanf(tilefp,"(%d,%d)=%g\n",&n,&m,&glue_float))
      { fprintf(stderr,"Reading tile file: expected glue def.\n"); exit(-1); } 
      rsc;
      // printf ("Glue float is g(%d,%d)=%g\n",n,m,glue_float);
      glue[n][m] = (double) glue_float;
      glue[m][n] = (double) glue_float;
   } ungetc(temp_char, tilefp); rsc;
   // printf("Bond strengths loaded (%d bond types)\n",num_bindings);
   // for (i=1;i<=num_bindings;i++) printf("%f ",strength[i]); printf("\n");

   rsc;
   while(fgets(&stringbuffer[0],256,tilefp)!=NULL) {
      return_code = parse_arg_line(&stringbuffer[0]); 
      if (return_code) {
	 exit(-1);
      }
      rsc;
   }

   fclose(tilefp);

   // better not free this: some strings are constants, some are alloc'ed.
   // for (i=1; i<=num_bindings;i++) free(btnames[i]); free(btnames);
}

void getargs(int argc, char **argv)
{
   int i; struct flake_param *fprm;
   struct timeval tv; 
   gettimeofday(&tv, NULL); srand48(tv.tv_usec); srandom(tv.tv_usec);
   /* NOTE: Disabled to allow compilation on 64-bit systems; doesn't seem to cause problems. (cge, 091028)
      if (sizeof(long) != 4) {
      printf("Error: sizeof long (%d) should be 4\n", (int)sizeof(long int));
      exit(-1);
      }
      */
   if (argc==2 && strcmp(argv[1],"--")==0) {
      printf("usage: xgrow tilefile [option=#]... \n");
      printf(" tilefile is an input file that specifies tiles\n");
      printf(" options:\n");
      printf("  block=  display block size, 1...10\n");
      printf("  size=   field side length (power-of-two) [default 256]\n");
      printf("  rand=   random number seed\n");
      printf("  T=      threshold T (relative to Gse) for irreversible Tile Assembly Model\n");
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
      printf("  anneal=g,t            anneal Gse from g to given Gse with time constant t\n"); 
      printf("  linan=h,s,C0,Cfin,dt  do a linear anneal, where delta G at each temperature is\n");
      printf("                        governed by (Delta) h and (Delta) s in kcals per mol.  \n");
      printf("                        Go from C0 > Cfin changing temps 1 degree every dt secs \n");
      printf("                        (in incremenets of 0.1C).  Ignores Gse value.\n");
      printf("  seed=i,j,n            seed tile type n at position i,j\n");
      printf("  tinybox=k             use dynamic flakes in a box containing k of each molecule on average.\n");
      printf("  addflakes=i,j,n:N@Gfc simulate N separate flakes\n");
      printf("  stripe=o[:p,w]*       width w stripe with p errors, offset o\n");
      printf("  wander                wandering `seed' designation\n");
      printf("  fission               can tile be removed if two flakes result?\n");
      printf("  no_fission             the answer is no [default]\n");
      printf("  chunk_fission         allow pairs & 2x2 blocks to dissociate as one (implies fission)\n");
      printf("  blast_rate_alpha      square kxk holes are blasted with this per-tile rate for 1x1 [default=0]\n"
	    "                         (rate relative to tile addition, i.e. scaled by total concentration)\n");
      printf("  blast_rate_beta        rate scales as 1/k^beta [default=4]\n");
      printf("  blast_rate_gamma       rate also scales as exp(-gamma*(k-1)) [default=0]\n");
      printf("  zero_bonds            can tiles be added if they bond with 0 strength?\n");
      printf("  no_zero_bonds          the answer is no [default]\n");
      printf("  min_strength=         set minimum bond energy below which attachments are considered incorrect\n");
      printf("  periodic              periodic boundary conditions\n");
      printf("  -linear               simulate linear A B tiles, write errs > stdout \n");
      printf("  -nw                   no X window (only if ?max set)\n");
      printf("  update_rate=          update display every so-many events\n");
      printf("  tracefile=            append datafile info (see below) EVERY so-many events\n");
      printf("  movie                 export MATLAB-format flake array information EVERY so-many events\n");
      printf("  tmax=                 quit after time t has passed\n");
      printf("  emax=                 quit after e events have occurred\n");
      printf("  smax=                 quit when the fragment or total size of fragments is size s\n");
      printf("  mmax=                quit when there are at least mm mismatches by xgrow's interpretation\n");
      printf("  smin=                 quit when the fragment or total size of fragments goes to or below size s\n");
      printf("  fsmax=                 quit when a single fragment reaches size s\n");
      printf("  untiltiles=           quit when all (numbered) tiles in the comma-delineated list are in the assembly\n");
      printf("  clean_cycles=         at end, remove how many layers of weakly attached tiles?\n"
	    "                        [default=0]\n");
      printf("  clean_X=              for cleaning, minimal ratio of off-rate to on-rate [default=1.0]\n");
      printf("  fill_cycles=          at end, add how many layers of strongly attached tiles?\n"
	    "                        [default=0]\n");
      printf("  fill_X=               for filling, minimal ratio of off-rate to on-rate [default=1.0]\n");
      printf("  error_radius=         when writing to file, #mismatches counts only those for which \n"
	    "                        all surounding tiles are present (after clean/fill) [default=0]\n");
      printf("  repair_unique_T=      alternative clean/fill, called Rx: remove mismatches, fill in interior sites \n"
	    "                        if there is a unique strength-T tile, then fill in by strongest tile\n");
      printf("  datafile=             append Gmc, Gse, ratek, time, size, #mismatched se, events, perimeter, dG, dG_bonds for each flake\n");
      printf("  arrayfile=            output MATLAB-format flake array information on exit (after cleaning)\n");
      printf("  exportfile=           on-request output of MATLAB-format flake array information\n");
      printf("                        [defaults to 'xgrow_export_output']\n");
      printf("  importfile=FILENAME   import all flakes from FILENAME.\n");
      printf("  importfile            import all flakes from xgrow_export_output.\n");
      printf("  pause                 start in paused state; wait for user to request simulation to start.\n");
      printf("  testing               run automated tests instead of a simulation.\n");
      exit (0);
   }
   if (argc == 1) {
      printf("* First argument must be a tile file!\nTry 'xgrow --' for help.\n");
      exit(1);
   }

   tmax=0; emax=0; smax=0; fsmax=0; smin=-1; mmax=0;
   wander=0; periodic=0; linear=0; fission_allowed=0; zero_bonds_allowed=0;
   Gfc=0; datafp=NULL; arrayfp=NULL;  largeflakefp=NULL; untiltilescountfp=NULL;
   Gmc=17; Gse=8.6; ratek = 1000000.0;  T=0;
   Gmch=30; Gseh=0; Ghyd=30; Gas=30; Gam=15; Gae=30; Gah=30; Gao=10;
   seed_i=250; seed_j=250; seed_n=1; hydro=0;


   // reset tile colors that haven't been assigned
   for (i=15; i<MAXTILETYPES; i++) tile_colors[i]=NULL;

   if      ( (sprintf(&tileset_name[0],"%s",argv[1]),tilefp = fopen(&tileset_name[0],"r"))!=NULL ) read_tilefile(tilefp); 
   else if ( (sprintf(&tileset_name[0],"%s.tiles",argv[1]),tilefp = fopen(&tileset_name[0],"r"))!=NULL ) read_tilefile(tilefp); 
   else if ( (sprintf(&tileset_name[0],"tilesets/%s",argv[1]),tilefp = fopen(&tileset_name[0],"r"))!=NULL ) read_tilefile(tilefp); 
   else if ( (sprintf(&tileset_name[0],"tilesets/%s.tiles",argv[1]),tilefp = fopen(&tileset_name[0],"r"))!=NULL ) read_tilefile(tilefp); 
   else {
      printf("* First argument must be a tile file!\nTry 'xgrow --' for help.\n");
      exit(0);   
   }

   for (i=2; i<argc; i++) {
      parse_arg_line(argv[i]);
   }
   if (tmax==0 && emax==0 && smax==0 && mmax==0 && fsmax==0 && smin==-1) {
     printf("No max setting: forcing UI mode.\n");
     XXX=1;
   }
   if (hydro && fission_allowed==2) {
      printf("* Current implementation does not allow chunk_fission and hydrolysis simultaneously.\n"); exit(0);
   }
   if (double_tiles || vdouble_tiles) {
      //if (fission_allowed != 0) {
      //printf("Double tiles cannot be used with fission or chunk_fission currently.\n");
      //exit(0);
      //}
      if (hydro) {
	 printf("Double tiles cannot be used with hydrolysis currently.\n");
	 exit(0);
      }
      if (fill_cycles || repair_unique) {
	 printf("Double tiles cannot be used with repair or fill cycles currently.\n");
	 exit(0);
      }
      if (blast_rate > 0) {
	 printf("Blasting cannot be used with double tiles currently.\n");
	 exit(0);
      }
   }

   for (size_P=5; (1<<size_P)<size; size_P++);
   size=(1<<size_P); 
   //if (XXX) {
   //   if (size*block > 800) block=800/size;
   //   if (block==0) { size=512; block=1; size_P=9; }
   //}

   if (blast_rate_alpha>0) { int kb;
      printf("blast_rate: alpha = %f, beta = %f, gamma = %f\n",blast_rate_alpha,blast_rate_beta,blast_rate_gamma);
      for (kb=1; kb<size; kb++)  
	 blast_rate += blast_rate_alpha * exp(-blast_rate_gamma*(kb-1)) / pow(kb*1.0,blast_rate_beta) ;
      printf("total blast_rate per site = %f\n",blast_rate);
   }


   NROWS=(size+NBDY*2);
   NCOLS=(size+NBDY*2);
   VOLUME=(NROWS*NCOLS);
   WINDOWWIDTH=(MAX(block*NCOLS,256)+PLAYLEFT+BOUNDWIDTH+30);
   WINDOWHEIGHT=(PLAYTOP+MAX(block*NROWS,256)+100);

   T=T*Gse;

   /* for consistency's sake, insert seed_i, seed_j, seed_n, Gfc into fparam */
   /* if no addflake commands were issued */
   if (fparam==NULL) {
      fprm = (struct flake_param *)malloc(sizeof(struct flake_param)); 
      fprm->seed_i=seed_i; fprm->seed_j=seed_j; fprm->seed_n=seed_n; 
      fprm->import_from=NULL;
      fprm->N= (tinybox == 0); 
      fprm->Gfc=Gfc;
      fprm->next_param=fparam;
      fparam=fprm;
   }

   // make sure all seed tiles are on the board
   for (fprm=fparam; fprm!=NULL; fprm=fprm->next_param) {
      while (fprm->seed_i>=size) fprm->seed_i/=2;
      while (fprm->seed_j>=size) fprm->seed_j/=2;
   }
   if (fparam != NULL) {
      seed_i=fparam->seed_i; seed_j=fparam->seed_j;
   }

   if (!XXX) {
     fprintf(stderr, " Starting simulation (1st seed=%d,%d,%d) on %d x %d board.\n",
	    seed_i, seed_j, seed_n, size, size);
   }

}

void write_datalines(FILE *out, char *text)
{ flake *fpp; int perimeter; double dG_bonds;  

   for (fpp=tp->flake_list; fpp!=NULL; fpp=fpp->next_flake) {
      if (strcmp(text,"")==0) fpp=fp;
      perimeter=calc_perimeter(fpp);
      dG_bonds = calc_dG_bonds(fpp);
      if (tp->hydro) fprintf(out, " %f %f %f %f %f %f %f %f %f ",
	    Gseh, Gmch, Ghyd, Gas, Gam, Gae, Gah, Gao, Gfc);
      fprintf(out, " %f %f %f %f %d %d %lld %d %f %f%s",
	    Gmc,tp->Gse,ratek,tp->t,fpp->tiles,fpp->mismatches,tp->events,
	    perimeter, fpp->G, dG_bonds,text);
      if (strcmp(text,"")==0) break;
   }
   fflush(out);
}

void write_largeflakedata(FILE *filep) {
   flake *fp;
   int large_flakes = 0;
   for (fp = tp->flake_list; fp != NULL; fp=fp->next_flake) {
      if (fp->tiles > tthresh) {
	 large_flakes++;
      }
   }
   fprintf(filep,"%d\n",large_flakes);
}

void write_untiltilescountdata(FILE *filep) {
   fprintf(filep,"%d\n",tp->untiltilescount);
}

void write_flake(FILE *filep, char *mode, flake *fp)
{ int n, row, col, tl; 

   if (filep!=NULL) {
      if (strcmp(mode,"flake")==0) n=export_flake_n++;
      else if (strcmp(mode,"movie")==0) n=export_movie_n++;
      else n=1;
      fprintf(filep,"\n%s{%d}={ ...\n",mode,n);
      fprintf(filep,"[ "); write_datalines(filep,"");
      fprintf(filep," ],...\n  [");
      for (tl=1; tl<=tp->N; tl++) 
	 fprintf(filep," %8.5f", (tp->conc[tl]>0)?-log(tp->conc[tl]):0 );
      fprintf(filep," ],...\n  [");
      for (row=0; row<size; row++) {
	 for (col=0; col<size; col++)
	    fprintf(filep, " %d", fp->Cell(row,col));
	 fprintf(filep, "; ...\n");
      }
      fprintf(filep," ] };\n\n");
      fflush(filep);
   }
}  

void export_flake(char *mode, flake *fp)
{
   if (export_fp==NULL) export_fp=fopen("xgrow_export_output","a+");
   printf("Writing flake #%d...\n",export_flake_n);
   write_flake(export_fp, mode, fp);
   fflush(export_fp);
}


/*    READ in by Shaun Lee
 * Count the number of flakes in an exported file.
 * Doesn't actually count the number of flakes, it just watches
 * the flake number and returns the last one it finds.
 */
int count_flakes(FILE *flake_file)
{
   int n, i, flake_size, row, garbage, end;
   /*
    * This assumes the parameters (which we want to ignore)
    * are less than 20 charaters long (Overkill.  14 should be enough)
    */
   char line[20]; int lnum=0;
   flake_size = 0;

   /* Run through once to make sure the flakes are the right size. */

   fscanf(flake_file, "\nflake{%d}={ ...\n[", &n); lnum+=2;   // the line number counting is very approximate.

   /* For debugging */
   fprintf(stderr, "Reading flake number: %d\n", n);   

   /* Run through the parameters, waiting for ],... to appear twice. */
   for(i = 0; i < 2; i++)
   {
      while (!feof(flake_file))
      {
	 fscanf(flake_file, "%s", line); 
	 if (strcmp(line, "],...")==0)
	    break;
      }
   }
   if (feof(flake_file)) { 
      fprintf(stderr,"Reached EOF on flake number %d without reading params.\n",n);
      exit(-1);      
   }

   /* Read in the '[' */
   fscanf(flake_file, "%s", line); lnum+=3;  // for params  (if written by xgrow)

   while (1)
   {
      fscanf(flake_file, "%s", line);  
      // printf("Just read in loop near line %d: %s\n", lnum, line); //****
      if (strcmp(line, "...")==0) /* then we've reached the end of the line */
      {
	 break; lnum++;
      }
      else
	 flake_size++;
   }

   if (flake_size > import_flake_size)
      import_flake_size = flake_size;

   /*
    * Now finish processing the first flake.
    */

   for(row = 2; row <= flake_size; row++)
   {
      //        fprintf(stderr, "Now on row: %d\n", row);
      for(i = 0; i <= flake_size; i++)
      {
	 fscanf(flake_file, "%d", &garbage);
	 /*	  fprintf(stderr, "%d:%d ", i, garbage); */
      }
      fscanf(flake_file, "%s", line); 
      /* Make sure that we've actually read the end of line. */
      if(strcmp(";", line) != 0)
      {
	 fprintf(stderr, "Error parsing input file near line %d. Expected ';' but read '%s'\n", lnum, line);
	 exit(-1);
      }
      fscanf(flake_file, "%s", line); lnum++;
      if(strcmp("...", line) != 0)
      {
	 fprintf(stderr, "Error parsing input file near line %d. Expected '...' but read '%s'\n", lnum,line);
	 exit(-1);
      }
   }

   row--; /* stupid off by one things. */

   if(row != flake_size)
   {
      fprintf(stderr, "Error near line %d: Flake dimensions do not match.\n%d rows and "
	    "%d columns\n", lnum, row, flake_size);
      exit(-1);
   }

   fscanf(flake_file, "%s", line); 
   assert(strcmp(line, "]")==0);
   fscanf(flake_file, "%s", line); lnum++;
   assert(strcmp(line, "};")==0);

   /* Watch for end of file as well as find the next flake number. */
   end = fscanf(flake_file, "\nflake{%d}={ ...\n[", &n); lnum+=2;

   /* Now go through the rest of the file. */
   while (end==1)
   {
      /* For debugging */
      printf("Reading flake number: %d\n", n); //****

      /* Run through the parameters, waiting for ],... to appear twice. */
      for(i = 0; i < 2; i++)
      {
	 while (!feof(flake_file))
	 {
	    fscanf(flake_file, "%s", line); 
	    if (strcmp(line, "],...")==0)
	       break;
	 }
      }
      if (feof(flake_file)) { 
	 fprintf(stderr,"Reached EOF on flake number %d without reading params.\n",n);
	 exit(-1);      
      }

      printf(" skipped param lines\n"); 
      fscanf(flake_file, "%s", line);  lnum+=3;
      assert(strcmp(line, "[")==0);
      for(row = 1; row <= flake_size; row++)
      {
	 for(i = 0; i < flake_size; i++)
	 {
	    fscanf(flake_file, "%s", line); 
	 }
	 fscanf(flake_file, "%s", line); lnum++;
	 assert(strcmp("...", line)==0);
      } 

      fscanf(flake_file, "%s", line); 
      assert(strcmp("]", line)==0);
      fscanf(flake_file, "%s", line); lnum++;
      assert(strcmp("};", line)==0);

      /* Watch for end of file as well as find the next flake number. */
      end = fscanf(flake_file, "\nflake{%d}={ ...\n[", &n);  lnum+=2;
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
void import_flake(flake *current_flake, FILE *flake_file, int flake_number)
{
   fpos_t flake_start;
   int tile_type, read_flake_number, translate_i, translate_j, i, j, end, flake_size, seed_set;
   char line[20]; int lnum=0;
   seed_set = flake_size = 0;

   /* Just to be sure. */
   rewind(flake_file);

   /* Find the flake we want. */
   fscanf(flake_file, "\nflake{%d}={ ...\n[", &read_flake_number); lnum++;  // the line number counting is very approximate.

   while (read_flake_number!=flake_number)
   {
      /* Skip through the parameters, waiting for ],... to appear twice. */
      for(i = 0; i < 2; i++)
      {
	 while (1)
	 {
	    /* Assumes the parameters are terminated by '],...' */
	    if (fscanf(flake_file, "%s", line)==0) 
	    {
	       fprintf(stderr,"Error near line %d: Expected ],... at end of parameters.\n",lnum+1);
	       exit(-1);
	    }
	    lnum++;
	    if (strcmp(line, "],...")==0)
	       break;
	 }
      }

      /* Skip through the data for this flake */
      while(1)
      {
	 if (fscanf(flake_file, "%s", line)==0)
	 {
	    fprintf(stderr,"Error near line %d: Expected }; at end of data.\n",lnum+1);
	    exit(-1);
	 }
	 lnum++;
	 if (strcmp(line, "};")==0)
	    break;
      }

      end = fscanf(flake_file, "\n\nflake{%d}={", &read_flake_number);

      if (!end)
      {
	 fprintf(stderr, "Flake %d not found in file (near line %d).\n", flake_number,lnum);
	 return;
      }
      lnum+=2;
   }

   /* Skip through the parameters, waiting for ],... to appear twice. */
   for(i = 0; i < 2; i++)
   {
      while (1)
      {
	 fscanf(flake_file, "%s", line); lnum++;
	 if (strcmp(line, "],...")==0)
	    break;
      }
   }
   /* Read in the starting '[' */
   fscanf(flake_file, "%s", line);  lnum++;
   assert(strcmp("[", line)==0);

   /*
    * Save the current position in the file stream so that we can
    * import data values after calculating the size of the import
    * flake.
    */
   fgetpos(flake_file, &flake_start);

   /* Find the width of the flake. */
   while (1)
   {
      fscanf(flake_file, "%s", line);
      if (strcmp(line, "...")==0)
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
   translate_i = (size - flake_size) / 2 + import_offset_i;
   translate_j = (size - flake_size) / 2 + import_offset_j;

   /* we already plopped down a seed tile.  erase it if the loaded assemble won't. */
   change_cell(current_flake, current_flake->seed_i, current_flake->seed_j, 0);

   for (i = 0; i < flake_size; i++)
   {
      for (j = 0; j < (flake_size - 1); j++)
      {
	 fscanf(flake_file, "%d", &tile_type);
	 change_cell(current_flake, translate_i + i, translate_j + j, tile_type);
      }
      fscanf(flake_file, "%d;", &tile_type);
      change_cell(current_flake, translate_i + i, translate_j + j, tile_type);

      fscanf(flake_file, "%s", line); lnum++;
      assert(strcmp(line, "...")==0);
   }

   /* Now choose a random active cell and set it as the flake's seed. */
   srand(time(0));
   i = flake_size * (((double)rand()) / ((double)RAND_MAX)) + translate_i;
   j = flake_size * (((double)rand()) / ((double)RAND_MAX)) + translate_j;
   while ((current_flake->Cell(i,j)) == 0 || tp->dt_left[current_flake->Cell(i,j)]) //FIXME: vdouble
   {
      i = flake_size * (((double)rand()) / ((double)RAND_MAX)) + translate_i;
      j = flake_size * (((double)rand()) / ((double)RAND_MAX)) + translate_j;
   }
   current_flake->seed_i = i;
   current_flake->seed_j = j;
   current_flake->seed_n = current_flake->Cell(i,j);

   /*    fprintf(stderr, "Just set seed for flake %d to (%d,%d), tile type %d\n", flake_number, i, j, current_flake->Cell(i,j)); */
   fp->events = 0;

   recalc_G(fp);
}


void closeargs()
{ 
   int i;  flake *fpp;

   // cleans all flakes  (removes "temporary" tiles on growth edge)
   for (fpp=tp->flake_list; fpp!=NULL; fpp=fpp->next_flake) 
      clean_flake(fpp,clean_X,clean_cycles); 

   // fills in all flakes  (adds tiles in holes or on growth edge that "ought" to be there)
   for (fpp=tp->flake_list; fpp!=NULL; fpp=fpp->next_flake) 
      fill_flake(fpp,fill_X,fill_cycles); 

   // alternative fix-up routine; redundant with clean/fill, but supposedly more reliable
   if (repair_unique) 
      for (fpp=tp->flake_list; fpp!=NULL; fpp=fpp->next_flake) 
	 repair_flake(fpp,repair_unique_T,Gse); 

   // recalc's all mismatch #'s if error_radius is given -- keep only those not near an empty cell.
   if (error_radius>.5) {
      for (fpp=tp->flake_list; fpp!=NULL; fpp=fpp->next_flake) 
	 error_radius_flake(fpp,error_radius); 
   }

   /* output information for *all* flakes */
   if (datafp!=NULL) {
      write_datalines(datafp,"\n"); fclose(datafp);
   } 

   if (largeflakefp!=NULL) {
      write_largeflakedata(largeflakefp); fclose(largeflakefp);
   }

   if (untiltilescountfp!=NULL) {
      write_untiltilescountdata(untiltilescountfp);
      fclose(untiltilescountfp);
   }
   if (arrayfp!=NULL) {
      export_flake_n=1;
      for (fpp=tp->flake_list; fpp!=NULL; fpp=fpp->next_flake) {
	 write_flake(arrayfp, "flake", fpp);
      }
      fclose(arrayfp);
   }
   if (export_fp!=NULL) fclose(export_fp);

   // free memory
   for (i=0;i<=tp->N;i++) free(tileb[i]); free(tileb);
   for (i=0;i<=tp->num_bindings;i++) free(glue[i]); free(glue);
   free(strength); free(stoic);

   while (fparam!=NULL) { fprm=fparam->next_param; free(fparam); fparam=fprm; }

   free_tube(tp);  // this frees all the flakes too
}

// (i,j) should not be empty. returns 0 if OK, 1 if mismatches with E or S or N or W, or unbound s.e.
#define errortile(i,j) ( (                                                      \
	 (tileb[fp->Cell(i,j)][1] != tileb[fp->Cell(i,(j)+1)][3] && \
	  tileb[fp->Cell(i,j)][1]*tileb[fp->Cell(i,(j)+1)][3] > 0 && \
	  (fp->tube->glue)[tileb[fp->Cell(i,j)][1]][tileb[fp->Cell(i,(j)+1)][3]] < min_strength) || \
	 (tileb[fp->Cell(i,j)][3] != tileb[fp->Cell(i,(j)-1)][1] &&           \
	  tileb[fp->Cell(i,j)][3]*tileb[fp->Cell(i,(j)-1)][1] > 0 &&   \
	  (fp->tube->glue)[tileb[fp->Cell(i,j)][3]][tileb[fp->Cell(i,(j)-1)][1]] < min_strength) || \
	 (tileb[fp->Cell(i,j)][2] != tileb[fp->Cell((i)+1,j)][0] &&           \
	  tileb[fp->Cell(i,j)][2]*tileb[fp->Cell((i)+1,(j))][0] > 0 &&   /* This one works!? */ \
	  (fp->tube->glue)[tileb[fp->Cell(i,j)][2]][tileb[fp->Cell((i)+1,j)][0]] < min_strength) || \
	 (tileb[fp->Cell(i,j)][0] != tileb[fp->Cell((i)-1,j)][2] &&           \
	  tileb[fp->Cell(i,j)][0]*tileb[fp->Cell((i)-1,(j))][2] > 0 &&   \
	  (fp->tube->glue)[tileb[fp->Cell(i,j)][0]][tileb[fp->Cell((i)-1,j)][2]] < min_strength) ) ? 1 : 0 )

#define getcolor(i,j) ( (fp->Cell(i,j)==0)? translate[0] : (              \
	 (err==1) ? ( errortile(i,j) ? errorcolor : goodcolor ) : (       \
	    (err==2) ? ( (fp->Cell(i,j)>fp->N/2) ?                           \
	       ( errortile(i,j) ? hydroerrcolor : hydrocolor ) :   \
	       ( errortile(i,j) ? errorcolor : goodcolor ) ) :     \
	    translate[fp->Cell(i,j)] ) ) )

/* NOTE: requires 2^P < NCOLS+2*NBDY */
void showpic(flake *fp, int err) /* display the field */  // err param is ignored!
{int row,col,i1,i2,color,j,j1,j2,blocktop=block;
   char *picture=(*spinimage).data; static int last_display_type=0;
   int new_display=(last_display_type!=(10*errors+errorc)); 
   last_display_type=(10*errors+errorc); // re-draw everything when colormap changes

   if (block>4) blocktop=block-1;  
   if (8==(*spinimage).depth) {
      if (block>1) /* I wish I knew how to do this faster */
	 for (row=0;row<size;row++)
	    for (col=0;col<size;col++) {
	       color = getcolor(row,col);
	       j=block*((col+NBDY)+block*NCOLS*(row+NBDY));
	       if (color!=picture[j] || new_display) {
		  for (i1=0;i1<blocktop;i1++) {
		     j1=i1*block*NCOLS+j;
		     for (i2=0;i2<blocktop;i2++)
			picture[j1+i2]=color;
		  }
		  if (block>4) {
		     int n=fp->Cell(row,col);
		     int Ccolm=lightcolor,Ccolp=lightcolor,Crowp=lightcolor,Crowm=lightcolor;
		     if (errors==1) {
			int ncolm=fp->Cell(row,col-1),  ncolp=fp->Cell(row,col+1);  
			int nrowm=fp->Cell(row-1,col),  nrowp=fp->Cell(row+1,col);  
			Ccolm=weakcolor; Ccolp=weakcolor; Crowp=weakcolor; Crowm=weakcolor;
			if (fp->tube->Gse_EW[ n ] [ ncolm ] > 1.5*Gse) Ccolm=strongcolor;
			if (fp->tube->Gse_EW[ ncolp ] [ n ] > 1.5*Gse) Ccolp=strongcolor;
			if (fp->tube->Gse_NS[ n ] [ nrowp ] > 1.5*Gse) Crowp=strongcolor;
			if (fp->tube->Gse_NS[ nrowm ] [ n ] > 1.5*Gse) Crowm=strongcolor;
			if (fp->tube->Gse_EW[ n ] [ ncolm ] < 0.5*Gse) Ccolm=nullcolor;
			if (fp->tube->Gse_EW[ ncolp ] [ n ] < 0.5*Gse) Ccolp=nullcolor;
			if (fp->tube->Gse_NS[ n ] [ nrowp ] < 0.5*Gse) Crowp=nullcolor;
			if (fp->tube->Gse_NS[ nrowm ] [ n ] < 0.5*Gse) Crowm=nullcolor;
		     }
		     if (n==0 && fp->Cell(row,col-1)==0) Ccolm=translate[0];
		     if (n==0 && fp->Cell(row,col+1)==0) Ccolp=translate[0];
		     if (n==0 && fp->Cell(row-1,col)==0) Crowm=translate[0];
		     if (n==0 && fp->Cell(row+1,col)==0) Crowp=translate[0];
		     for (i1=0;i1<blocktop;i1++) // col-1 side
			picture[i1*block*NCOLS+j-1]=Ccolm;
		     for (i1=0;i1<blocktop;i1++) // col+1 side
			picture[i1*block*NCOLS+j+blocktop]=Ccolp;
		     for (i1=0;i1<blocktop;i1++) // row+1 side
			picture[blocktop*block*NCOLS+j+i1]=Crowp;
		     for (i1=0;i1<blocktop;i1++) // row-1 side
			picture[-block*NCOLS+j+i1]=Crowm;
		  }
	       }
	    }
      else { 
	 for (row=0,j=NCOLS*NBDY+NBDY;row<size;row++,j+=2*NBDY) 
	    for (col=0;col<size;col++,j++)
	       picture[j]=getcolor(row,col);
      }
   } else {/* depth is not == 8, use xputpixel (this is really ugly) */
      if (block>1) /* I wish I knew how to do this faster */
	 for (row=0;row<size;row++)
	    for (col=0;col<size;col++) {
	       color=getcolor(row,col);
	       if (color!=XGetPixel(spinimage,j1=block*(col+NBDY),j2=block*(row+NBDY)) 
		     || new_display) {
		  for (i2=0;i2<blocktop;i2++)
		     for (i1=0;i1<blocktop;i1++)
			XPutPixel(spinimage,j1+i1,j2+i2,color);
		  if (block>4) {
		     for (i1=0;i1<blocktop;i1++) {
			int n=fp->Cell(row,col);
			int Ccolm=lightcolor,Ccolp=lightcolor,Crowp=lightcolor,Crowm=lightcolor;
			if (errors==1) {
			   int ncolm=fp->Cell(row,col-1),  ncolp=fp->Cell(row,col+1);  
			   int nrowm=fp->Cell(row-1,col),  nrowp=fp->Cell(row+1,col);  
			   Ccolm=weakcolor; Ccolp=weakcolor; Crowp=weakcolor; Crowm=weakcolor;
			   if (fp->tube->Gse_EW[ n ] [ ncolm ] > 1.5*Gse) Ccolm=strongcolor;
			   if (fp->tube->Gse_EW[ ncolp ] [ n ] > 1.5*Gse) Ccolp=strongcolor;
			   if (fp->tube->Gse_NS[ n ] [ nrowp ] > 1.5*Gse) Crowp=strongcolor;
			   if (fp->tube->Gse_NS[ nrowm ] [ n ] > 1.5*Gse) Crowm=strongcolor;
			   if (fp->tube->Gse_EW[ n ] [ ncolm ] < 0.5*Gse) Ccolm=nullcolor;
			   if (fp->tube->Gse_EW[ ncolp ] [ n ] < 0.5*Gse) Ccolp=nullcolor;
			   if (fp->tube->Gse_NS[ n ] [ nrowp ] < 0.5*Gse) Crowp=nullcolor;
			   if (fp->tube->Gse_NS[ nrowm ] [ n ] < 0.5*Gse) Crowm=nullcolor;
			}
			if (n==0 && fp->Cell(row,col-1)==0) Ccolm=translate[0];
			if (n==0 && fp->Cell(row,col+1)==0) Ccolp=translate[0];
			if (n==0 && fp->Cell(row-1,col)==0) Crowm=translate[0];
			if (n==0 && fp->Cell(row+1,col)==0) Crowp=translate[0];
			XPutPixel(spinimage,j1-1,j2+i1,Ccolm); // col-1 side
			XPutPixel(spinimage,j1+blocktop,j2+i1,Ccolp); // col+1 side
			XPutPixel(spinimage,j1+i1,j2+blocktop,Crowp); // row+1 side
			XPutPixel(spinimage,j1+i1,j2-1,Crowm); // row-1 side
		     }
		  }
	       }
	    }
      else
	 for (row=0;row<size;row++)
	    for (col=0;col<size;col++) {
	       color = getcolor(row,col);
	       XPutPixel(spinimage,col+NBDY,row+NBDY,color);
	    }
   }
   XPutImage(display,playground,gc,spinimage,0,0,0,0,block*NCOLS,block*NROWS); 
   return;
}

#define getcolordij(row,col,di,dj) \
   (((row+(di))>=0 && (row+(di))<size && (col+(dj))>=0 && (col+(dj))<size) ? \
    getcolor(row+di,col+dj) : 0)

/* NOTE: requires 2^P < NCOLS+2*NBDY */
void add_sample_pic(flake *fp, int err) /* add sample the field */
{int row,col,i1,i2,di=0,ddi,j1,j2,dj=0,ddj,blocktop=block, m;
   int color, oldcolor, n_tries=10000,n, collision=0, anything=0;
   m= 2*(err>1); // used in getcolor macro

   // we will try n_tries random offsets before giving up on 
   // avoiding collision.
   // this gets a di, dj
   XDrawImageString(display,samplebutton,gcr,0,font_height,"  SAMPLING   ",13);
   XPutImage(display,playground,gc,spinimage,0,0,0,0,block*NCOLS,block*NROWS); 

   for (n=1; n<n_tries && (collision==1 || anything==0); n++) {
      di= random()%(2*size-1) - size;
      dj= random()%(2*size-1) - size;
      //   printf("trying to place flake %d at %d %d\n", fp, di, dj);
      collision=0; anything=0;
      for (row=0;row<size;row++)
	 for (col=0;col<size;col++) { 
	    color = 0;
	    for (ddi=-1; ddi<=1; ddi++) for (ddj=-1; ddj<=1; ddj++)
	       color = MAX( getcolordij(row,col,di+ddi,dj+ddj), color ); 
	    oldcolor=XGetPixel(spinimage,j1=block*(col+NBDY),j2=block*(row+NBDY));
	    if (oldcolor>0 && color>0) collision=1;   
	    if (getcolordij(row,col,di,dj) >0) anything++;
	 }
      anything = anything && (anything==fp->tiles);
   }
   // if (n==n_tries)  printf("---GAVE UP---\n"); else printf("!!!SUCCESS!!!\n");

   if ( collision==0 ) {
      if (block>4) blocktop=block-1;  
      // always do this the slow way!
      for (row=0;row<size;row++)
	 for (col=0;col<size;col++) { 
	    oldcolor=XGetPixel(spinimage,j1=block*(col+NBDY),j2=block*(row+NBDY));
	    color = MAX( getcolordij(row,col,di,dj), oldcolor ); 
	    if (color!=oldcolor) {
	       for (i2=0;i2<blocktop;i2++)
		  for (i1=0;i1<blocktop;i1++)
		     XPutPixel(spinimage,j1+i1,j2+i2,color);
	       if (block>4) {
		  for (i1=0;i1<blocktop;i1++) {
		     XPutPixel(spinimage,j1-1,j2+i1,lightcolor); // col-1 side
		     XPutPixel(spinimage,j1+blocktop,j2+i1,lightcolor); // col+1 side
		     XPutPixel(spinimage,j1+i1,j2+blocktop,lightcolor); // row+1 side
		     XPutPixel(spinimage,j1+i1,j2-1,lightcolor); // row-1 side
		  }
	       }
	    }
	 }
   }
   XPutImage(display,playground,gc,spinimage,0,0,0,0,block*NCOLS,block*NROWS); 
   XDrawImageString(display,samplebutton,gcr,0,font_height,"   sample    ",13);
   return;
}


/* NOTE: requires 2^P < NCOLS+2*NBDY */
void showphase() /* replace tiles by phase diagram T=1 & T=2 lines */
{int row,col,color,i1,i2,j1,j2,blocktop=block;
   if (block>4) blocktop=block-1;
   for (row=0;row<size;row++)
      for (col=0;col<size;col++) {
	 color = ((size-row)/2==col || (size-row)==col) ?
	    translate[7]:translate[0];
	 j1=block*(col+NBDY); j2=block*(row+NBDY);
	 for (i2=0;i2<blocktop;i2++)
	    for (i1=0;i1<blocktop;i1++)
	       XPutPixel(spinimage,j1+i1,j2+i2,color);
      }
   XPutImage(display,playground,gc,spinimage,0,0,0,0,block*NCOLS,block*NROWS); 
   return;
}


/* fix up the pause button */
void setpause(int value)
{paused=value; if (paused==0) sampling=0;
   if (paused) 
      XDrawImageString(display,pausebutton,gcr,0,font_height,"  run/PAUSE  ",13);
   else
      XDrawImageString(display,pausebutton,gcr,0,font_height,"  RUN/pause  ",13);
}

/* fix up the colors button */
void settilecolor(int value)
{errorc=value;
   if (hydro) {
      if (errorc==2) 
	 XDrawImageString(display,colorbutton,gcr,0,font_height,"tile/err/HYD",12);
      else if (errorc==1) 
	 XDrawImageString(display,colorbutton,gcr,0,font_height,"tile/ERR/hyd",12);
      else if (errorc==0)
	 XDrawImageString(display,colorbutton,gcr,0,font_height,"TILE/err/hyd",12);
   } else {
      if (errorc==1) 
	 XDrawImageString(display,colorbutton,gcr,0,font_height,"  tile/ERR  ",12);
      else if (errorc==0)
	 XDrawImageString(display,colorbutton,gcr,0,font_height,"  TILE/err  ",12);
   }
}

void setsidecolor(int value)
{errors=value;
   if (errors==1) 
      XDrawImageString(display,sidebutton,gcr,0,font_height,"  box/BONDS ",12);
   else if (errors==0)
      XDrawImageString(display,sidebutton,gcr,0,font_height,"  BOX/bonds ",12);
}

/* fix up the export button */
void setexport(int value)
{export_mode=value;
   if (export_mode==2) 
      if (export_movie==0)
	 XDrawImageString(display,exportbutton,gcr,0,font_height,"export[MOVIE]",13);
      else
	 switch (random()%3) {
	    case 0:
	       XDrawImageString(display,exportbutton,gcr,0,font_height,"ExPoRt[MOVIE]",13);
	       break;
	    case 1:
	       XDrawImageString(display,exportbutton,gcr,0,font_height,"eXpOrT[MOVIE]",13);
	       break;
	    case 2:
	       XDrawImageString(display,exportbutton,gcr,0,font_height,"EXPORT[MOVIE]",13);
	       break;
	 }
   else if (export_mode==1) 
      XDrawImageString(display,exportbutton,gcr,0,font_height,"export [ALL] ",13);
   else if (export_mode==0)
      XDrawImageString(display,exportbutton,gcr,0,font_height,"export [ONE] ",13);
}

/* fix up the seed mode button */
void setwander(int value)
{wander=value;
   if (wander) 
      XDrawImageString(display,seedbutton,gcr,0,font_height,"fixed/WANDER",12);
   else
      XDrawImageString(display,seedbutton,gcr,0,font_height,"FIXED/wander",12);
}

/* fix up the fission mode button */
void setfission(int value)
{fission_allowed=value;
   if (hydro && value==2) fission_allowed=0;
   if (fission_allowed>0) 
      if (fission_allowed==1)
	 XDrawImageString(display,fissionbutton,gcr,0,font_height," fission  OK ",13);
      else
	 XDrawImageString(display,fissionbutton,gcr,0,font_height,"chunk fission",13);
   else
      XDrawImageString(display,fissionbutton,gcr,0,font_height," NO  fission ",13);
}


/* this fixes the window up whenever it is uncovered */
void repaint()
{int i=0;
   XDrawString(display,quitbutton,    gcr,0,font_height,"    quit     ",13);
   XDrawString(display,restartbutton, gcr,0,font_height,"   restart   ",13);
   XDrawString(display,cleanbutton,   gcr,0,font_height,"clean/fill/Rx",13);
   XDrawString(display,samplebutton,  gcr,0,font_height,"   sample    ",13);
   XDrawString(display,flakebutton,   gcr,0,font_height,"next/big/prev",13);
   XDrawString(display,tempbutton,    gcr,0,font_height," cool   heat ",13);
   setexport(export_mode);
   setpause(paused);
   settilecolor(errorc); setsidecolor(errors);
   setwander(wander);
   setfission(fission_allowed);

   /* write various strings */
   if (fp) {
      sprintf(stringbuffer,"flake %d (%d by %d%s, seed %d @ (%d,%d))",
	    fp->flake_ID, (1<<fp->P),(1<<fp->P), periodic?", periodic":"",
	    fp->seed_n, fp->seed_i, fp->seed_j);
      XDrawImageString(display,window,gc,5,(++i)*font_height,
	    stringbuffer,strlen(stringbuffer));
      sprintf(stringbuffer,"%lld events, %d tiles, %d mismatches",
	    fp->events, fp->tiles, fp->mismatches);
      XDrawImageString(display,window,gc,5,(++i)*font_height,
	    stringbuffer,strlen(stringbuffer));
   }

   sprintf(stringbuffer,"([DX] = %g uM, T = %5.3f C, 5-mer s.e.)    ",
	 1000000.0*20.0*exp(-Gmc),  4000/(Gse/5+11)-273.15);
   XDrawImageString(display,window,gc,5,(++i)*font_height,
	 stringbuffer,strlen(stringbuffer));

   if (T>0) 
      sprintf(stringbuffer,"Gmc=%4.1f  Gse=%4.1f  k=%6.0f   T=%4.1f",
	    Gmc,tp->Gse,ratek,T/tp->Gse);
   else 
      sprintf(stringbuffer,"Gmc=%4.1f  Gse=%4.1f  k=%6.0f   ",
	    Gmc,tp->Gse,ratek);
   XDrawImageString(display,window,gc,5,(++i)*font_height,
	 stringbuffer,strlen(stringbuffer));

   if (tp->hydro) {
      sprintf(stringbuffer,"Gmch=%4.1f  Gseh=%4.1f  Ghyd=%4.1f",
	    Gmch,Gseh,Ghyd);
      XDrawString(display,window,gc,5,(++i)*font_height,
	    stringbuffer,strlen(stringbuffer));
      sprintf(stringbuffer,"Gas=%4.1f Gam=%4.1f Gae=%4.1f Gah=%4.1f Gao=%4.1f",
	    Gas,Gam,Gae,Gah,Gao);
      XDrawString(display,window,gc,5,(++i)*font_height,
	    stringbuffer,strlen(stringbuffer));
   }
   if (fp) {
      sprintf(stringbuffer,"t = %12.3f sec; G = %12.3f      ",tp->t, 
	    ( wander || tp->conc[fp->seed_n]==0 ) ? fp->G : (fp->G+log(tp->conc[fp->seed_n])));           
   }
   XDrawImageString(display,window,gc,5,(++i)*font_height,
	 stringbuffer,strlen(stringbuffer));
   sprintf(stringbuffer, "%lld events (%llda,%lldd,%lldh,%lldf), %lld tiles total %s      ",
	 tp->events, tp->stat_a, tp->stat_d, tp->stat_h, tp->stat_f,
	 tp->stat_a-tp->stat_d+tp->num_flakes, tp->ewrapped?"[wrapped]":"");
   XDrawImageString(display,window,gc,5,(++i)*font_height,
	 stringbuffer,strlen(stringbuffer));

   if (fp && fp->flake_conc>0) { int tl;
      sprintf(stringbuffer, "Gfc=%4.1f; Gmc=[", -log(fp->flake_conc)); 
      for (tl=1; tl<=tp->N && strlen(stringbuffer)<230; tl++)
	 sprintf(stringbuffer+strlen(stringbuffer)," %4.1f",
	       (tp->conc[tl]>0)?-log(tp->conc[tl]):0 );
      sprintf(stringbuffer+strlen(stringbuffer)," ]");
      XDrawImageString(display,window,gc,5,(++i)*font_height,
	    stringbuffer,strlen(stringbuffer));
   }

   XDrawString(display,window,gc,WINDOWWIDTH-150,WINDOWHEIGHT-65
	 ," left: identify",15); 
   XDrawString(display,window,gc,WINDOWWIDTH-150,WINDOWHEIGHT-45
	 ,"middle:puncture",15); 
   XDrawString(display,window,gc,WINDOWWIDTH-150,WINDOWHEIGHT-25
	 ,"right: Gmc Gse ",15); 

   XDrawString(display,window,gc,WINDOWWIDTH-190,WINDOWHEIGHT-5
	 ,"EW, RS, CGE '98-'15",19); 

   if (!sampling && fp) showpic(fp,errorc); 
   else XPutImage(display,playground,gc,spinimage,0,0,0,0,block*NCOLS,block*NROWS); 
}

/* a lot of this is taken from the basicwin program in the
   Xlib Programming Manual */
void openwindow(int argc, char **argv)
{char *window_name;
   char *icon_name="xgrow";
   Pixmap icon_pixmap;
   char *display_name=NULL;
   XEvent report;
   XColor xcolor,colorcell;
   Colormap cmap;
   char *buffer;
   int i,j;
# define icon_bitmap_width 16
# define icon_bitmap_height 16
   static char icon_bitmap_bits[] = {
      0x1f, 0xf8, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0xf8,
      0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};

   window_name=(char *) malloc(256);
   sprintf(window_name,"xgrow: %s",tileset_name);

   /* open up the display */
   if ((display=XOpenDisplay(display_name))==NULL)
   {fprintf(stderr,"%s: cannot connect to X server %s\n",
	 progname,XDisplayName(display_name));
   exit(-1);
   }

   screen=DefaultScreen(display);
   depth=DefaultDepth(display,screen);
   cmap=DefaultColormap(display,screen);
   /* color? This is not the right way to do it, but .... */
   if (1==depth) 
   {fprintf(stderr,"Sorry but this program needs a color monitor.\n");
      exit(-1);
   }
   black=BlackPixel(display,screen);
   white=WhitePixel(display,screen);
   if (XAllocNamedColor(display,cmap,"firebrick",&colorcell,&xcolor))
      darkcolor=colorcell.pixel;
   if (XAllocNamedColor(display,cmap,"wheat",&colorcell,&xcolor))
      lightcolor=colorcell.pixel;
   if (XAllocNamedColor(display,cmap,"red",&colorcell,&xcolor))
      nullcolor=colorcell.pixel;
   if (XAllocNamedColor(display,cmap,"wheat",&colorcell,&xcolor))
      weakcolor=colorcell.pixel;
   if (XAllocNamedColor(display,cmap,"green",&colorcell,&xcolor))
      strongcolor=colorcell.pixel;
   if (XAllocNamedColor(display,cmap,"red",&colorcell,&xcolor))
      errorcolor=colorcell.pixel;
   if (XAllocNamedColor(display,cmap,"green",&colorcell,&xcolor))
      goodcolor=colorcell.pixel;
   if (XAllocNamedColor(display,cmap,"gold",&colorcell,&xcolor))
      hydrocolor=colorcell.pixel;
   if (XAllocNamedColor(display,cmap,"pink",&colorcell,&xcolor))
      hydroerrcolor=colorcell.pixel;
   if (XAllocNamedColor(display,cmap,"black",&colorcell,&xcolor))
      translate[0]=colorcell.pixel;

   for (i=1;i<MAXTILETYPES;i++) 
      if (tile_colors[i]!=NULL &&
	    XAllocNamedColor(display,cmap,tile_colors[i],&colorcell,&xcolor))
	 translate[i]=colorcell.pixel;
      else translate[i]=translate[(i>14) ? ((i-1)%14+1) : 0]; 

   /* make the main window */
   window=XCreateSimpleWindow(display,RootWindow(display,screen),
	 100,0,WINDOWWIDTH,WINDOWHEIGHT,4,black,lightcolor);

   /* make the icon */
   icon_pixmap=XCreateBitmapFromData(display,window,
	 icon_bitmap_bits,icon_bitmap_width,icon_bitmap_height);

   size_hints.flags=PPosition | PSize | PMinSize;
   size_hints.min_width=WINDOWWIDTH;
   size_hints.min_height=WINDOWHEIGHT;
#ifdef X11R3
   size_hints.x=x;
   size_hints.y=y;
   size_hints.width=WINDOWWIDTH;
   size_hints.height=WINDOWHEIGHT;
   XSetStandardProperties(display,window,window_name,icon_name,
	 icon_pixmap,argv,argc,&size_hints);
#else
   {XWMHints wm_hints;
      XClassHint class_hints;
      XTextProperty windowName, iconName;
      if (XStringListToTextProperty(&window_name,1,&windowName)==0)
      {fprintf(stderr,"%s: structure allocation for windowName failed.\n"
	    ,progname);
      exit(-1);
      }
      if (XStringListToTextProperty(&icon_name,1,&iconName)==0)
      {fprintf(stderr,"%s: structure allocation for iconName failed.\n"
	    ,progname);
      exit(-1);
      }
      wm_hints.initial_state=NormalState;
      wm_hints.input=True;
      wm_hints.icon_pixmap=icon_pixmap;
      wm_hints.flags=StateHint|IconPixmapHint|InputHint;
      class_hints.res_name=progname;
      class_hints.res_class="Basicwin";
      XSetWMProperties(display,window,&windowName,&iconName,
	    argv,argc,&size_hints,&wm_hints,&class_hints);
   }
#endif


   /* make the buttons */
   quitbutton=XCreateSimpleWindow(display,window,
	 WINDOWWIDTH-140,WINDOWHEIGHT-124,120,20,2,black,darkcolor);
   restartbutton=XCreateSimpleWindow(display,window,
	 WINDOWWIDTH-140,WINDOWHEIGHT-150,120,20,2,black,darkcolor);
   cleanbutton=XCreateSimpleWindow(display,window,
	 WINDOWWIDTH-140,WINDOWHEIGHT-176,120,20,2,black,darkcolor);
   seedbutton=XCreateSimpleWindow(display,window,
	 WINDOWWIDTH-140,WINDOWHEIGHT-202,120,20,2,black,darkcolor);
   fissionbutton=XCreateSimpleWindow(display,window,
	 WINDOWWIDTH-140,WINDOWHEIGHT-228,120,20,2,black,darkcolor);
   tempbutton=XCreateSimpleWindow(display,window,
	 WINDOWWIDTH-140,WINDOWHEIGHT-254,120,20,2,black,darkcolor);
   flakebutton=XCreateSimpleWindow(display,window,
	 WINDOWWIDTH-140,WINDOWHEIGHT-280,120,20,2,black,darkcolor);
   samplebutton=XCreateSimpleWindow(display,window,
	 WINDOWWIDTH-140,WINDOWHEIGHT-306,120,20,2,black,darkcolor);
   exportbutton=XCreateSimpleWindow(display,window,
	 WINDOWWIDTH-140,WINDOWHEIGHT-332,120,20,2,black,darkcolor);
   pausebutton=XCreateSimpleWindow(display,window,
	 WINDOWWIDTH-140,WINDOWHEIGHT-358,120,20,2,black,darkcolor);
   colorbutton=XCreateSimpleWindow(display,window,
	 WINDOWWIDTH-140,WINDOWHEIGHT-384,120,20,2,black,darkcolor);
   sidebutton=XCreateSimpleWindow(display,window,
	 WINDOWWIDTH-140,WINDOWHEIGHT-410,120,20,2,black,darkcolor);
   playground=XCreateSimpleWindow(display,window,
	 PLAYLEFT,PLAYTOP,block*NCOLS,block*NROWS,2,translate[4],white);

   /* pick the events to look for */
   event_mask=ExposureMask|ButtonPressMask|StructureNotifyMask;
   XSelectInput(display,window,event_mask);
   event_mask=ButtonPressMask; 
   /* note that with this simple mask if one just covers a button
      it will not get redrawn.  I wonder if anyone will notice?  If I put
      the exposuremask in here, things flash irritatingly on being uncovered. */
   XSelectInput(display,quitbutton,event_mask);
   XSelectInput(display,pausebutton,event_mask);
   XSelectInput(display,restartbutton,event_mask);
   XSelectInput(display,cleanbutton,event_mask);
   XSelectInput(display,colorbutton,event_mask);
   XSelectInput(display,sidebutton,event_mask);
   XSelectInput(display,flakebutton,event_mask);
   XSelectInput(display,seedbutton,event_mask);
   XSelectInput(display,fissionbutton,event_mask);
   XSelectInput(display,tempbutton,event_mask);
   XSelectInput(display,samplebutton,event_mask);
   XSelectInput(display,exportbutton,event_mask);
   event_mask=ButtonReleaseMask|ButtonPressMask|PointerMotionHintMask
      |ButtonMotionMask; 
   XSelectInput(display,playground,event_mask);

   /* pick font: 9x15 is supposed to almost always be there */
   if ((font=XLoadQueryFont(display,"9x15"))==NULL) {
      if ((font=XLoadQueryFont(display,"6x13"))==NULL)
      {
	 fprintf(stderr,"%s: Cannot open 9x15 or 6x13 font\n",progname);
	 exit(-1);
      }
   }
   font_height=font->ascent+font->descent;

   /* make graphics contexts: 
      gc for black on white (actually, lightcolor)
      gccolor for background and buttons 
      gcr for reverse video  */

   gc=XCreateGC(display,window,0,NULL);
   XSetFont(display,gc,font->fid);
   XSetForeground(display,gc,black);
   XSetBackground(display,gc,lightcolor); 
   /* speed up? */
   XSetPlaneMask(display,gc,black|white|translate[0]|translate[1]|
	 translate[2]|translate[3]|translate[16]|translate[4]);

   gcr=XCreateGC(display,window,0,NULL); 
   XSetFont(display,gcr,font->fid);
   XSetForeground(display,gcr,lightcolor);
   XSetBackground(display,gcr,darkcolor);

   gccolor=XCreateGC(display,window,0,NULL); 
   XSetFont(display,gccolor,font->fid);
   XSetForeground(display,gccolor,darkcolor);
   XSetBackground(display,gccolor,lightcolor); 

   /* show the window and buttons */
   XMapWindow(display,window);
   XMapWindow(display,quitbutton);
   XMapWindow(display,pausebutton);
   XMapWindow(display,restartbutton);
   XMapWindow(display,cleanbutton);
   XMapWindow(display,colorbutton);
   if (block>4) XMapWindow(display,sidebutton);
   if (~(fparam->N==1 && fparam->next_param==NULL))  
      XMapWindow(display,flakebutton); 
   XMapWindow(display,seedbutton);
   XMapWindow(display,fissionbutton);
   XMapWindow(display,tempbutton);
   XMapWindow(display,samplebutton);
   XMapWindow(display,exportbutton);
   XMapWindow(display,playground);


   /* make image structure */
   /* wait for playground to be displayed before proceeding */
   i=1; /* a flag */
   while (i)
   {XNextEvent(display,&report); 
      switch (report.type)
      {case Expose:
	 if (report.xexpose.window!=playground) i=0;
	 default:
	 break;
      }
   }

   /* FIXME: I Don't know why 6 is the right value here!? 8 doesn't work */
   buffer = malloc( block*block*NROWS*NCOLS*depth/6 );
   spinimage=XCreateImage(display, CopyFromParent, depth, ZPixmap, 0, buffer,  block*NCOLS, block*NROWS, 32, 0);
   XInitImage(spinimage);

   if (NULL==spinimage)
   {fprintf(stderr,"trouble creating image structure\n");
      exit(-1);
   } 

   /* make sure everything get written first time */
   for (i=0;i<block*NROWS;i++) 
      for (j=0;j<block*NCOLS;j++)
	 XPutPixel(spinimage,i,j,translate[0]);

}

void cleanup()
{XUnloadFont(display,font->fid);
   XFreeGC(display,gc); 
   XFreeGC(display,gcr); 
   XFreeGC(display,gccolor); 
   XCloseDisplay(display);
   XDestroyImage(spinimage); 
   exit(1);
} 


int main(int argc, char **argv)
{
   int x,y,b,i,j;    int clear_x=0,clear_y=0;
   int mousing=0; int stat=0;
   double new_Gse, new_Gmc;
   XEvent report;
   progname=argv[0];

   for (i=0;i<MAXTILETYPES;i++) {
      translate[i]=0; /* translate[] from tiles to colors */
   }

   getargs(argc, argv);

   if (import_flake_size > size)
   {
      fprintf(stderr, "Error: The input flake size is larger than the "
	    "simulation flake size.\n");
      exit(-1);
   }

   if (hydro) { /* automatically double the number of tiles */
      tileb=(int**) realloc(tileb,sizeof(int*)*(2*N+1));
      for (i=1;i<=N;i++) {
	 tileb[i+N]= (int*) calloc(sizeof(int),4);
	 for (j=0;j<4;j++) tileb[i+N][j] = tileb[i][j];
      }
      stoic=(double*) realloc(stoic,sizeof(double)*(2*N+1));
      for (i=1;i<=N;i++) stoic[i+N]= stoic[i];
      N = N*2;
   }

   if (linear) {
      linear_simulate(ratek,Gmc,Gse,tmax,emax,smax,mmax);
      return 0;
   }

   if (testing) {
      tp = init_tube(size_P,N,num_bindings);   
      set_params(tp,tileb,strength,glue,stoic,0,initial_rc,updates_per_RC,
	    anneal_h,anneal_s,startC,endC,seconds_per_C,
	    dt_right, dt_left, dt_down, dt_up, hydro,ratek,
	    Gmc,Gse,Gmch,Gseh,Ghyd,Gas,Gam,Gae,Gah,Gao,T,tinybox, seed_i, seed_j, Gfc);
#ifdef TESTING_OK
      run_xgrow_tests(tp,Gmc,Gse,seed_i,seed_j,seed_n,size);
#endif
      return 0;
   }

   if (XXX) openwindow(argc,argv);

   /* printf("xgrow: tile set read, beginning simulation\n"); */

   /* set initial state */
   tp = init_tube(size_P,N,num_bindings);   
   set_params(tp,tileb,strength,glue,stoic,anneal_g,anneal_t,updates_per_RC,anneal_h,anneal_s,startC,endC,seconds_per_C,dt_right, dt_left, dt_down, dt_up, hydro,ratek,
	 Gmc,Gse,Gmch,Gseh,Ghyd,Gas,Gam,Gae,Gah,Gao,T,tinybox,seed_i,seed_j,Gfc);

   fprm=fparam;


   /* initialize flakes */
   while (fprm!=NULL)
   {
      int fn;
      for (fn=1; fn <= fprm->N; fn++)
      {
	 if (tp->dt_left[fprm->seed_n]) {
	    fprm->seed_n = tp->dt_left[fprm->seed_n]; // FIXME: vdouble
	    fprm->seed_j--;
	 }
	 insert_flake(fp=init_flake(size_P,N,
		  fprm->seed_i,fprm->seed_j,fprm->seed_n,fprm->Gfc), tp);
	 if (tp->dt_right[fprm->seed_n]) {
	    change_cell (fp,seed_i,seed_j+1,tp->dt_right[fprm->seed_n]);
	    fp->seed_is_double_tile = 1;
	 }
	 assert (!tp->dt_left[fprm->seed_n]);

	 if (fprm->import_from != NULL)
	 {
	    fprintf(stderr, "WARNING: In imported flakes, the seed position is chosen randomly.\n");
	    import_flake(fp, fprm->import_from, fn);
	 }
      }
      fprm=fprm->next_param;
   }

   //   print_tree(tp->flake_tree,0,'*'); 

   if (stripe_args!=NULL) {
      /* STRIPE OPTION HAS HARDCODED TILESET NONSENSE -- A BUG */
      /* SHOULD REPLACE THIS OPTION BY READING INITIAL FLAKE FROM FILE */
      int i,j,k,w; double p; char *s=stripe_args;
      char XOR[2][2]={ {4,7}, {6,5} }; /* XOR[S][E] */ char c,cc;
      i=size-1; j = atoi(s)%size; 
      for (k=0; k<size; k++) 
	 change_cell(fp, (i-k+size)%size, (j+k)%size, 4+random()%4); 
      fp->seed_i=i; fp->seed_j=j; fp->seed_n=fp->Cell(i,j);
      s = strchr(s,':');
      while (s!=NULL) {
	 p = atof(s+1); s = strchr(s,',');
	 if (s!=NULL) {
	    w = atoi(s+1); s = strchr(s,':');
	    for (;w>0;w--) {
	       i=(i-1+size)%size;
	       for (k=0; k<size; k++) {
		  cc=c= XOR[(fp->Cell((i-k+1+size)%size,(j+k)%size)-4)/2]
		     [(fp->Cell((i-k+size)%size,(j+k+1)%size)-4)/2];
		  if (drand48()<p) do cc=4+random()%4; while (cc==c);
		  change_cell(fp, (i-k+size)%size, (j+k)%size, cc);
		  tp->events--; /* don't count these as events */
		  /* ERROR: stats are also modified!!! */
	       }
	    }
	 }
      }
      /* no corner or boundary tiles for stripe simulations */
      tp->conc[0] -= tp->conc[1]; tp->conc[0]+=(tp->conc[1]=exp(-35));
      tp->conc[0] -= tp->conc[2]; tp->conc[0]+=(tp->conc[2]=exp(-35));
      tp->conc[0] -= tp->conc[3]; tp->conc[0]+=(tp->conc[3]=exp(-35));
   }


   // printf("flake initialized, size_P=%d, size=%d\n",size_P,size);

   new_Gse=Gse; new_Gmc=Gmc;
   if (tracefp!=NULL) write_datalines(tracefp,"\n");


   // printf("tmax=%f  emax=%d  smax=%d\n",tmax,emax,smax);

   if (XXX) repaint();

   /* loop forever, looking for events */
   while((tmax==0 || tp->t < tmax) && 
	 (emax==0 || tp->events < emax) &&
	 (smax==0 || tp->stat_a-tp->stat_d < smax) &&
	 (mmax==0 || tp->stat_m < mmax) &&
	 (smin==-1 || tp->stat_a-tp->stat_d > smin) &&
	 (fsmax==0 || tp->largest_flake_size < fsmax) &&
	 (tp->seconds_per_C == 0 || tp->currentC > tp->endC) &&
	 !(untiltiles && tp->all_present)) { 

      Gse=tp->Gse;  // keep them sync'd in case "anneal" is ongoing.
      if (!XXX) {
	 simulate(tp,update_rate,tmax,emax,smax,fsmax,smin,mmax);
	 if (tracefp!=NULL) write_datalines(tracefp,"\n");
	 if (export_mode==2 && export_movie==1) export_flake("movie",fp);
      } else {
	 if (0==paused && 0==mousing && !XPending(display)) {
	    simulate(tp,update_rate,tmax,emax,smax,fsmax,smin,mmax);
	    fp = tp->flake_list;
	    assert (!fp || !tp->tinybox ||
		  ((!fp->seed_is_double_tile && fp->tiles > 1) || fp->tiles > 2));
	    if (tracefp!=NULL) write_datalines(tracefp,"\n");
	    if (export_mode==2 && export_movie==1) export_flake("movie",fp);
	    if (fp && fp->flake_conc>0) recalc_G(fp);
	    // make sure displayed G is accurate for conc's
	    // hopefully this won't slow things down too much.
	    stat++; if (stat==1) { stat=0; repaint(); }
	 }
	 if (paused|mousing|XPending(display)) {
	    XNextEvent(display,&report); 
	    switch (report.type)
	    {case Expose:
	       if (report.xexpose.count!=0) break; /* more in queue, wait for them */
	       repaint();  
	       break;
	       case ConfigureNotify:
	       break; 
	       case MotionNotify:
	       if (report.xbutton.window==playground) {
		  int newx,newy;
		  Window root,child;
		  unsigned int keys_buttons;
		  int window_x,window_y;
		  if (tp->hydro) break; /* don't know how to reset params */
		  x=report.xbutton.x/block;
		  y=report.xbutton.y/block;
		  b=report.xbutton.button;
		  if (mousing==0) { int i,j,t;
		     i=MIN(MAX(0,y-2),size-1); j=MIN(MAX(0,x-2),size-1); t=fp->Cell(i,j);
		     sprintf(stringbuffer,"([DX] = %g uM, T = %5.3f C, 5-mer s.e.)    "
			   "tile #%d = {%d %d %d %d} at (%d,%d)           ",
			   1000000.0*20.0*exp(-Gmc),  4000/(Gse/5+11)-273.15,
			   t,tp->tileb[t][0],tp->tileb[t][1],tp->tileb[t][2],tp->tileb[t][3],i,j);
		     XDrawImageString(display,window,gc,5,2*font_height,
			   stringbuffer,strlen(stringbuffer));
		  }
		  if (mousing==3) {
		     /* was sketch(x,y,b); now change Gse & Gmc */
		     new_Gse=(30.0*x)/size; new_Gmc=30-(30.0*y)/size;
		     /* draw current Gse, Gmc values */
		     sprintf(stringbuffer,"Gmc=%4.1f->%4.1f  Gse=%4.1f->%4.1f",
			   Gmc,new_Gmc,Gse,new_Gse);
		     XDrawImageString(display,window,gc,5,3*font_height,
			   stringbuffer,strlen(stringbuffer));
		     /* this is necessary in order to get more than one MotionNotify
			event; I don't know why. */
		  }
		  if (!XQueryPointer(display,playground,
			   &root,&child,&window_x,&window_y,
			   &newx,&newy,&keys_buttons))
		  {mousing=0; printf("Weird X feature\n"); }
	       }
	       break;
	       case ButtonRelease:
	       if (mousing==3) { // change Gmc, Gse w/ visual GUI
		  printf("Changing Gmc -> %f, Gse -> %f \n", new_Gmc, new_Gse);
		  reset_params(tp, Gmc, Gse, new_Gmc, new_Gse,Gseh);
		  if (Gfc>0) Gfc+=(new_Gmc-Gmc); 
		  fprm=fparam;
		  while (fprm!=NULL) {  /* fix up all flake info, for restarting */
		     if (fprm->Gfc > 0) fprm->Gfc += (new_Gmc-Gmc); 
		     fprm=fprm->next_param;
		  } 
		  Gse=new_Gse; Gmc=new_Gmc; tp->Gse=Gse; tp->Gmc=Gmc;
		  showpic(fp,errorc); 
	       } else if (mousing==1) {
		  /* clear a region, i.e., "puncture" */
		  int i,j,mi,Mi,mj,Mj;
		  x=report.xbutton.x/block;
		  y=report.xbutton.y/block;
		  b=report.xbutton.button; 
		  if (clear_y>y) { mi=MIN(MAX(0,y-2),size-1); Mi=MIN(MAX(0,clear_y-2),size-1); }
		  else           { mi=MIN(MAX(0,clear_y-2),size-1); Mi=MIN(MAX(0,y-2),size-1); }
		  if (clear_x>x) { mj=MIN(MAX(0,x-2),size-1); Mj=MIN(MAX(0,clear_x-2),size-1); }
		  else           { mj=MIN(MAX(0,clear_x-2),size-1); Mj=MIN(MAX(0,x-2),size-1); }
		  if (mi != Mi && mj != Mj) { int fa=fission_allowed;
		     if (mi<=fp->seed_i && fp->seed_i<=Mi &&
			   mj<=fp->seed_j && fp->seed_j<=Mj) { int si=fp->seed_i,sj=fp->seed_j;
			// better move the seed out of the way, if possible
			for (i=mi; i<=Mi; i++) if (fp->Cell(i,mj-1)>0) { si=i; sj=(mj-1+size)%size; }
			for (i=mi; i<=Mi; i++) if (fp->Cell(i,Mj+1)>0) { si=i; sj=(Mj+1)%size; }
			for (j=mj; j<=Mj; j++) if (fp->Cell(mi-1,j)>0) { si=(mi-1+size)%size; sj=j; }
			for (j=mj; j<=Mj; j++) if (fp->Cell(Mi+1,j)>0) { si=(Mi+1)%size; sj=j; }
			change_seed(fp,si,sj);
		     }
		     fission_allowed=1; 
		     for (i=mi; i<=Mi; i++)
			for (j=mj; j<=Mj; j++) {
			   if (fp->Cell(i,j)>0) { 
			      change_cell(fp, i, j, 0); flake_fission(fp,i,j); 
			   }
			}
		     fission_allowed=fa; 
		  }
		  showpic(fp,errorc); 
	       }
	       mousing=0; 
	       break;
	       case ButtonPress:
	       if (report.xbutton.window==quitbutton) {
		  // force regular immediate exit by setting smax to current size
		  smax=tp->stat_a-tp->stat_d;
		  // (was cleanup(); but this circumvented datafile output, etc.)
	       } else if (report.xbutton.window==pausebutton) {
		  setpause(1-paused); repaint(); 
	       } else if (report.xbutton.window==restartbutton) {
		  free_tube(tp); 
		  tp = init_tube(size_P,N,num_bindings);   
		  set_params(tp,tileb,strength,glue,stoic,anneal_g,anneal_t,updates_per_RC,
			anneal_h, anneal_s, startC, endC, seconds_per_C,
			dt_right, dt_left, dt_down, dt_up, hydro,ratek,
			Gmc,Gse,Gmch,Gseh,Ghyd,Gas,Gam,Gae,Gah,Gao,T,tinybox,seed_i,seed_j,
			Gfc);
		  fprm=fparam; 
           while (fprm!=NULL)
           {
              int fn;
              for (fn=1; fn <= fprm->N; fn++)
              {
             if (tp->dt_left[fprm->seed_n]) {
                fprm->seed_n = tp->dt_left[fprm->seed_n]; // FIXME: vdouble
                fprm->seed_j--;
             }
             insert_flake(fp=init_flake(size_P,N,
                  fprm->seed_i,fprm->seed_j,fprm->seed_n,fprm->Gfc), tp);
             if (tp->dt_right[fprm->seed_n]) {
                change_cell (fp,seed_i,seed_j+1,tp->dt_right[fprm->seed_n]);
                fp->seed_is_double_tile = 1;
             }
             assert (!tp->dt_left[fprm->seed_n]);

             if (fprm->import_from != NULL)
             {
                fprintf(stderr, "WARNING: In imported flakes, the seed position is chosen randomly.\n");
                import_flake(fp, fprm->import_from, fn);
             }
              }
              fprm=fprm->next_param;
           }
		  repaint();
	       } else if (report.xbutton.window==colorbutton) { // show tiles or error or hyd
		  settilecolor(hydro ? (errorc+1)%3 : (errorc+1)%2); repaint(); 
	       } else if (report.xbutton.window==sidebutton) { // show box or null/weak/strong
		  setsidecolor((errors+1)%2); repaint(); 
	       } else if (report.xbutton.window==seedbutton) {
		  setwander(1-wander); repaint(); 
	       } else if (report.xbutton.window==fissionbutton) {
		  fission_allowed = (fission_allowed+1)%3;
		  setfission(fission_allowed); 
		  update_all_rates(tp);  
		  repaint(); 
	       } else if (report.xbutton.window==cleanbutton) {
		  x=report.xbutton.x;
		  y=report.xbutton.y;
		  b=report.xbutton.button;
		  if (x<50) {  // do a clean_flake cycle
		     printf("Cleaning 1 cycle, clean_X=%f\n",clean_X);
		     clean_flake(fp,clean_X,1); 
		  } else if (x<90) { // do a fill_flake cycle
		     printf("Filling 1 cycle, fill_X=%f\n",fill_X);
		     fill_flake(fp,fill_X,1); 
		  } else {
		     printf("Repairing, uniqueness for T=%f\n",repair_unique_T);
		     repair_flake(fp,repair_unique_T,Gse); 
		  }
		  repaint();           
	       } else if (report.xbutton.window==exportbutton) {
		  x=report.xbutton.x;
		  y=report.xbutton.y;
		  b=report.xbutton.button;
		  if (x>70) {  // change from ONE to ALL to MOVIE
		     setexport((export_mode+1)%3); 
		  } else if (export_mode==0) { 
		     // output to file (unless MOVIE mode already)
		     export_flake("flake", fp);
		  } else if (export_mode==1) { 
		     flake *tfp=tp->flake_list; 
		     while (tfp->next_flake != NULL) {
			export_flake("flake", tfp);
			tfp=tfp->next_flake;
		     }
		  } else if (export_mode==2) {
		     // turn movie mode on & off
		     export_movie= !export_movie;
		  }
		  repaint();           
	       } else if (report.xbutton.window==samplebutton) {
		  flake *tfp; int n,num_big=0;
		  x=report.xbutton.x;
		  y=report.xbutton.y;
		  b=report.xbutton.button;
		  // stop simulation if you're sampling.
		  setpause(1);  sampling=1;
		  // pick a (if possible non-monomer) sample and add to the field
		  for (tfp=tp->flake_list; tfp!=NULL; tfp=tfp->next_flake) {
		     if (tfp->tiles>1) num_big++;
		  } 
		  n = random()%(num_big>0 ? num_big : tp->num_flakes);
		  tfp=tp->flake_list;
		  while (num_big>0 && tfp->tiles==1) { tfp=tfp->next_flake; } 
		  for (i=0;i<n;i++) {
		     tfp=tfp->next_flake;
		     while (num_big>0 && tfp->tiles==1) { tfp=tfp->next_flake; } 
		  }
		  add_sample_pic(tfp,errorc); 
		  repaint(); 
	       } else if (report.xbutton.window==flakebutton) { 
		  flake *tfp=tp->flake_list; 
		  x=report.xbutton.x;
		  y=report.xbutton.y;
		  b=report.xbutton.button; 
		  // cycle through flakes.
		  if (x>80) {
		     fp=fp->next_flake;
		     if (fp==NULL) fp=tp->flake_list;
		  } else if (x<40) {
		     while (tfp && tfp->next_flake != NULL && tfp->next_flake != fp) 
			tfp=tfp->next_flake;
		     if (tfp==NULL) fp=tp->flake_list; else fp=tfp;
		  } else {  // find the biggest flake.
		     fp=tfp;
		     if (fp) {
			while (tfp->next_flake != NULL) {
			   tfp=tfp->next_flake;
			   if (tfp->tiles > fp->tiles) fp=tfp; 
			}
		     }
		  } 
		  //             print_tree(tp->flake_tree,0,'*');  
		  sampling=0; repaint();
	       } else if (report.xbutton.window==tempbutton) { // change Gse w/ button
		  x=report.xbutton.x;  
		  y=report.xbutton.y;
		  b=report.xbutton.button;
		  if (tp->hydro) break; /* don't know how to reset params */
		  if (x>60) new_Gse=Gse-0.1; else new_Gse=Gse+0.1;
		  reset_params(tp, Gmc, Gse, new_Gmc, new_Gse,Gseh);
		  Gse=new_Gse; tp->Gse=Gse; repaint();
	       } else if (report.xbutton.window==playground) // we're in ButtonPress
	       {x=report.xbutton.x/block;
		  y=report.xbutton.y/block;
		  b=report.xbutton.button;
		  if (b==3) {
		     if (tp->hydro) break; /* don't know how to reset params */
		     new_Gse=(30.0*x)/size; new_Gmc=30-(30.0*y)/size;
		     /* draw current Gse, Gmc values */
		     sprintf(stringbuffer,"Gmc=%4.1f->%4.1f  Gse=%4.1f->%4.1f",
			   Gmc,new_Gmc,Gse,new_Gse);
		     XDrawImageString(display,window,gc,5,3*font_height,
			   stringbuffer,strlen(stringbuffer));
		     /* later: if down button, draw T=1/T=2 diagram */
		     mousing=3; showphase();
		  } else if (b==2) { // "puncture"
		     mousing=1; clear_x=x; clear_y=y;  /* prepare to clear a region */
		  } else if (b==1) { // "identify"
		     int i,j,t;
		     i=MIN(MAX(0,y-2),size-1); j=MIN(MAX(0,x-2),size-1); t=fp->Cell(i,j);
		     sprintf(stringbuffer,"([DX] = %g uM, T = %5.3f C, 5-mer s.e.)    "
			   "tile #%d = {%d %d %d %d} at (%d,%d)           ",
			   1000000.0*20.0*exp(-Gmc),  4000/(Gse/5+11)-273.15,
			   t,tp->tileb[t][0],tp->tileb[t][1],tp->tileb[t][2],tp->tileb[t][3],i,j);
		     XDrawImageString(display,window,gc,5,2*font_height,
			   stringbuffer,strlen(stringbuffer));
		  }
	       }
	       else {
		  /*  what here ?  */
	       }
	       break;
	       default:
	       break;
	    } /* end of switch */
	 } /* end of if XPending */
      }} /* end of while(...) if...else */

   closeargs();
   return 0;
} /* end of main */



