/*  xgrow
           
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

  Compiling:  see rosemake
    
*/

# include <X11/Xlib.h>
# include <X11/Xutil.h>
# include <X11/Xos.h>
# include <X11/Xatom.h> 
# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <unistd.h>
# include <malloc.h>
# include <math.h>
# include <assert.h>

# include "grow.h"

#define NUPDATES 10000
#define NSTATS   1

 /* lattice dimensions (plus two for boundaries): */
 /* NCOLS should be a multiple of bytes per long word (BPW) */
         /* THIS is for X bitmap stuff, so add total 4 for boundary */
# define NBDY 2
# define NROWS (256+NBDY*2)
# define NCOLS (256+NBDY*2)
# define VOLUME (NROWS*NCOLS)

# define WINDOWWIDTH (block*NCOLS+PLAYLEFT+BOUNDWIDTH+20)
# define WINDOWHEIGHT (PLAYTOP+block*NROWS+10)
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
# define MAXTILETYPES 256

long int translate[MAXTILETYPES]; /* for converting colors */
int paused=0, errorc=0;
static char *progname;
char stringbuffer[100];

/* various window stuff */
Display *display;
int screen;
Window window,quitbutton,pausebutton,playground,
        fillbutton,colorbutton,randombutton,seedbutton;
GC gc, gcr, gccolor;
XImage *spinimage=NULL;
XFontStruct *font=NULL;
int font_height;
XSizeHints size_hints;
long event_mask;
int depth;
long int darkcolor,lightcolor,black,white;


/* simulation parameters */
flake *fp; double Gse, Gmc, ratek;
double Gseh, Gmch, Ghyd, Gas, Gam, Gae, Gah, Gao, Gfc;

int block=1; /* default to small blocks; calling with argument changes this */
int wander, periodic, deplete, linear; 
FILE *datafp, *arrayfp, *tilefp;
double* strength;
int N, num_bindings, units_length;
int** units; double* stoic;
int hydro;
double tmax; int emax, smax;
int seed_i,seed_j,seed_n;
char *stripe_args=NULL;
int XXX=1;


void read_tilefile(FILE *tilefp) 
{ 
 float strength_float, stoic_float; int i,j;

 fscanf(tilefp,"tile edges matches {{N E S W}*}\n");

 fscanf(tilefp,"num tile types=%d\n",&N);
 fscanf(tilefp,"num binding types=%d\n", &num_bindings);

 fscanf(tilefp,"tile edges=");
 units_length = N+1;
 units = (int**) calloc(sizeof(int*),units_length);
 stoic = (double*) calloc(sizeof(double),units_length);
 fscanf(tilefp,"{\n");
 units[0] = (int*) calloc(sizeof(int),4);
 for (j=0;j<4;j++) {
   units[0][j] = 0;
 }
 stoic[0]=0;
 for (i=1;i<units_length;i++) {
   units[i] = (int*) calloc(sizeof(int),4);
   fscanf(tilefp,"{");
   for (j=0;j<4;j++) {
     fscanf(tilefp,"%d",&units[i][j]);
   }
   fscanf(tilefp,"}"); 
   if (fscanf(tilefp,"[%g]",&stoic_float)) 
      stoic[i]=stoic_float; else stoic[i]=1.0;
   fscanf(tilefp,"\n");
 }
 fscanf(tilefp,"}\n");

 fscanf(tilefp,"binding strengths=\n");
 strength = (double*) calloc(sizeof(double),num_bindings+1);
 fscanf(tilefp,"{");
 strength[0]=0;         /* bond type 0 ("null") always has strength 0 */
 for (i=1;i<=num_bindings;i++) {
   fscanf(tilefp,"%g",&strength_float);
   strength[i]=(double)strength_float;
 }
 fscanf(tilefp,"}\n");
 fclose(tilefp);
}

void getargs(int argc, char **argv)
{
 int i; 
 struct timeval tv; 
 gettimeofday(&tv, NULL); srand48(tv.tv_usec); srandom(tv.tv_usec);

 if (sizeof(long) != 4) {
   printf("sizeof long (%d) should be 4\n", sizeof(long int));
   exit(-1);
 }

 if (argc==2 && strcmp(argv[1],"--")==0) {
   printf("usage: xgrow tilefile [option=#]... \n");
   printf(" tilefile is an input file that specifies tiles\n");
   printf(" options:\n");
   printf("  block= display block size, 1...4\n");
   printf("  rand=  random number seed\n");
   printf("  k=     hybridization rate constant (/sec)\n");
   printf("  Gmc=   initiation free energy  (units kT)\n");
   printf("  Gse=   interaction free energy per binding\n");
   printf("  Gmch=  initiation free energy  for hydrolyzed units\n");
   printf("  Gseh=  interaction free energy for hydrolyzed units\n");
   printf("  Ghyd=  free energy of hydrolysis\n");
   printf("  Gas=   activation energy for spontaneous hydrolysis\n");
   printf("  Gam=   activation energy for mismatched sticky ends\n");
   printf("  Gae=   activation energy for unmatched sticky ends\n");
   printf("  Gah=   activation energy for hydrolyzed neighbors\n");
   printf("  Gao=   delta a. e. for output vs input-triggers hydrolysis\n");
   printf("  Gfc=   log concentration of flakes (otherwise no depletion)\n");
/* printf("  anneal=g/t        anneal Gse to g with time constant t\n"); */
   printf("  seed=i,j,n        seed tile type n at position i,j\n");
   printf("  stripe=o[:p,w]*   width w stripe with p errors, offset o\n");
   printf("  wander            wandering `seed' designation\n");
   printf("  periodic          periodic boundary conditions\n");
   printf("  -linear           simulate linear A B tiles, write errs > stdout \n");
   printf("  -nw               no X window (only if ?max set)\n");
   printf("  tmax=             quit after time t has passed\n");
   printf("  emax=             quit after e events have occurred\n");
   printf("  smax=             quit when the fragment is size s\n");
   printf("  datafile=         append Gmc, Gse, time, size, #mismatched se\n");
   printf("  arrayfile=        output matrix of final tiles\n");
   exit (0);
 }

 if ( (tilefp = fopen(&argv[1][0],"r"))!=NULL )
     { read_tilefile(tilefp); }
 else {
   printf("  First argument must be a tile file!\n");
   exit(0);   
 }
 

 tmax=0; emax=0; smax=0;
 wander=0; periodic=0; deplete=0; linear=0;
 Gfc=0; datafp=NULL; arrayfp=NULL; 
 Gmc=17; Gse=8.6; ratek = 1000000.0;
 Gmch=30; Gseh=0; Ghyd=30; Gas=30; Gam=15; Gae=30; Gah=30; Gao=10;
 seed_i=250; seed_j=250; seed_n=1; hydro=0;

 for (i=2; i<argc; i++) {
   if (strncmp(argv[i],"block=",6)==0) 
     { block=atoi(&argv[i][6]); if (block>4) block=4; if (block<1) block=1; }
   if (strncmp(argv[i],"rand=",5)==0) 
     { srand48(atoi(&argv[i][5])); srandom(atoi(&argv[i][5])); }
   if (strncmp(argv[i],"k=",2)==0) ratek=atof(&argv[i][2]);
   if (strncmp(argv[i],"Gmc=",4)==0) Gmc=atof(&argv[i][4]);
   if (strncmp(argv[i],"Gse=",4)==0) Gse=atof(&argv[i][4]);
   if (strncmp(argv[i],"Gmch=",5)==0) {hydro=1; Gmch=atof(&argv[i][5]);}
   if (strncmp(argv[i],"Gseh=",5)==0) {hydro=1; Gseh=atof(&argv[i][5]);}
   if (strncmp(argv[i],"Ghyd=",5)==0) {hydro=1; Ghyd=atof(&argv[i][5]);}
   if (strncmp(argv[i],"Gas=",4)==0) {hydro=1; Gas=atof(&argv[i][4]);}
   if (strncmp(argv[i],"Gam=",4)==0) {hydro=1; Gam=atof(&argv[i][4]);}
   if (strncmp(argv[i],"Gae=",4)==0) {hydro=1; Gae=atof(&argv[i][4]);}
   if (strncmp(argv[i],"Gah=",4)==0) {hydro=1; Gah=atof(&argv[i][4]);}
   if (strncmp(argv[i],"Gao=",4)==0) {hydro=1; Gao=atof(&argv[i][4]);}
   if (strncmp(argv[i],"Gfc=",4)==0) {deplete=1; Gfc=atof(&argv[i][4]);}
   if (strcmp(argv[i],"periodic")==0) periodic=1;
   if (strcmp(argv[i],"wander")==0) wander=1;
   if (strncmp(argv[i],"seed=",5)==0) {
      char *p=(&argv[i][5]);
      seed_i=atoi(p);
      if ((p=strchr(p,','))!=NULL) {
         seed_j=atoi(p+1);
         if ((p=strchr(p+1,','))!=NULL) seed_n=atoi(p+1);
      }
   }
   if (strncmp(argv[i],"stripe=",7)==0) 
      { stripe_args=(&argv[i][7]); periodic=1; wander=1; }
   if (strcmp(argv[i],"-nw")==0) XXX=0;
   if (strcmp(argv[i],"-linear")==0) linear=1;
   if (strncmp(argv[i],"tmax=",5)==0) tmax=atof(&argv[i][5]);
   if (strncmp(argv[i],"emax=",5)==0) emax=atoi(&argv[i][5]);
   if (strncmp(argv[i],"smax=",5)==0) smax=atoi(&argv[i][5]);
   if (strncmp(argv[i],"datafile=",9)==0) datafp=fopen(&argv[i][9], "a");
   if (strncmp(argv[i],"arrayfile=",10)==0) arrayfp=fopen(&argv[i][10], "w");
 }
 if (tmax==0 && emax==0 && smax==0) XXX=1;
}

void closeargs()
{ 
  int row,col,i;
  int size = (1<<fp->P); 


  clean_flake(fp);

  for (i=0;i<=fp->N;i++) free(units[i]); free(units);
  free(strength); free(stoic);

  if (datafp!=NULL) {
     if (fp->hydro) fprintf(datafp, " %f %f %f %f %f %f %f %f %f ",
       Gseh, Gmch, Ghyd, Gas, Gam, Gae, Gah, Gao, Gfc);
     fprintf(datafp, " %f %f %f %f %d %d %d\n",
       Gmc,Gse,ratek,fp->t,fp->stat_a-fp->stat_d,fp->mismatches,fp->events);
     fclose(datafp);
  } 
  if (arrayfp!=NULL) {
     fprintf(arrayfp,"\n");
     for (row=0; row<size; row++) {
        for (col=0; col<size; col++)
           fprintf(arrayfp, " %d", fp->Cell(row,col));
        fprintf(arrayfp, "\n");
     }
     fclose(arrayfp);
  }
  free_flake(fp);
}
  
  

#define errortile(i,j) ((fp->Cell(i,j)==0) ? 0: (                    \
         (units[fp->Cell(i,j)][1] != units[fp->Cell(i,(j)+1)][3] ||  \
          units[fp->Cell(i,j)][2] != units[fp->Cell((i)+1,j)][0]) ? 2:3)) 

#define getcolor(i,j) translate[ (err ?                               \
         (errortile(i,j)+(fp->hydro?m*(fp->Cell(i,j)>fp->N/2):0)) :   \
                        (fp->Cell(i,j))) ]

/* NOTE: requires 2^P < NCOLS+2*NBDY */
void showpic(flake *fp, int err) /* display the field */
{int row,col,i1,i2,color,j,j1,j2,blocktop=block, size, m;
 char *picture=(*spinimage).data;
 size = (1<<fp->P); 
 m= 2*(err>1);
 if (block>4) blocktop=block-1;
 if (8==(*spinimage).depth) {
  if (block>1) /* I wish I knew how to do this faster */
    for (row=0;row<size;row++)
      for (col=0;col<size;col++) {
        color = getcolor(row,col);
        j=block*((col+NBDY)+block*NCOLS*(row+NBDY));
        if (color!=picture[j])
	  for (i1=0;i1<blocktop;i1++) {
            j1=i1*block*NCOLS+j;
            for (i2=0;i2<blocktop;i2++)
               picture[j1+i2]=color;
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
       if (color!=XGetPixel(spinimage,j1=block*(col+NBDY),j2=block*(row+NBDY)))
         for (i2=0;i2<blocktop;i2++)
           for (i1=0;i1<blocktop;i1++)
            XPutPixel(spinimage,j1+i1,j2+i2,color);
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


/* fix up the pause button */
void setpause(int value)
{paused=value;
 if (paused) 
   XDrawImageString(display,pausebutton,gcr,0,font_height,"  run/PAUSE  ",13);
 else
   XDrawImageString(display,pausebutton,gcr,0,font_height,"  RUN/pause  ",13);
}

/* fix up the colors button */
void setcolor(int value)
{errorc=value;
 if (errorc==2) 
   XDrawImageString(display,colorbutton,gcr,0,font_height,"tile/err/HYD",12);
 else if (errorc==1) 
   XDrawImageString(display,colorbutton,gcr,0,font_height,"tile/ERR/hyd",12);
 else if (errorc==0)
   XDrawImageString(display,colorbutton,gcr,0,font_height,"TILE/err/hyd",12);
}

/* fix up the colors button */
void setwander(int value)
{wander=value;
 if (wander) 
   XDrawImageString(display,seedbutton,gcr,0,font_height,"fixed/WANDER",12);
 else
   XDrawImageString(display,seedbutton,gcr,0,font_height,"FIXED/wander",12);
}


/* this fixes the window up whenever it is uncovered */
void repaint()
{int i=0;
 XDrawString(display,quitbutton,  gcr,0,font_height,"    quit     ",13);
 XDrawString(display,fillbutton,  gcr,0,font_height,"    clear    ",13);
 XDrawString(display,randombutton,gcr,0,font_height,"random square",13);
 setpause(paused);
 setcolor(errorc);
 setwander(wander);

 /* write various strings */
 sprintf(stringbuffer,"%d by %d lattice: %d tiles, %d mismatches   ",
       (1<<fp->P),(1<<fp->P),fp->stat_a-fp->stat_d,fp->mismatches);
 XDrawImageString(display,window,gc,5,(++i)*font_height,
               stringbuffer,strlen(stringbuffer));
 if (wander) 
  sprintf(stringbuffer,"%s boundary.  seed i,j = %d,%d; n = %d   ",
        periodic?"periodic":"empty", fp->seed_i,fp->seed_j,fp->seed_n);
 else
  sprintf(stringbuffer, "%s boundary.", periodic?"periodic":"empty");
 XDrawImageString(display,window,gc,5,(++i)*font_height,
               stringbuffer,strlen(stringbuffer));

 sprintf(stringbuffer,"Gmc=%4.1f  Gse=%4.1f  k=%6.0f",
               Gmc,Gse,ratek);
 XDrawString(display,window,gc,5,(++i)*font_height,
               stringbuffer,strlen(stringbuffer));

 if (fp->hydro) {
   sprintf(stringbuffer,"Gmch=%4.1f  Gseh=%4.1f  Ghyd=%4.1f",
               Gmch,Gseh,Ghyd);
   XDrawString(display,window,gc,5,(++i)*font_height,
               stringbuffer,strlen(stringbuffer));
   sprintf(stringbuffer,"Gas=%4.1f Gam=%4.1f Gae=%4.1f Gah=%4.1f Gao=%4.1f",
               Gas,Gam,Gae,Gah,Gao);
   XDrawString(display,window,gc,5,(++i)*font_height,
               stringbuffer,strlen(stringbuffer));
 }

 sprintf(stringbuffer,"t = %12.3f sec; G = %12.3f      ",fp->t, fp->G);
 XDrawImageString(display,window,gc,5,(++i)*font_height,
               stringbuffer,strlen(stringbuffer));
 sprintf(stringbuffer, "%d events (%da,%dd,%dh,%df)     ",
        fp->events, fp->stat_a, fp->stat_d, fp->stat_h, fp->stat_f);
 XDrawImageString(display,window,gc,5,(++i)*font_height,
               stringbuffer,strlen(stringbuffer));

 if (deplete) {
  sprintf(stringbuffer, "Gfc=%4.1f; Gmc(%d)=%4.1f Gmc(%d)=%4.1f Gmc(%d)=%4.1f",
        Gfc, 1, (fp->conc[1]>0)?-log(fp->conc[1]):0,
             N/2, (fp->conc[N/2]>0)?-log(fp->conc[N/2]):0,
             N, (fp->conc[N]>0)?-log(fp->conc[N]):0);
  XDrawImageString(display,window,gc,5,(++i)*font_height,
               stringbuffer,strlen(stringbuffer));
 }

 XDrawString(display,window,gc,WINDOWWIDTH-120,WINDOWHEIGHT-5
       ,"EW '98-'02",10); 

 showpic(fp,errorc); 
}
 
/* a lot of this is taken from the basicwin program in the
   Xlib Programming Manual */
void openwindow(int argc, char **argv)
{char *window_name="xgrow";
 char *icon_name="xgrow";
 Pixmap icon_pixmap;
 char *display_name=NULL;
 XEvent report;
 XColor xcolor,colorcell;
 Colormap cmap;
 int i,j;
# define icon_bitmap_width 16
# define icon_bitmap_height 16
 static char icon_bitmap_bits[] = {
   0x1f, 0xf8, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0xf8,
   0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0xff, 0xff,
   0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};

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
 if (XAllocNamedColor(display,cmap,"firebrick",&colorcell,&xcolor)) {
   darkcolor=colorcell.pixel;
 }
 if (XAllocNamedColor(display,cmap,"wheat",&colorcell,&xcolor))
              lightcolor=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"black",&colorcell,&xcolor))
              translate[0]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"blue",&colorcell,&xcolor))
              translate[1]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"red",&colorcell,&xcolor))
              translate[2]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"green",&colorcell,&xcolor))
              translate[3]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"yellow",&colorcell,&xcolor))
              translate[4]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"gold",&colorcell,&xcolor))
              translate[5]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"purple",&colorcell,&xcolor))
              translate[6]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"white",&colorcell,&xcolor))
              translate[7]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"dark blue",&colorcell,&xcolor)) {
   translate[8]=colorcell.pixel;
 }
 if (XAllocNamedColor(display,cmap,"dark red",&colorcell,&xcolor))
              translate[9]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"dark green",&colorcell,&xcolor))
              translate[10]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"wheat",&colorcell,&xcolor))
              translate[11]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"orange",&colorcell,&xcolor))
              translate[12]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"cyan",&colorcell,&xcolor))
              translate[13]=colorcell.pixel;
 if (XAllocNamedColor(display,cmap,"light grey",&colorcell,&xcolor))
              translate[14]=colorcell.pixel;
  /* fill out the color table for future uses... only 14 colors!! */
 /* MUST MODIFY THIS FOR GENERAL TILE SETS? */
 for(i=15;i<256;i++)
   translate[i]=translate[(i-1)%14+1]; 

   /* make the main window */
 window=XCreateSimpleWindow(display,RootWindow(display,screen),
   0,0,WINDOWWIDTH,WINDOWHEIGHT,4,black,lightcolor);

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
    WINDOWWIDTH-140,WINDOWHEIGHT-50,120,20,2,black,darkcolor);
 pausebutton=XCreateSimpleWindow(display,window,
    WINDOWWIDTH-140,WINDOWHEIGHT-76,120,20,2,black,darkcolor);
 fillbutton=XCreateSimpleWindow(display,window,
    WINDOWWIDTH-140,WINDOWHEIGHT-102,120,20,2,black,darkcolor);
 colorbutton=XCreateSimpleWindow(display,window,
    WINDOWWIDTH-140,WINDOWHEIGHT-128,120,20,2,black,darkcolor);
 randombutton=XCreateSimpleWindow(display,window,
    WINDOWWIDTH-140,WINDOWHEIGHT-154,120,20,2,black,darkcolor);
 seedbutton=XCreateSimpleWindow(display,window,
    WINDOWWIDTH-140,WINDOWHEIGHT-180,120,20,2,black,darkcolor);
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
 XSelectInput(display,fillbutton,event_mask);
 XSelectInput(display,colorbutton,event_mask);
 XSelectInput(display,randombutton,event_mask);
 XSelectInput(display,seedbutton,event_mask);
 event_mask=ButtonReleaseMask|ButtonPressMask|PointerMotionHintMask
                   |ButtonMotionMask; 
 XSelectInput(display,playground,event_mask);

/* pick font: 9x15 is supposed to almost always be there */
 if ((font=XLoadQueryFont(display,"9x15"))==NULL)
  {fprintf(stderr,"%s: Cannot open 9x15 font\n",progname);
   exit(-1);
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
 XMapWindow(display,fillbutton);
 XMapWindow(display,colorbutton);
 XMapWindow(display,randombutton);
 XMapWindow(display,seedbutton);
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
 
  spinimage=XGetImage((Display *) display, (Drawable) playground,
            0,0,block*NCOLS,block*NROWS,
            AllPlanes,ZPixmap);
  if (NULL==spinimage)
        {fprintf(stderr,"trouble creating image structure\n");
         exit(-1);
        } 
  /* make sure everything get written first time */
  //  for (i=0;i<block*block*VOLUME;i++) spinimage->data[i]=translate[0];

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
{unsigned int width, height;
 int x,y,b,i,j;
 int stat=0;
 XEvent report;
 progname=argv[0];

 for (i=0;i<MAXTILETYPES;i++) {
   translate[i]=0;
 }

 getargs(argc, argv);

 if (hydro) { /* automatically double the number of tiles */
   units=(int**) realloc(units,sizeof(int*)*(2*N+1));
   for (i=1;i<=N;i++) {
     units[i+N]= (int*) calloc(sizeof(int),4);
     for (j=0;j<4;j++) units[i+N][j] = units[i][j];
   }
   stoic=(double*) realloc(stoic,sizeof(double)*(2*N+1));
   for (i=1;i<=N;i++) stoic[i+N]= stoic[i];
   N = N*2;
 }

 if (linear) {
    linear_simulate(ratek,Gmc,Gse,tmax,emax,smax);
    return 0;
 }

 if (XXX) openwindow(argc,argv);

 printf("xgrow: tile set read, beginning simulation\n");
 
 if (DEBUG==2) {
   /* set initial state: 2^3 grid */
   fp = init_flake(3,N,num_bindings);   
   set_params(fp,units,strength,stoic,hydro,ratek,Gmc,Gse,Gmch,Gseh,Ghyd,Gas,Gam,Gae,Gah,Gao,Gfc);   
   fp->seed_i=6; fp->seed_j=6;  
 } else {
   /* set initial state: 2^8 grid */
   fp = init_flake(8,N,num_bindings);   
   set_params(fp,units,strength,stoic,hydro,ratek,Gmc,Gse,Gmch,Gseh,Ghyd,Gas,Gam,Gae,Gah,Gao,Gfc);
   if (stripe_args==NULL) {
      fp->seed_i=seed_i; fp->seed_j=seed_j; fp->seed_n=seed_n; 
   } else {
      int i,j,k,w; double p; int size=(1<<fp->P); char *s=stripe_args;
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
                  fp->events--; /* don't count these as events */
                  /* ERROR: stats are also modified!!! */
	       }
	    }
	 }
      }
      /* no corner or boundary tiles for stripe simulations */
      fp->conc[0] -= fp->conc[1]; fp->conc[0]+=(fp->conc[1]=exp(-35));
      fp->conc[0] -= fp->conc[2]; fp->conc[0]+=(fp->conc[2]=exp(-35));
      fp->conc[0] -= fp->conc[3]; fp->conc[0]+=(fp->conc[3]=exp(-35));
   }
 }

 /* loop forever, looking for events */
 while((tmax==0 || fp->t < tmax) && 
       (emax==0 || fp->events < emax) &&
       (smax==0 || fp->stat_a-fp->stat_d < smax)) { 
  if (!XXX) 
     simulate(fp,NUPDATES,tmax,emax,smax);
  else {
   if (0==paused && !XPending(display)) {
     simulate(fp,NUPDATES,tmax,emax,smax);
     stat++; if (stat==NSTATS) { stat=0; repaint(); }
   }
   if (paused|XPending(display))
    {XNextEvent(display,&report); 
     switch (report.type)
      {case Expose:
        if (report.xexpose.count!=0) break; /* more in queue, wait for them */
        repaint();  
        break;
       case ConfigureNotify:
        width=report.xconfigure.width;
        height=report.xconfigure.height;
        if ((width<size_hints.min_width)||(height<size_hints.min_height))
            {fprintf(stderr,"%s: window too small to proceed.\n",progname);
             cleanup();
            } 
        break; 
       case ButtonPress:
        if (report.xbutton.window==quitbutton)
            { cleanup(); } 
        else if (report.xbutton.window==pausebutton) 
            { setpause(1-paused); repaint(); }
        else if (report.xbutton.window==fillbutton)
            {free_flake(fp); fp=init_flake(8,N,num_bindings); 
             set_params(fp,units,strength,stoic,hydro,ratek,Gmc,Gse,Gmch,Gseh,Ghyd,Gas,Gam,Gae,Gah,Gao,Gfc);
             fp->seed_i=seed_i; fp->seed_j=seed_j; fp->seed_n=seed_n; 
             repaint();
            }
        else if (report.xbutton.window==colorbutton)
            { setcolor((errorc+1)%3); repaint(); }
        else if (report.xbutton.window==seedbutton)
            { setwander(1-wander); repaint(); }
        else if (report.xbutton.window==randombutton) {
             for (x=NCOLS/2-10;x<NCOLS/2+10;x++)
              for (y=NROWS/2-10;y<NROWS/2+10;y++)
/*               change_cell(fp,x,y,lrand48()%(fp->N+1)); */
               change_cell(fp,x,y,4+((x+y)%2));  /* bar code */
             repaint();
        } else if (report.xbutton.window==playground) 
            {x=report.xbutton.x/block;
             y=report.xbutton.y/block;
             b=report.xbutton.button;
             /* sketch(x,y,b); */
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


