
# include <stdio.h>
# include <stdlib.h>
# include <malloc.h>
# include <math.h>
# include <assert.h>

# include "grow.h"

/* index is 8 bits N NE E SE S SW W NW where                         */
/*  N  E  S  W    refer to where there is a tile present, and        */
/*   NE SE SW NW  indicate whether the corner tile connects its      */
/*                neighbors (given that all three are present)       */
/* then ring[index]=1 if there is 1 (or fewer) connected groups.     */
unsigned char ring[256];
#define ROTATE(i)          ((((i)&1)<<7) + ((i)>>1))
#define ROTATE_CLEAR(i)    (               ((i)>>1))

/* sets up data structures -- cell field, hierarchical rates, etc    */
flake *init_flake(unsigned char P, unsigned char N, int num_bindings)
{
  int i,j,n,m,p;
  int size = (1<<P);
  flake *fp = (flake *)malloc(sizeof(flake));

  fp->P = P; fp->N = N; fp->num_bindings = num_bindings;

  fp->units = (int**) calloc(sizeof(int*),N+1);
  for (i=0;i<=N;i++) {
   fp->units[i] = (int*) calloc(sizeof(int),4);
   for (j=0;j<4;j++) fp->units[i][j]=0;
  }
  fp->strength = (double*) calloc(sizeof(double),num_bindings+1);

  fp->cell = (unsigned char **)calloc(sizeof(char *),2+size);
  for (i=0;i<2+size;i++) 
     fp->cell[i]=(unsigned char *)calloc(sizeof(char),2+size);
  fp->rate = (double ***)calloc(sizeof(char **),P+1);
  fp->empty = (int ***)calloc(sizeof(int **),P+1);
  for (p=0;p<=P;p++) {
     size = (1<<p);
     fp->rate[p]=(double **)calloc(sizeof(double *),size);
     fp->empty[p]=(int **)calloc(sizeof(int *),size);
     for (i=0;i<size;i++) {
       fp->rate[p][i]=(double *)calloc(sizeof(double),size);
       fp->empty[p][i]=(int *)calloc(sizeof(int),size);
       for (j=0;j<size;j++) fp->rate[p][i][j]=0;
     }
  }

  fp->flake_conc=0;
  fp->conc = (double *)calloc(sizeof(double),N+1);
  for (n=0;n<N+1;n++) fp->conc[n]=0;
  fp->Gcb  = (double *)calloc(sizeof(double),N+1);
  for (n=0;n<N+1;n++) fp->Gcb[n]=0;
  fp->Gse_EW = (double **)calloc(sizeof(double *),N+1);
  for (n=0;n<N+1;n++) fp->Gse_EW[n]=(double *)calloc(sizeof(double),N+1);
  for (n=0;n<N+1;n++) for (m=0;m<N+1;m++) fp->Gse_EW[n][m]=0;
  fp->Gse_NS = (double **)calloc(sizeof(double *),N+1);
  for (n=0;n<N+1;n++) fp->Gse_NS[n]=(double *)calloc(sizeof(double),N+1);
  for (n=0;n<N+1;n++) for (m=0;m<N+1;m++) fp->Gse_NS[n][m]=0;

  fp->G=0; fp->events=0; fp->t=0;
  fp->seed_i=0; fp->seed_j=0; fp->seed_n=1;
  fp->stat_a=fp->stat_d=fp->stat_h=fp->stat_f=0;
  fp->mismatches=0;

  fp->rv  = (double *)calloc(sizeof(double),N+1);
  fp->Fnext = (int *)calloc(sizeof(int),size*size);
  for (n=0;n<size*size;n++) fp->Fnext[n]=-1;
  fp->Fgroup = (int *)calloc(sizeof(int),size*size);

  /* set_params() will have to put reasonable values in place */
  /* note that empty and rate are correct, because there are no tiles yet */
  /* NOTE this routine is not "proper" -- assumes all allocs are OK */

  return fp;
}

void free_flake(flake *fp)
{
  int i,n,p;
  int size = (1<<fp->P);
  
  for (i=0;i<2+size;i++) free(fp->cell[i]); free(fp->cell);

  for (p=0;p<=fp->P;p++) {
     size = (1<<p);
     for (i=0;i<size;i++) {
       free(fp->rate[p][i]);
       free(fp->empty[p][i]);
     }
     free(fp->rate[p]);
     free(fp->empty[p]);
  }
  free(fp->rate);
  free(fp->empty);

  free(fp->conc);
  free(fp->Gcb);
  for (n=0;n<fp->N+1;n++) free(fp->Gse_EW[n]); free(fp->Gse_EW);
  for (n=0;n<fp->N+1;n++) free(fp->Gse_NS[n]); free(fp->Gse_NS);

  free(fp->rv); free(fp->Fgroup); free(fp->Fnext);
  
  for (n=0;n<fp->N+1;n++) free(fp->units[n]); free(fp->units);
  free(fp->strength); 

  /* NOTE because init_flake is not safe to out-of-mem, this could die */
}

/* set up info for Sierpinski rules.  */
/* later: add hydrolysis rules */
/* later still: read this from file */
void set_params(flake *fp, int** units, double* strength, double* stoic,
 int hydro, double k, double Gmc, double Gse,
 double Gmch, double Gseh, double Ghyd, 
 double Gas, double Gam, double Gae, double Gah, double Gao,
 double Gfc)
{
   int i,j,n,m;

   /* make our own copy of tile set and strengths, so the calling program
      can do as it pleases with it's units & strength information */
   for (i=0;i<=fp->N;i++) for (j=0;j<4;j++) fp->units[i][j] = units[i][j];
   for (i=0;i<=fp->num_bindings;i++) fp->strength[i] = strength[i]; 

   fp->hydro = hydro;

   fp->k = k;
   fp->kas = exp(-Gas);  fp->kao = exp(-Gao);
   fp->kam = exp(-Gam);  fp->kae = exp(-Gae); fp->kah = exp(-Gah);
   
   /* set fp->conc from Gmc... and from Gmch for hydrolysis rules */
   fp->conc[0]=0;  
   for (n=1; n <= fp->N; n++) 
      fp->conc[0]+=
        (fp->conc[n]=exp(-((n>fp->N/2&&fp->hydro)?Gmch:Gmc))*stoic[n]);

   fp->G = -log(fp->conc[1]); /* for corner seed tile by itself, total G=0 */
   if (Gfc>0) fp->flake_conc = exp(-Gfc);

   /* set fp->Gcb from Ghyd */
   for (n=0; n <= fp->N; n++) fp->Gcb[n]=0;
   if (fp->hydro) for (n=fp->N/2+1; n <= fp->N; n++) fp->Gcb[n]=Ghyd;

   /* set Gse_EW Gse_NS from Gse, Gseh rules */
   /* borrowed from old calcG() routine; uses (fp->units)[] and fp->strength[] */
   for (n=1; n<=fp->N; n++)
      for (m=1; m<=fp->N; m++) {
         fp->Gse_EW[n][m] = ((fp->units)[n][3]==(fp->units)[m][1]) *
              (fp->strength)[(fp->units)[m][1]] * (fp->hydro?(((n > fp->N/2) || (m > fp->N/2))?Gseh:Gse):Gse);
         fp->Gse_NS[n][m] = ((fp->units)[n][2]==(fp->units)[m][0]) *
              (fp->strength)[(fp->units)[m][0]] * (fp->hydro?(((n>fp->N/2 || m>fp->N/2))?Gseh:Gse):Gse);
      }

   /* make sure ring[] has entries */
   ring[0]=1; ring[255]=1;
   for (i=1; i<255; i++) {
      n=i;
      while ((n&1)==1) n=ROTATE(n);
      while ((n&1)==0) n=ROTATE(n);  
      /* now we have the beginning of a block of 1s */
      while ((n&1)==1) n=ROTATE_CLEAR(n);
      /* if there was only one block of 1s, 
         then they've all been erased now */
      ring[i] = (n==0);
   }
}


/* macro definition of summed sticky end bond energy                    */
/* computes energy IF Cell(i,j) were n, given its current neighbors     */
/* assumes "fp" arg is a simple variable, but others can be expressions */
/* note that n != 0  and assumes 0 <= i,j < (1<<fp->P)                  */
#define Gse(fp,i,j,n) (                          \
        fp->Gse_EW[ n ] [ fp->Cell(i,(j)-1) ] +  \
        fp->Gse_EW[ fp->Cell(i,(j)+1) ] [ n ] +  \
        fp->Gse_NS[ n ] [ fp->Cell((i)+1,j) ] +  \
        fp->Gse_NS[ fp->Cell((i)-1,j) ] [ n ] )

/* similar definition to count the number of ses that are mismatched    */
/* also,  n != 0   and assumes 0 <= i,j < (1<<fp->P)                    */

#define Mism(fp,i,j,n) (                           \
 ((fp->units)[n][1] != (fp->units)[fp->Cell(i,(j)+1)][3] &&    \
  (fp->units)[n][1]*(fp->units)[fp->Cell(i,(j)+1)][3] > 0) +   \
 ((fp->units)[n][3] != (fp->units)[fp->Cell(i,(j)-1)][1] &&    \
  (fp->units)[n][3]*(fp->units)[fp->Cell(i,(j)-1)][1] > 0) +   \
 ((fp->units)[n][2] != (fp->units)[fp->Cell((i)+1,j)][0] &&    \
  (fp->units)[n][2]*(fp->units)[fp->Cell((i)+1,j)][0] > 0) +   \
 ((fp->units)[n][0] != (fp->units)[fp->Cell((i)-1,j)][2] &&    \
  (fp->units)[n][0]*(fp->units)[fp->Cell((i)-1,j)][2] > 0) )

/* gives concentration-independent rates                                   */
/* for non-empty cells i,j, computes rate rv[n] to convert to type n       */
/*  so rv[0] gives the off-rate  (not the sum)                             */
/*  n=0 (empty) is the most common case.  the returned value is sum rv[n]  */
/* for empty cells i,j, returns 0                                          */
/* rv must already exist, of size fp->N+1                                  */
/* 0 <= i,j < 2^P                                                          */
double calc_rates(flake *fp, int i, int j, double *rv)
{
  int n,mi,ei,hi,mo,eo,ho; double r, sumr;
  unsigned char N,E,S,W;

  if (rv!=NULL) for (n=0;n<=fp->N;n++) rv[n]=0;
  n = fp->Cell(i,j);
  if (n==0) return 0;                           /* no off-rate for empties   */
  if (i==fp->seed_i && j==fp->seed_j) return 0; /* seed site doesn't go away */

  r = fp->k * exp(-Gse(fp,i,j,n)); 
  if (rv!=NULL) rv[0]=r; sumr = r;

  /*  hydrolysis model assumes tiles 1...N/2 non hydrolyzed, N/2+1...N hydro */
  if (n > 0 && fp->hydro) {
    /*
       for "hydrolysis" of tile n -> n+N/2, for n = 1...N/2 
       and for "dehydrolysis"   n -> n-N/2, for n = N/2+1...N
       we need to know Gse, Gseh, Ghyd, kao = exp(-Gao),
            kas = exp(-Gas), kae = exp(-Gae), kam = exp(-Gam), kah = exp(-Gah)
         and to calc mi (# input bindings mismatched)
                     ei (# input bindings empty)
         NOTE mismatch => both tiles have se types 1..3
              empty    => tile input se have types 1..3; output has no se
    */
    /* empty and mismatched inputs */
    S=fp->Cell(i+1,j); E=fp->Cell(i,j+1);
    ei = ((fp->units)[n][2] != 0 && (fp->units)[S][0]==0) + 
         ((fp->units)[n][1] != 0 && (fp->units)[E][3]==0);
    mi = ((fp->units)[n][2] != (fp->units)[S][0] && (fp->units)[n][2]*(fp->units)[S][0] !=0) +
         ((fp->units)[n][1] != (fp->units)[E][3] && (fp->units)[n][1]*(fp->units)[E][3] !=0);
    hi = (E>fp->N/2) + (S>fp->N/2);

    /* empty and mismatched outputs */
    N=fp->Cell(i-1,j); W=fp->Cell(i,j-1);
    eo = ((fp->units)[n][0] != 0 && (fp->units)[N][2]==0) + 
         ((fp->units)[n][3] != 0 && (fp->units)[W][1]==0);
    mo = ((fp->units)[n][0] != (fp->units)[N][2] && (fp->units)[n][0]*(fp->units)[N][2] !=0) +
         ((fp->units)[n][3] != (fp->units)[W][1] && (fp->units)[n][3]*(fp->units)[W][1] !=0);
    ho = (N>fp->N/2) + (W>fp->N/2);

    r = fp->k * (fp->kas + mi * fp->kam + ei * fp->kae + hi * fp->kah +
                fp->kao * (mo * fp->kam + eo * fp->kae + ho * fp->kah));
    if (n<=fp->N/2) {
      if (rv!=NULL) rv[n+fp->N/2]=r; sumr += r;
    } else {
      r = r * exp(-fp->Gcb[n]+fp->Gcb[n-fp->N/2]-Gse(fp,i,j,n)+Gse(fp,i,j,n-fp->N/2));
      if (rv!=NULL) rv[n-fp->N/2]=r; sumr += r;
    }
  } 

  return sumr;
}


/* figure the delta_rate for this cell, and propagate up            */
/* Cell(i,j) has just changed to n, but old rates and empty status  */
/*   has not yet changed to reflect this.  Since the cell and its   */
/*   neighbors must all be changed, we only guarantee that THIS     */
/*   cell's hierarchical status has not been updated yet.           */
/*   ii,jj is a neighbor cell; possibly a border cell.              */
/*   Note: "empty" means "empty and neighbor to a tile".            */
void update_rates(flake *fp, int ii, int jj)
{
   int p; int size=(1<<fp->P);
   double oldrate, newrate, oldempty, newempty;

   if (periodic) { ii=(ii+size)%size; jj=(jj+size)%size; }

   if (!(ii < 0 || ii >= size || jj < 0 || jj >= size)) {
      oldrate  = fp->rate[fp->P][ii][jj];
      oldempty = fp->empty[fp->P][ii][jj];
      newempty = (fp->Cell(ii,jj)==0) &&
                 ( fp->Cell(ii+1,jj) || fp->Cell(ii,jj+1) ||
                   fp->Cell(ii-1,jj) || fp->Cell(ii,jj-1) );
      newrate  = calc_rates(fp, ii, jj, NULL);
      for (p=fp->P; p>=0; p--) {
         fp->rate[p][ii][jj] += newrate-oldrate;
         fp->empty[p][ii][jj] += newempty-oldempty;
         ii = (ii>>1); jj = (jj>>1);
      }
   }
}


/* convert Cell(i,j) to type n.  */
/* ALL modifications to the cell array go through this interface!   */
/* we only modify cells with 0 <= i,j < (1<<fp->P)                  */
/* but since it might be called out of range in FILL, we fix it up. */
void change_cell(flake *fp, int i, int j, unsigned char n)
{
   int size=(1<<fp->P);
   if (periodic) { i=(i+size)%size; j=(j+size)%size; }

   if (fp->Cell(i,j)==n) return;
   if (fp->Cell(i,j)==0) {                         /* tile addition */
     fp->G += -log(fp->conc[n]) - Gse(fp,i,j,n);   
     fp->conc[0] -= (fp->conc[n]>fp->flake_conc)?fp->flake_conc:fp->conc[n]; 
     fp->conc[n] -= (fp->conc[n]>fp->flake_conc)?fp->flake_conc:fp->conc[n]; 
     fp->stat_a++;
     fp->mismatches += Mism(fp,i,j,n);
   } else if (n==0) {                              /* tile loss */
     fp->conc[0] += fp->flake_conc; 
     fp->conc[fp->Cell(i,j)] += fp->flake_conc; 
     fp->G += log(fp->conc[fp->Cell(i,j)]) + Gse(fp,i,j,fp->Cell(i,j));
     fp->stat_d++;
     fp->mismatches -= Mism(fp,i,j,fp->Cell(i,j));
   } else {                                        /* tile hydrolysis */
     fp->G += Gse(fp,i,j,fp->Cell(i,j)) - Gse(fp,i,j,n) +
              fp->Gcb[fp->Cell(i,j)] - fp->Gcb[n]; 
     fp->stat_h++;
     /* by our rules, hydrolyzed tiles have same se types as non-hyd. */
   }
   fp->events++;

   fp->Cell(i,j)=n; 
   if (periodic) { int size=(1<<fp->P);
     if (i==0)      fp->Cell(size,j)=n;
     if (i==size-1) fp->Cell(-1,j)=n;
     if (j==0)      fp->Cell(i,size)=n;
     if (j==size-1) fp->Cell(i,-1)=n;
   }

   update_rates(fp,i,j);
   update_rates(fp,i+1,j);
   update_rates(fp,i-1,j);
   update_rates(fp,i,j+1);
   update_rates(fp,i,j-1);
}

/* use rates & empty & conc to choose a cell to change,      */
/* and call calc_rates to identify what change to make.      */
/* report choice, but don't act on it.                       */
double choose_cell(flake *fp, int *ip, int *jp, int *np)
{
  double dt,sum,cum,r,kc,k00,k01,k10,k11;
  int p,i,j,di,dj,n,oops;

  kc = fp->k*fp->conc[0];
  sum = fp->rate[0][0][0] + kc*fp->empty[0][0][0];
  dt = -log(drand48())/sum;

  i=0; j=0;
  for (p=0; p<fp->P; p++) { /* choosing subquadrant from within p:i,j */
     sum = fp->rate[p][i][j] + kc*fp->empty[p][i][j]; 
     k00 = fp->rate[p+1][2*i][2*j]+kc*fp->empty[p+1][2*i][2*j];
     k10 = fp->rate[p+1][2*i+1][2*j]+kc*fp->empty[p+1][2*i+1][2*j];
     k01 = fp->rate[p+1][2*i][2*j+1]+kc*fp->empty[p+1][2*i][2*j+1];
     k11 = fp->rate[p+1][2*i+1][2*j+1]+kc*fp->empty[p+1][2*i+1][2*j+1];
     /* avoid possible round-off error... but still check for it */
     /* sum =?= k00+k01+k10+k11; */
     r = drand48() * (k00+k01+k10+k11);  cum = 0; oops=0;
     d2printf("%f / %f for choosing %d from %d: %d %d\n",r,sum,p+1,p,i,j);
     if (r < (cum += k00)) { di=0; dj=0; } else
     if (r < (cum += k10)) { di=1; dj=0; } else
     if (r < (cum += k01)) { di=0; dj=1; } else
     if (r < (cum += k11)) { di=1; dj=1; } else { di=1; dj=1; oops=1; }
     if (oops)
     {
         /* lots of round-off error has hosed this sum */
         dprintf("Four rose petals (%f) make no rose (%f) [at %d: %d,%d]!!!\n",
             cum,sum,p,i,j);
     }
     /* always fix-up any numerical error that could have accumulated here */
     fp->rate[p][i][j] = fp->rate[p+1][2*i][2*j]+fp->rate[p+1][2*i][2*j+1]+
                         fp->rate[p+1][2*i+1][2*j]+fp->rate[p+1][2*i+1][2*j+1];
     i=2*i+di; j=2*j+dj;
  }
  *ip=i; *jp=j;

  if (fp->Cell(i,j) == 0) {   /* choose on-event for type 1...N            */
     /* depletion seems to get conc[0] out of wack --- no, that was a bug
     if (fp->flake_conc>0) 
       for(fp->conc[0]=0,n=1;n<=fp->N;n++) fp->conc[0]+=fp->conc[n]; */
     r = drand48() * fp->conc[0];  cum = 0;
     for (n=1; n<=fp->N; n++) if (r < (cum += fp->conc[n])) break; 
     if (n>fp->N) 
       { printf("Eye of the needle!!! %f =!= %f\n",fp->conc[0],cum); n=0; }
  } else {                    /* choose off-event 0 or conversion to 1...N */
     sum = calc_rates(fp,i,j,fp->rv);
     if (sum==0) printf("Zero-sum game!!!\n");
     r = drand48() * sum;  cum = 0;
     for (n=0; n<=fp->N; n++) if (r < (cum += fp->rv[n])) break; 
     if (n>fp->N) { printf("A rose is not a rose!!!\n"); n=0; }
  }
  *np = n;
  
  return dt;
}

/*************** the fill routine for checking connectedness ***********/

/* A tile is defined to be HCONNECTED to a neighbor if they both have   */
/* non-zero se types  (who cares what the Gse is).                      */
/* This is hypothetical on i,j being tile n != 0.                       */
/* CONNECTED is the non-hypothetical version.                           */
#define HCONNECTED_N(fp,i,j,n) \
    ((fp->units)[n][0]!=0 && (fp->units)[fp->Cell((i)-1,j)][2]!=0)
#define HCONNECTED_E(fp,i,j,n) \
    ((fp->units)[n][1]!=0 && (fp->units)[fp->Cell(i,(j)+1)][3]!=0)
#define HCONNECTED_S(fp,i,j,n) \
    ((fp->units)[n][2]!=0 && (fp->units)[fp->Cell((i)+1,j)][0]!=0)
#define HCONNECTED_W(fp,i,j,n) \
    ((fp->units)[n][3]!=0 && (fp->units)[fp->Cell(i,(j)-1)][1]!=0)
#define HCONNECTED(fp,i,j,n) \
    ( HCONNECTED_N(fp,i,j,n) || HCONNECTED_E(fp,i,j,n) || \
      HCONNECTED_S(fp,i,j,n) || HCONNECTED_W(fp,i,j,n) )
#define CONNECTED_N(fp,i,j) HCONNECTED_N(fp,i,j,fp->Cell(i,j))
#define CONNECTED_E(fp,i,j) HCONNECTED_E(fp,i,j,fp->Cell(i,j))
#define CONNECTED_S(fp,i,j) HCONNECTED_S(fp,i,j,fp->Cell(i,j))
#define CONNECTED_W(fp,i,j) HCONNECTED_W(fp,i,j,fp->Cell(i,j))
#define CONNECTED(fp,i,j) HCONNECTED(fp,i,j,fp->Cell(i,j))
/* all of these can only be used for 0 <= i,j < size                     */

/* Cell i,j has just dissociated.  Previously, every cell was connected  */
/* to the seed cell. To maintain this property, we do a fill from each   */
/* neighbor, and we see if they connect to each other or to the seed.    */
/* A queue is used, so that the fill stays as local as possible.         */
/* The queue is in the array fp->Fnext[n], storing linked lists such     */
/* that no cell is in more than one linked list (careful w/ Fpush!), and */
/* the fill memory is in fp->Fgroup[n] identifying the group class.      */
/* Fnext & Fgroup are assumed to contain all "-1" & "0" before & after.  */
/* NOTE: Qn() packs (i,j) into n; so boundary positions can't be packed. */
/* Thus Fpush() implements periodic boundary conditions, which are only  */
/* used if the boundary is indeed periodic (ie non-zero).                */
#define Qi(n) ((n)/size)
#define Qj(n) ((n)%size)
#define Qn(i,j) ((((i)+size)%size)*size+(((j)+size)%size))
#define Fempty(g) (head[g]==-1)
#define Fpush(g,i,j) {                                      \
  if (Fempty(g)) head[g]=tail[g]=Qn(i,j);                   \
  else { fp->Fnext[tail[g]]=Qn(i,j); tail[g]=Qn(i,j); }     \
  if ((((i)+size)%size)==fp->seed_i &&                      \
      (((j)+size)%size)==fp->seed_j) seeded[ming[g]]=1;     \
  fp->Fgroup[Qn(i,j)]=g;}
#define Fpull(g,i,j) { int oldh=head[g];                    \
  i=Qi(head[g]); j=Qj(head[g]);                             \
  head[g]=fp->Fnext[head[g]]; fp->Fnext[oldh]=-1;           \
  if (Fempty(g)) tail[g]=-1; }
#define Fmerge(g1,g2) {                                                   \
  if (ming[g1]!=ming[g2]) {                                               \
     int gi,newg = MIN(ming[g1],ming[g2]),oldg = MAX(ming[g1],ming[g2]);  \
     for (gi=1;gi<5;gi++) if (ming[gi]==oldg) ming[gi]=newg;              \
     seeded[newg] = (seeded[newg]||seeded[oldg]); ngroups--;              \
  } }


void flake_fission(flake *fp, int ii, int jj)
{ 
                 /* groups are 0=clear, 1=E, 2=S, 3=W, 4=N.           */
  int head[5],   /* first pointer.  PULL grabs this guy to process    */
      tail[5],   /* NULL pointer. PUSH adds new guys here.            */
      ming[5],   /* the min group # of groups that are connected.     */
                 /* always = to the smallest # of all merged groups.  */
                 /* only groups in range(ming) have active Qs.        */
      seeded[5]; /* is the group connected to the seed?               */
                 /* only required to be valid for range(ming)         */
  int ngroups=0; /* number of connected groups (w/Q or w/o)           */
  int size = (1<<fp->P);
  int i,j,g,gg,implicit,active;
 
  for (g=0;g<5;g++) { head[g]=tail[g]=-1; ming[g]=-1; seeded[g]=0; }
  ming[0]=0;

  if (DEBUG==2) { int n; 
     for (i=j=n=0;n<size*size;n++) 
         { i+=(fp->Fgroup[n] != 0); j+=(fp->Fnext[n] != -1); }
     if (i>0 || j>0) printf("%d groupies and %d wannabies on entry\n",i,j);
  }

  /* start filling the neighbors until they all connect, or until       */
  /* the only unfinished group is known implicitly to be the seed group */
  /* (they could also all just finish, but that's exceedingly unlikely) */
  i=ii; j=jj;
  if (fp->Cell(i,j+1)!=0) { ngroups++; ming[1]=1; Fpush(1,i,j+1) }
  if (fp->Cell(i+1,j)!=0) { ngroups++; ming[2]=2; Fpush(2,i+1,j) }
  if (fp->Cell(i,j-1)!=0) { ngroups++; ming[3]=3; Fpush(3,i,j-1) }
  if (fp->Cell(i-1,j)!=0) { ngroups++; ming[4]=4; Fpush(4,i-1,j) }
  do {
    for (g=1; g<5; g++) {
      if (!Fempty(g)) {
        Fpull(g,i,j)
        if (CONNECTED_E(fp,i,j)) {
          if ((gg=fp->Fgroup[Qn(i,j+1)])==0) { Fpush(ming[g],i,j+1) }
          else { Fmerge(g,gg) }
        }
        if (CONNECTED_S(fp,i,j)) {
          if ((gg=fp->Fgroup[Qn(i+1,j)])==0) { Fpush(ming[g],i+1,j) }
          else { Fmerge(g,gg) }
        }
        if (CONNECTED_W(fp,i,j)) {
          if ((gg=fp->Fgroup[Qn(i,j-1)])==0) { Fpush(ming[g],i,j-1) }
          else { Fmerge(g,gg) }
        }
        if (CONNECTED_N(fp,i,j)) {
          if ((gg=fp->Fgroup[Qn(i-1,j)])==0) { Fpush(ming[g],i-1,j) }
          else { Fmerge(g,gg) }
        }
      }
    }
    active = (!Fempty(1))+(!Fempty(2))+(!Fempty(3))+(!Fempty(4));
    implicit = ( seeded[1]+seeded[2]+seeded[3]+seeded[4]==0 && active==1 );
  } while (active>0 && !implicit && ngroups>1);

  if (ngroups==1) { /* all groups merged to one. clean up; we're fine */
    while (!Fempty(1)) { Fpull(1,i,j) }
    while (!Fempty(2)) { Fpull(2,i,j) }
    while (!Fempty(3)) { Fpull(3,i,j) } 
    while (!Fempty(4)) { Fpull(4,i,j) }
    /* fp->Fnext is now all -1, as are head & tail */
  } else {
    /* more than one group.  who has the seed? */
    /* either all Qs are now empty, or */
    /* one group may not have finished. it has seed if others don't */
    if (implicit) {
       for (g=1; Fempty(g); g++); /* now g is the non-empty group */ 
       seeded[g]=1; while (!Fempty(g)) { Fpull(g,i,j) }
    }  /* fp->Fnext is now all -1, as are head & tail */

    /* now re-zero non-seeeded Fgroups, and dissociate tiles*/
    i=ii; j=jj;
    if (fp->Cell(i,j+1)!=0 && seeded[ming[1]]==0) 
          { change_cell(fp,i,j+1,0); Fpush(0,i,j+1) }
    if (fp->Cell(i+1,j)!=0 && seeded[ming[2]]==0) 
          { change_cell(fp,i+1,j,0); Fpush(0,i+1,j) }
    if (fp->Cell(i,j-1)!=0 && seeded[ming[3]]==0) 
          { change_cell(fp,i,j-1,0); Fpush(0,i,j-1) }
    if (fp->Cell(i-1,j)!=0 && seeded[ming[4]]==0) 
          { change_cell(fp,i-1,j,0); Fpush(0,i-1,j) }
    while (!Fempty(0)) {
       Fpull(0,i,j)
       if (fp->Fgroup[Qn(i,j+1)]>0 && seeded[ming[fp->Fgroup[Qn(i,j+1)]]]==0) 
          { change_cell(fp,i,j+1,0); Fpush(0,i,j+1) }
       if (fp->Fgroup[Qn(i+1,j)]>0 && seeded[ming[fp->Fgroup[Qn(i+1,j)]]]==0)
          { change_cell(fp,i+1,j,0); Fpush(0,i+1,j) }
       if (fp->Fgroup[Qn(i,j-1)]>0 && seeded[ming[fp->Fgroup[Qn(i,j-1)]]]==0) 
          { change_cell(fp,i,j-1,0); Fpush(0,i,j-1) }
       if (fp->Fgroup[Qn(i-1,j)]>0 && seeded[ming[fp->Fgroup[Qn(i-1,j)]]]==0)
          { change_cell(fp,i-1,j,0); Fpush(0,i-1,j) }
    }
    /* again fp->Fnext is now all -1, as are head & tail */
  }    

  /* now re-zero Fgroup for stuff that sticks around */
  i=ii; j=jj;
  if (fp->Cell(i,j+1)!=0) { Fpush(0,i,j+1) }
  if (fp->Cell(i+1,j)!=0) { Fpush(0,i+1,j) }
  if (fp->Cell(i,j-1)!=0) { Fpush(0,i,j-1) }
  if (fp->Cell(i-1,j)!=0) { Fpush(0,i-1,j) }
  while (!Fempty(0)) {
     Fpull(0,i,j)
     if (fp->Fgroup[Qn(i,j+1)]>0) { Fpush(0,i,j+1) }
     if (fp->Fgroup[Qn(i+1,j)]>0) { Fpush(0,i+1,j) }
     if (fp->Fgroup[Qn(i,j-1)]>0) { Fpush(0,i,j-1) }
     if (fp->Fgroup[Qn(i-1,j)]>0) { Fpush(0,i-1,j) }
  }
  /* now fp->Fnext is all -1 and fp->Fgroup is all 0 */
  fp->stat_f += ngroups-1;
}



/* remove all tiles whose off-rate is faster than the on-rate */
void clean_flake(flake *fp)
{
  int i,j; double kc;
  int size = (1<<fp->P);

  kc = fp->k*fp->conc[0]; /* on-rate */

  /* first memorize, then remove, to avoid changing rates during removal */
  for (i=0; i<size; i++)
    for (j=0; j<size; j++)
      fp->Fgroup[i+size*j] = (fp->rate[fp->P][i][j] > kc);
  for (i=0; i<size; i++)
    for (j=0; j<size; j++)
      if (fp->Fgroup[i+size*j]) change_cell(fp, i, j, 0);
  
}

/* simulates 'events' events */
void simulate(flake *fp, int events, double tmax, int emax, int smax)
{
   int i,j,n,p,oldn; double dt;
   unsigned char ringi; 

   emax = (emax==0 || fp->events+events<emax)?(fp->events+events):emax;

   change_cell(fp, fp->seed_i, fp->seed_j, fp->seed_n);

   while (fp->events < emax && 
          (tmax==0 || fp->t < tmax) && 
          (smax==0 || fp->stat_a-fp->stat_d < smax)) {

    if (DEBUG==2) {
     for (i=0; i<8; i++) {
       for (j=0;j<8;j++) printf("%5d",fp->Cell(i,j)); printf("\n");
     } printf("\n");
     for (i=0; i<8; i++) {
       for (p=3; p>=0; p--){
         for (j=0;i<(1<<p) && j<(1<<p);j++) 
           printf("%5d",fp->empty[p][i][j]); printf(" ");
       } printf("\n");
     } printf("\n");
     for (i=0; i<8; i++) {
       for (p=3; p>=0; p--){
         for (j=0;i<(1<<p) && j<(1<<p);j++) 
           printf("%5.1f",log(fp->rate[p][i][j])); printf(" ");
       } printf("\n");
     } printf("\n");
    }

     /* let the designated seed site wander around */
     /* must do this very frequently, else treadmilling would get stuck */
     if (wander && (fp->events%10)==0) {  
      int new_i, new_j, new_n, old_i, old_j, old_n, size=(1<<fp->P);
      new_i = fp->seed_i-1+random()%3;
      new_j = fp->seed_j-1+random()%3;
      if (periodic  || (new_i>=0 && new_i<size && new_j>=0 && new_j<size)) {
         if (periodic) { new_i=(new_i+size)%size; new_j=(new_j+size)%size; }
         if ((new_n = fp->Cell(new_i,new_j)) != 0) { 
	    /* OK, remove both fellows, then replace them both. */
            /* this way, rates[] reflects that the seed can't go away. */
            old_i=fp->seed_i; old_j=fp->seed_j; old_n=fp->seed_n;
            change_cell(fp,old_i,old_j,0); 
            change_cell(fp,new_i,new_j,0); 
            fp->seed_n = new_n; fp->seed_i=new_i; fp->seed_j=new_j;
            change_cell(fp,old_i,old_j,old_n); 
            change_cell(fp,new_i,new_j,new_n); 
            fp->events-=4; /* don't count this bookkeeping! */
            /* ERROR: stats are also modified!!! */
	 }
      }
     }

      dt = choose_cell(fp, &i, &j, &n); oldn = fp->Cell(i,j);
      if (i != fp->seed_i || j != fp->seed_j) {
        fp->t += dt;
        /* to add, must have se contact */
        if (oldn!=0 || HCONNECTED(fp,i,j,n)) change_cell(fp,i,j,n);
        if (n==0) {  /* dissociation: check connectedness */
           ringi = ((fp->Cell(i-1,j)!=0)<<7) +
             ((CONNECTED_W(fp,i-1,j+1) && CONNECTED_S(fp,i-1,j+1))<<6) +
                   ((fp->Cell(i,j+1)!=0)<<5) +
             ((CONNECTED_N(fp,i+1,j+1) && CONNECTED_W(fp,i+1,j+1))<<4) +
                   ((fp->Cell(i+1,j)!=0)<<3) +
             ((CONNECTED_E(fp,i+1,j-1) && CONNECTED_N(fp,i+1,j-1))<<2) +
                   ((fp->Cell(i,j-1)!=0)<<1) +
             ((CONNECTED_S(fp,i-1,j-1) && CONNECTED_E(fp,i-1,j-1))<<0);
           if (ring[ringi]==0) { /* couldn't quickly confirm... */
              flake_fission(fp,i,j);
           }
	} 
      } else dprintf("can't move seed!\n");
      d2printf("%d,%d -> %d\n",i,j,n);
   }
} 


/* for testing analytic solution to 2-tile 1D polymerization */
/* simulates until some limit is reached (all must be given) */
void linear_simulate(double ratek, double Gmc, double Gse, 
                     double tmax, int emax, int smax)
{
   double r,t; int e, s;
   unsigned char *tile;
   double k, rf, rr1, rr2, pr, pf;
   int i, errs; 

   if (tmax==0 || emax==0 || smax==0) return;

   tile = calloc(smax,sizeof(char));
   /* tile 1 = "A", tile 2 = "B", and we will start with "A" */

   rf = ratek*exp(-Gmc); rr1 = ratek*exp(-Gse); rr2 = ratek*exp(-2*Gse);
   e=0; t=0; s=1; tile[0]=1;

   while (e < emax && t < tmax && s < smax) {
     e++;
     if (s==1) {
        k = 2*rf; pr = 0; pf = 0.5; 
     } else if (tile[s-1]==tile[s-2]) {
        k = 2*rf+rr2; pr = rr2/k; pf=rf/k;
     } else {
        k = 2*rf+rr1; pr = rr1/k; pf=rf/k;
     }
     t += -log(drand48())/k;
    
     r = drand48();
     if (r<pr) s--; else if (r<pr+pf) tile[s++]=1; else tile[s++]=2;
   }

   for (errs=0,i=1; i<s; i++) if (tile[i-1] != tile[i]) errs++;
   fprintf(stdout,"%f %f %f   %f %d %d %d\n",
       Gmc,Gse,ratek, t, s, errs, e);

   free(tile);
}
