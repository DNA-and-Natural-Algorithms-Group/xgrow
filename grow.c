
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

int num_flakes=0;

/* sets up data structures for a flake -- cell field, hierarchical rates... */
flake *init_flake(unsigned char P, unsigned char N, 
  int seed_i, int seed_j, int seed_n, double Gfc)
{
  int i,j,p;
  int size = (1<<P);
  flake *fp = (flake *)malloc(sizeof(flake));

  fp->P = P; fp->N = N; 

  //  printf("Making flake %d x %d, %d tiles, seed=%d,%d,%d @ %6.2f\n",
  //         size,size,N,seed_i,seed_j,seed_n,Gfc);

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

  fp->flake_conc= (Gfc>0)?exp(-Gfc):0;
  fp->G=0; fp->mismatches=0; fp->tiles=0; fp->events=0;
  fp->seed_i=seed_i; fp->seed_j=seed_j; fp->seed_n=seed_n;
  fp->flake_ID = 0;  // until it's in a tube

  fp->next_flake=NULL; fp->tree_node=NULL; fp->tube=NULL;

  /* note that empty and rate are correct, because there are no tiles yet */

  change_cell(fp,seed_i,seed_j,seed_n);  
  /* not yet in tube: doesn't update tube stats */

  /* NOTE this routine is not "proper" -- assumes all allocs are OK */

  return fp;
}

/* returns the next flake in the list */
flake *free_flake(flake *fp)
{
  int i,p; flake *fpn;
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

  fpn=fp->next_flake; free(fp); 

  return fpn;
  /* NOTE because init_flake is not safe to out-of-mem, this could die */
}

/* for debugging purposes */
void print_tree(flake_tree *ftp, int L, char s)
{ int i;
 for (i=0;i<L;i++) printf(" ");
  printf("node %d %c (%d): empty %d, rate %g: L(%d) R(%d) U(%d)\n",
    L,s,(int)ftp,ftp->empty,ftp->rate,
    (int)ftp->left, (int)ftp->right, (int)ftp->up);
  if (ftp->left!=NULL) print_tree(ftp->left,L+1,'L');
  if (ftp->right!=NULL) print_tree(ftp->right,L+1,'R');
  if (ftp->left==NULL) { 
    for (i=0;i<=L;i++) printf(" ");
    printf("flake %d: %d tiles: tree_node(%d)\n",
	  (int)ftp->fp->flake_ID, ftp->fp->tiles, (int)ftp->fp->tree_node);
  }
}

/* gets rid of the tree nodes, without disturbing the flakes themselves */
void free_tree(flake_tree *ftp)
{
  if (ftp->left!=NULL) free_tree(ftp->left);
  if (ftp->right!=NULL) free_tree(ftp->right);
  free(ftp);
}

/* sets up data structures for tube -- tile set, params, scratch, stats  */
tube *init_tube(unsigned char P, unsigned char N, int num_bindings)
{
  int i,j,n,m;
  int size = (1<<P);
  tube *tp = (tube *)malloc(sizeof(tube));

  tp->P = P; tp->N = N; tp->num_bindings = num_bindings;
  tp->hydro=0;  tp->num_flakes=0;

  tp->units = (int**) calloc(sizeof(int*),N+1);
  for (i=0;i<=N;i++) {
   tp->units[i] = (int*) calloc(sizeof(int),4);
   for (j=0;j<4;j++) tp->units[i][j]=0;
  }
  tp->strength = (double*) calloc(sizeof(double),num_bindings+1);

  tp->conc = (double *)calloc(sizeof(double),N+1);
  for (n=0;n<N+1;n++) tp->conc[n]=0;
  tp->Gcb  = (double *)calloc(sizeof(double),N+1);
  for (n=0;n<N+1;n++) tp->Gcb[n]=0;
  tp->Gse_EW = (double **)calloc(sizeof(double *),N+1);
  for (n=0;n<N+1;n++) tp->Gse_EW[n]=(double *)calloc(sizeof(double),N+1);
  for (n=0;n<N+1;n++) for (m=0;m<N+1;m++) tp->Gse_EW[n][m]=0;
  tp->Gse_NS = (double **)calloc(sizeof(double *),N+1);
  for (n=0;n<N+1;n++) tp->Gse_NS[n]=(double *)calloc(sizeof(double),N+1);
  for (n=0;n<N+1;n++) for (m=0;m<N+1;m++) tp->Gse_NS[n][m]=0;

  tp->events=0; tp->t=0;
  tp->stat_a=tp->stat_d=tp->stat_h=tp->stat_f=0;

  tp->rv  = (double *)calloc(sizeof(double),N+1);
  tp->Fnext = (int *)calloc(sizeof(int),size*size);
  for (n=0;n<size*size;n++) tp->Fnext[n]=-1;
  tp->Fgroup = (int *)calloc(sizeof(int),size*size);

  tp->flake_list=NULL;
  tp->flake_tree=NULL;

  /* set_params() will have to put reasonable values in place */
  /* NOTE this routine is not "proper" -- assumes all allocs are OK */

  return tp;
}

void free_tube(tube *tp)
{
  int n; flake *fp;
  
  free(tp->conc);
  free(tp->Gcb);
  for (n=0;n<tp->N+1;n++) free(tp->Gse_EW[n]); free(tp->Gse_EW);
  for (n=0;n<tp->N+1;n++) free(tp->Gse_NS[n]); free(tp->Gse_NS);

  free(tp->rv); free(tp->Fgroup); free(tp->Fnext);
  
  for (n=0;n<tp->N+1;n++) free(tp->units[n]); free(tp->units);
  free(tp->strength); 

  free_tree(tp->flake_tree);  

  fp=tp->flake_list;
  while (fp!=NULL) { fp=free_flake(fp); } 

  free(tp);
  /* NOTE because init_flake is not safe to out-of-mem, this could die */
}


/* set up info for tile set, in flake data struc  */
/* fp->seed_n should have a defined value before entering set_params */
void set_params(tube *tp, int** units, double* strength, double* stoic,
 int hydro, double k, double Gmc, double Gse,
 double Gmch, double Gseh, double Ghyd, 
 double Gas, double Gam, double Gae, double Gah, double Gao, double T)
{
   int i,j,n,m;

   /* make our own copy of tile set and strengths, so the calling program
      can do as it pleases with it's units & strength information */
   for (i=0;i<=tp->N;i++) for (j=0;j<4;j++) tp->units[i][j] = units[i][j];
   for (i=0;i<=tp->num_bindings;i++) tp->strength[i] = strength[i]; 

   tp->hydro = hydro; 
   tp->T=T;

   tp->k = k;
   tp->kas = exp(-Gas);  tp->kao = exp(-Gao);
   tp->kam = exp(-Gam);  tp->kae = exp(-Gae); tp->kah = exp(-Gah);
   
   /* set tp->conc from Gmc... and from Gmch for hydrolysis rules */
   tp->conc[0]=0;  
   for (n=1; n <= tp->N; n++) 
      tp->conc[0]+=
        (tp->conc[n]=exp(-((n>tp->N/2&&tp->hydro)?Gmch:Gmc))*stoic[n]);

   /* set tp->Gcb from Ghyd */
   for (n=0; n <= tp->N; n++) tp->Gcb[n]=0;
   if (tp->hydro) for (n=tp->N/2+1; n <= tp->N; n++) tp->Gcb[n]=Ghyd;

   /* set Gse_EW Gse_NS from Gse, Gseh rules */
   /* uses (tp->units)[] and tp->strength[] */
   for (n=1; n<=tp->N; n++)
      for (m=1; m<=tp->N; m++) {
         tp->Gse_EW[n][m] = ((tp->units)[n][3]==(tp->units)[m][1]) *
              (tp->strength)[(tp->units)[m][1]] * 
                (tp->hydro?(((n > tp->N/2) || (m > tp->N/2))?Gseh:Gse):Gse);
         tp->Gse_NS[n][m] = ((tp->units)[n][2]==(tp->units)[m][0]) *
              (tp->strength)[(tp->units)[m][0]] * 
                (tp->hydro?(((n>tp->N/2 || m>tp->N/2))?Gseh:Gse):Gse);
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

/* recalculate flake energy & rates from scratch                    */
/* assume locations not adjacent to a tile correctly have zero rate */
/* BUG: if conc's go to zero, log returns nan                       */
void recalc_G(flake *fp)
{
  int n,i,j,size=(1<<fp->P);  tube *tp=fp->tube;

   /* don't count the seed tile concentration */
   fp->G=log(tp->conc[fp->seed_n]);  
   fp->mismatches=0; fp->tiles=0;
   /* add up all tile's entropy, bond energy, and hydrolysis energy */
   /* while we're at it, make sure 'mismatches' is correct          */
   /* and re-evaluate all off-rate & hydrolysis rates               */
   for (i=0;i<size;i++)
     for(j=0;j<size;j++) {
       if ((n=fp->Cell(i,j))>0) {
         fp->G += -log(tp->conc[n]) - Gse(fp,i,j,n)/2.0 - tp->Gcb[n];
         fp->mismatches += Mism(fp,i,j,n); fp->tiles++;
         update_rates(fp,i,j);
       } else if 
          ( fp->Cell(i+1,j) || fp->Cell(i,j+1) ||
            fp->Cell(i-1,j) || fp->Cell(i,j-1) ) {
         update_rates(fp,i,j);
       }
     }
   fp->mismatches/=2;
   update_tube_rates(fp);
}

/* calculate the dG of the flake, excluding concentration effects */
double calc_dG_bonds(flake *fp)
{
  int n,i,j,size=(1<<fp->P);  tube *tp=fp->tube;
  double dG=0;

   /* add up bond energy and hydrolysis energy only */
   for (i=0;i<size;i++)
     for(j=0;j<size;j++) {
       if ((n=fp->Cell(i,j))>0) {
         dG += - Gse(fp,i,j,n)/2.0 - tp->Gcb[n];
       }
     }
   return dG; 
}

/* calculate the perimeter of the flake */
int calc_perimeter(flake *fp)
{
  int n,i,j,size=(1<<fp->P); 
  int perimeter=0;

   /* add up number of empty cells next to this one */
   for (i=0;i<size;i++)
     for(j=0;j<size;j++) {
       if ((n=fp->Cell(i,j))>0) {
         perimeter += (fp->Cell(i-1,j)==0)+(fp->Cell(i+1,j)==0)+
                      (fp->Cell(i,j-1)==0)+(fp->Cell(i,j+1)==0);
       }
     }
   return perimeter;
}



void reset_params(tube *tp, double old_Gmc, double old_Gse, 
		  double new_Gmc, double new_Gse)
{  int n,m;
 flake *fp;

 if (!(tp->hydro)) {         /* not clear what to do for hydro rules */

   /* reset bond strengths */
   for (n=1; n<=tp->N; n++)
      for (m=1; m<=tp->N; m++) {
         tp->Gse_EW[n][m] = ((tp->units)[n][3]==(tp->units)[m][1]) *
              (tp->strength)[(tp->units)[m][1]] * new_Gse;
         tp->Gse_NS[n][m] = ((tp->units)[n][2]==(tp->units)[m][0]) *
              (tp->strength)[(tp->units)[m][0]] * new_Gse;
      }

   /* changing Gmc when Gfc>0 indicates either 
      dilution ( in which case all conc including Gfc decrease proportionally )
      concentration ( in which case same conc increase proportionally )
      Thus effects of increasing & decreasing are reversible
 
      adding tiles ( in which case all conc except Gfc increase additively ) 
      is currently not an option
   */
  tp->conc[0]=0;  
  for (n=1; n <= tp->N; n++) 
     tp->conc[0]+= (tp->conc[n]*=exp(-(new_Gmc-old_Gmc)));

  for (fp=tp->flake_list; fp!=NULL; fp=fp->next_flake) {
    fp->flake_conc*=exp(-(new_Gmc-old_Gmc));
    //    printf("\nPrior Params recalc_G(#%d)\n",fp->flake_ID);
    //             print_tree(tp->flake_tree,0,'*');  
    recalc_G(fp);
    //    printf("\nReset Params recalc_G(#%d)\n",fp->flake_ID);
    //             print_tree(tp->flake_tree,0,'*');  

  }
 }
}

/* put flake in the binary tree, somewhat balancing rates       */
/* (more principled would be to build Huffman tree bottom-up)   */
/* Note: tp->conc[0] & tp->k must be valid, to calculate rates. */
/* fp->rate and fp->empty must be non-zero for the same reason, */
/* hence recalc_G is used to update rates based on tube params. */
void insert_flake(flake *fp, tube *tp)
{
  int empty; double rate,kc; flake_tree *ftp, *ftpL, *ftpR;

  if (fp->N != tp->N || fp->P != tp->P) {
    printf("flake and tube incompatible!!\n"); exit(1);
  }

  fp->tube=tp; recalc_G(fp); 

  kc    = tp->k*tp->conc[0];
  rate  = fp->rate[0][0][0];
  empty = fp->empty[0][0][0];

  if (tp->flake_tree==NULL) {
    ftp = (flake_tree *)malloc(sizeof(flake_tree));
    ftp->left=ftp->right=NULL; ftp->fp=fp; ftp->up=NULL;
    ftp->empty=empty; ftp->rate=rate;
    tp->flake_tree=ftp;
    fp->tree_node=ftp;
  } else {
    ftp = tp->flake_tree;
    while (1) {
      if (ftp->fp!=NULL || (ftp->rate+kc*ftp->empty < rate+kc*empty)) {
         ftpL = (flake_tree *)malloc(sizeof(flake_tree));
         ftpL->rate=ftp->rate; ftpL->empty=ftp->empty; 
         ftpL->fp=ftp->fp; ftpL->left=ftp->left; ftpL->right=ftp->right;
         ftpL->up=ftp;
         if (ftp->fp!=NULL) ftp->fp->tree_node=ftpL;

         ftpR = (flake_tree *)malloc(sizeof(flake_tree));
         ftpR->rate=rate; ftpR->empty=empty; 
         ftpR->fp=fp; ftpR->left=ftpR->right=NULL;
         ftpR->up=ftp;
         fp->tree_node=ftpR;

         ftp->rate+=rate; ftp->empty+=empty;
         ftp->left=ftpL; ftp->right=ftpR; ftp->fp=NULL;
         
         break;
      } else if (ftp->left->rate+kc*ftp->left->empty < 
                 ftp->right->rate+kc*ftp->right->empty) {
        ftp->rate+=rate; ftp->empty+=empty;
	ftp=ftp->left;
      } else {
        ftp->rate+=rate; ftp->empty+=empty;
        ftp=ftp->right;
      }
    }
  }  
  fp->next_flake=tp->flake_list;
  tp->flake_list=fp;
  fp->flake_ID=++tp->num_flakes;
}


/* gives concentration-independent rates                                   */
/* for non-empty cells i,j, computes rate rv[n] to convert to type n       */
/*  so rv[0] gives the off-rate  (not the sum)                             */
/*  n=0 (empty) is the most common case.  the returned value is sum rv[n]  */
/* for empty cells i,j, returns 0                                          */
/* rv must already exist, of size fp->N+1                                  */
/* 0 <= i,j < 2^P                                                          */
double calc_rates(flake *fp, int i, int j, double *rv)
{
  int n,mi,ei,hi,mo,eo,ho; double r, sumr; tube *tp=fp->tube;
  unsigned char N,E,S,W;

  if (rv!=NULL) for (n=0;n<=fp->N;n++) rv[n]=0;
  if (tp==NULL) return 0;
  if (tp->T>0) return 0;   /* no off-rates: irreversible Tile Assembly Model */
  n = fp->Cell(i,j);
  if (n==0) return 0;                           /* no off-rate for empties   */
  if (i==fp->seed_i && j==fp->seed_j) return 0; /* seed site doesn't go away */

  r = tp->k * exp(-Gse(fp,i,j,n)); 
  if (rv!=NULL) rv[0]=r; sumr = r;

  /*  hydrolysis model assumes tiles 1...N/2 non hydrolyzed, N/2+1...N hydro */
  if (n > 0 && tp->hydro) {
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
    ei = ((tp->units)[n][2] != 0 && (tp->units)[S][0]==0) + 
         ((tp->units)[n][1] != 0 && (tp->units)[E][3]==0);
    mi = ((tp->units)[n][2] != 
           (tp->units)[S][0] && (tp->units)[n][2]*(tp->units)[S][0] !=0) +
         ((tp->units)[n][1] != 
           (tp->units)[E][3] && (tp->units)[n][1]*(tp->units)[E][3] !=0);
    hi = (E>fp->N/2) + (S>fp->N/2);

    /* empty and mismatched outputs */
    N=fp->Cell(i-1,j); W=fp->Cell(i,j-1);
    eo = ((tp->units)[n][0] != 0 && (tp->units)[N][2]==0) + 
         ((tp->units)[n][3] != 0 && (tp->units)[W][1]==0);
    mo = ((tp->units)[n][0] != 
           (tp->units)[N][2] && (tp->units)[n][0]*(tp->units)[N][2] !=0) +
         ((tp->units)[n][3] != 
           (tp->units)[W][1] && (tp->units)[n][3]*(tp->units)[W][1] !=0);
    ho = (N>fp->N/2) + (W>fp->N/2);

    r = tp->k * (tp->kas + mi * tp->kam + ei * tp->kae + hi * tp->kah +
                tp->kao * (mo * tp->kam + eo * tp->kae + ho * tp->kah));
    if (n<=fp->N/2) {
      if (rv!=NULL) rv[n+fp->N/2]=r; sumr += r;
    } else {
      r = r * exp(-tp->Gcb[n]+tp->Gcb[n-fp->N/2]
                  -Gse(fp,i,j,n)+Gse(fp,i,j,n-fp->N/2));
      if (rv!=NULL) rv[n-fp->N/2]=r; sumr += r;
    }
  } 

  return sumr;
}


/* figure the delta_rate for this cell, and propagate up.           */
/* note: as long as hierarchy is accurate w/r to contents,          */
/* this will make cell i,j & it's contributions correct,            */
/* so it can be used / needs to be done after param changes.        */
/* Mainly used after tile association / dissociation / hydrolysis:  */
/* Cell(i,j) has just changed to n, but old rates and empty status  */
/*   has not yet changed to reflect this.  Since the cell and its   */
/*   neighbors must all be changed, we only guarantee that THIS     */
/*   cell's hierarchical status has not been updated yet.           */
/*   ii,jj is cell i,j or a neighbor cell; possibly a border cell.  */
/*   Note: "empty" means "empty and neighbor to a tile".            */
/*   In irreversible Tile Assembly Model, "empty" means that there  */
/*   is(are) a tile type(s) that could be added at the location.    */
void update_rates(flake *fp, int ii, int jj)
{
   int n,p; int size=(1<<fp->P); tube *tp=fp->tube;
   double oldrate, newrate; int oldempty, newempty;

   if (periodic) { ii=(ii+size)%size; jj=(jj+size)%size; }

   if (!(ii < 0 || ii >= size || jj < 0 || jj >= size)) {
      oldrate  = fp->rate[fp->P][ii][jj];
      oldempty = fp->empty[fp->P][ii][jj];
      newempty = (fp->Cell(ii,jj)==0) &&
                ( fp->Cell(ii+1,jj) || fp->Cell(ii,jj+1) ||
                  fp->Cell(ii-1,jj) || fp->Cell(ii,jj-1) );
      if (newempty && tp!=NULL && tp->T>0) { 
        /* calculate how many tile types could make >= T bonds */
	for (n=1;n<=fp->N;n++)
          if (Gse(fp,ii,jj,n)>=tp->T) newempty++;
        newempty--;  newempty=(newempty>0);  /* only care if exists */
      }
      newrate  = calc_rates(fp, ii, jj, NULL);
      for (p=fp->P; p>=0; p--) {
         fp->rate[p][ii][jj] += newrate-oldrate;
         fp->empty[p][ii][jj] += newempty-oldempty;
	 if (p<fp->P)
     /* always fix-up any numerical error that could have accumulated here */
     fp->rate[p][ii][jj] = fp->rate[p+1][2*ii][2*jj]   +
                           fp->rate[p+1][2*ii][2*jj+1] +
                           fp->rate[p+1][2*ii+1][2*jj] +
                           fp->rate[p+1][2*ii+1][2*jj+1];

         ii = (ii>>1); jj = (jj>>1);
      }
   }
}

void update_tube_rates(flake *fp)
{
  flake_tree *ftp=fp->tree_node; 
  double oldrate, newrate; int oldempty, newempty;
 
  if (ftp==NULL) return;
 
  oldrate = ftp->rate; oldempty=ftp->empty;
  newrate = fp->rate[0][0][0]; newempty=fp->empty[0][0][0];

  while (ftp!=NULL) {
    ftp->rate+=newrate-oldrate;
    ftp->empty+=newempty-oldempty;
    ftp=ftp->up;
  }
  
}


/* convert Cell(i,j) to type n.                                       */
/* update all hierarchical rates and empty counts, in flake and tube. */
/* ALL modifications to the cell array go through this interface!     */
/* we only modify cells with 0 <= i,j < (1<<fp->P)                    */
/* but since it might be called out of range in FILL, we fix it up.   */
/* BUG: changes in concentration should change G for every tile, but  */
/* we don't update G; also, if conc[n]==0, nan is result.             */
void change_cell(flake *fp, int i, int j, unsigned char n)
{
   int size=(1<<fp->P);  tube *tp=fp->tube; 
   if (periodic) { i=(i+size)%size; j=(j+size)%size; }

   if (fp->Cell(i,j)==n) return;  // nothing to change!

   if (tp!=NULL) { /* flake has been added to a tube */
    if (fp->Cell(i,j)==0) {                         /* tile addition */
     if (tp->conc[n]<=fp->flake_conc) return; // conc's can't got to zero!
     fp->G += -log(tp->conc[n]) - Gse(fp,i,j,n);   
     tp->conc[n] -= fp->flake_conc; 
     tp->conc[0] -= fp->flake_conc;
     tp->stat_a++; fp->tiles++; 
     fp->mismatches += Mism(fp,i,j,n);
    } else if (n==0) {                              /* tile loss */
     tp->conc[0] += fp->flake_conc; 
     tp->conc[fp->Cell(i,j)] += fp->flake_conc; 
     fp->G += log(tp->conc[fp->Cell(i,j)]) + Gse(fp,i,j,fp->Cell(i,j));
     tp->stat_d++; fp->tiles--; 
     fp->mismatches -= Mism(fp,i,j,fp->Cell(i,j));
    } else {                                        /* tile hydrolysis */
     fp->G += Gse(fp,i,j,fp->Cell(i,j)) - Gse(fp,i,j,n) +
              tp->Gcb[fp->Cell(i,j)] - tp->Gcb[n]; 
     tp->stat_h++; 
     /* by our rules, hydrolyzed tiles have same se types as non-hyd. */
    }
    tp->events++; fp->events++;
   }

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

   if (tp!=NULL) update_tube_rates(fp);
}


/* use rates & empty & conc to choose a cell to change,      */
/* and call calc_rates to identify what change to make.      */
/* report choice, but don't act on it.                       */
void choose_cell(flake *fp, int *ip, int *jp, int *np)
{
  double sum,cum,r,kc,k00,k01,k10,k11;
  int p,i,j,di=1,dj=1,n,oops;   tube *tp=fp->tube;


  kc = tp->k*tp->conc[0];
  sum = fp->rate[0][0][0] + kc*fp->empty[0][0][0];

  i=0; j=0;  r=drand48();  // we'll re-use this random number for all levels
  for (p=0; p<fp->P; p++) { /* choosing subquadrant from within p:i,j */
     k00 = fp->rate[p+1][2*i][2*j]+kc*fp->empty[p+1][2*i][2*j];
     k10 = fp->rate[p+1][2*i+1][2*j]+kc*fp->empty[p+1][2*i+1][2*j];
     k01 = fp->rate[p+1][2*i][2*j+1]+kc*fp->empty[p+1][2*i][2*j+1];
     k11 = fp->rate[p+1][2*i+1][2*j+1]+kc*fp->empty[p+1][2*i+1][2*j+1];
     sum = (k00+k01+k10+k11);  
     /* avoid possible round-off error... but still check for it */
     d2printf("%f / %f for choosing %d from %d: %d %d\n",r,sum,p+1,p,i,j);
     do {
      r = r*sum;  oops=0;
      if ( (r-=k00) < 0) { di=0; dj=0; r=(r+k00)/k00; } else
      if ( (r-=k10) < 0) { di=1; dj=0; r=(r+k10)/k10; } else
      if ( (r-=k01) < 0) { di=0; dj=1; r=(r+k01)/k01; } else
      if ( (r-=k11) < 0) { di=1; dj=1; r=(r+k11)/k11; } else 
                         { r=drand48(); oops=1; }
     } while (oops);
     /* always fix-up any numerical error that could have accumulated here */
     fp->rate[p][i][j] = fp->rate[p+1][2*i][2*j]+fp->rate[p+1][2*i][2*j+1]+
                         fp->rate[p+1][2*i+1][2*j]+fp->rate[p+1][2*i+1][2*j+1];
     i=2*i+di; j=2*j+dj;
  }
  *ip=i; *jp=j;
  // upon exit, we should still have a good random number r

  if (fp->Cell(i,j) == 0) {   /* choose on-event for type 1...N            */
     do {
       r = r * tp->conc[0];  cum = 0;  oops=0;
       for (n=1; n<=fp->N; n++) if (r < (cum += tp->conc[n])) break; 
       if (n>fp->N) { // apparently conc[0] is not the sum of conc[n], oops
	 printf("Eye of the needle!!! %f =!= %f\n",tp->conc[0],cum); 
         r=drand48(); oops=1; 
         tp->conc[0]=0; for (n=1; n <= tp->N; n++) tp->conc[0]+=tp->conc[n];
       }
     } while (oops);
  } else {                    /* choose off-event 0 or conversion to 1...N */
    if (tp->hydro) {
       sum = calc_rates(fp,i,j,tp->rv);
       if (sum==0) printf("Zero-sum game!!!\n");
       r = r * sum;  cum = 0;
       for (n=0; n<=fp->N; n++) if (r < (cum += tp->rv[n])) break; 
       if (n>fp->N) { printf("A rose is not a rose!!!\n"); n=0; }
    } else {
       n=0;  // always an off-event, unless hydrolysis rules are used.
    }
  }
  *np = n;
}

flake *choose_flake(tube *tp)
{
  double r,kc,kL,kR;  int oops;   
  flake_tree *ftp=tp->flake_tree; 

  kc = tp->k*tp->conc[0];

  r=drand48();  // we'll re-use this random number for all levels
  while (ftp->fp==NULL) {
     kL = ftp->left->rate+kc*ftp->left->empty;
     kR = ftp->right->rate+kc*ftp->right->empty;
     /* always fix-up any numerical error that could have accumulated here */
     ftp->rate = ftp->left->rate+ftp->right->rate;
     ftp->empty = ftp->left->empty+ftp->right->empty; // shouldn't be necessary
     do {
      r = r*(kL+kR);  oops=0;
      if ( (r-=kL) < 0) { ftp=ftp->left; r=(r+kL)/kL; } else
      if ( (r-=kR) < 0) { ftp=ftp->right; r=(r+kR)/kR; } else
                        { r=drand48(); oops=1; }
     } while (oops);
  }
  // upon exit, ftp is now a leaf; ftp->fp is our chosen flake
  return ftp->fp;
}

/*************** the fill routine for checking connectedness ***********/

/* A tile is defined to be HCONNECTED to a neighbor if they both have   */
/* non-zero se types  (who cares what the Gse is).                      */
/* This is hypothetical on i,j being tile n != 0.                       */
/* CONNECTED is the non-hypothetical version.                           */
#define HCONNECTED_N(fp,i,j,n) \
    ((fp->tube->units)[n][0]!=0 && (fp->tube->units)[fp->Cell((i)-1,j)][2]!=0)
#define HCONNECTED_E(fp,i,j,n) \
    ((fp->tube->units)[n][1]!=0 && (fp->tube->units)[fp->Cell(i,(j)+1)][3]!=0)
#define HCONNECTED_S(fp,i,j,n) \
    ((fp->tube->units)[n][2]!=0 && (fp->tube->units)[fp->Cell((i)+1,j)][0]!=0)
#define HCONNECTED_W(fp,i,j,n) \
    ((fp->tube->units)[n][3]!=0 && (fp->tube->units)[fp->Cell(i,(j)-1)][1]!=0)
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
/* The queue is in the array tp->Fnext[n], storing linked lists such     */
/* that no cell is in more than one linked list (careful w/ Fpush!), and */
/* the fill memory is in tp->Fgroup[n] identifying the group class.      */
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
  else { tp->Fnext[tail[g]]=Qn(i,j); tail[g]=Qn(i,j); }     \
  if ((((i)+size)%size)==fp->seed_i &&                      \
      (((j)+size)%size)==fp->seed_j) seeded[ming[g]]=1;     \
  tp->Fgroup[Qn(i,j)]=g;}
#define Fpull(g,i,j) { int oldh=head[g];                    \
  i=Qi(head[g]); j=Qj(head[g]);                             \
  head[g]=tp->Fnext[head[g]]; tp->Fnext[oldh]=-1;           \
  if (Fempty(g)) tail[g]=-1; }
#define Fmerge(g1,g2) {                                                   \
  if (ming[g1]!=ming[g2]) {                                               \
     int gi,newg = MIN(ming[g1],ming[g2]),oldg = MAX(ming[g1],ming[g2]);  \
     for (gi=1;gi<5;gi++) if (ming[gi]==oldg) ming[gi]=newg;              \
     seeded[newg] = (seeded[newg]||seeded[oldg]); ngroups--;              \
  } }


/* here, we loose the part that's unconnected to the seed, rather     */
/* than saving it and creating a new flake                            */
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
  tube *tp=fp->tube;
 
  for (g=0;g<5;g++) { head[g]=tail[g]=-1; ming[g]=-1; seeded[g]=0; }
  ming[0]=0;

  if (DEBUG==2) { int n; 
     for (i=j=n=0;n<size*size;n++) 
         { i+=(tp->Fgroup[n] != 0); j+=(tp->Fnext[n] != -1); }
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
          if ((gg=tp->Fgroup[Qn(i,j+1)])==0) { Fpush(ming[g],i,j+1) }
          else { Fmerge(g,gg) }
        }
        if (CONNECTED_S(fp,i,j)) {
          if ((gg=tp->Fgroup[Qn(i+1,j)])==0) { Fpush(ming[g],i+1,j) }
          else { Fmerge(g,gg) }
        }
        if (CONNECTED_W(fp,i,j)) {
          if ((gg=tp->Fgroup[Qn(i,j-1)])==0) { Fpush(ming[g],i,j-1) }
          else { Fmerge(g,gg) }
        }
        if (CONNECTED_N(fp,i,j)) {
          if ((gg=tp->Fgroup[Qn(i-1,j)])==0) { Fpush(ming[g],i-1,j) }
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
    }  /* tp->Fnext is now all -1, as are head & tail */

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
       if (tp->Fgroup[Qn(i,j+1)]>0 && seeded[ming[tp->Fgroup[Qn(i,j+1)]]]==0) 
          { change_cell(fp,i,j+1,0); Fpush(0,i,j+1) }
       if (tp->Fgroup[Qn(i+1,j)]>0 && seeded[ming[tp->Fgroup[Qn(i+1,j)]]]==0)
          { change_cell(fp,i+1,j,0); Fpush(0,i+1,j) }
       if (tp->Fgroup[Qn(i,j-1)]>0 && seeded[ming[tp->Fgroup[Qn(i,j-1)]]]==0) 
          { change_cell(fp,i,j-1,0); Fpush(0,i,j-1) }
       if (tp->Fgroup[Qn(i-1,j)]>0 && seeded[ming[tp->Fgroup[Qn(i-1,j)]]]==0)
          { change_cell(fp,i-1,j,0); Fpush(0,i-1,j) }
    }
    /* again tp->Fnext is now all -1, as are head & tail */
  }    

  /* now re-zero Fgroup for stuff that sticks around */
  i=ii; j=jj;
  if (fp->Cell(i,j+1)!=0) { Fpush(0,i,j+1) }
  if (fp->Cell(i+1,j)!=0) { Fpush(0,i+1,j) }
  if (fp->Cell(i,j-1)!=0) { Fpush(0,i,j-1) }
  if (fp->Cell(i-1,j)!=0) { Fpush(0,i-1,j) }
  while (!Fempty(0)) {
     Fpull(0,i,j)
     if (tp->Fgroup[Qn(i,j+1)]>0) { Fpush(0,i,j+1) }
     if (tp->Fgroup[Qn(i+1,j)]>0) { Fpush(0,i+1,j) }
     if (tp->Fgroup[Qn(i,j-1)]>0) { Fpush(0,i,j-1) }
     if (tp->Fgroup[Qn(i-1,j)]>0) { Fpush(0,i-1,j) }
  }
  /* now tp->Fnext is all -1 and tp->Fgroup is all 0 */
  tp->stat_f += ngroups-1;
}



/* remove all tiles whose off-rate is 'X' times faster than its on-rate. */
/* repeat 'iters' times. */
void clean_flake(flake *fp, double X, int iters)
{
  int i,j; double kc;  tube *tp=fp->tube;
  int size = (1<<fp->P); int it;

  kc = tp->k*tp->conc[0]; /* on-rate */

  /* first memorize, then remove, to avoid changing rates during removal */
  for (it=0; it<iters; it++) {
    for (i=0; i<size; i++)
      for (j=0; j<size; j++)
        tp->Fgroup[i+size*j] = 
          (fp->rate[fp->P][i][j] > X * tp->k * tp->conc[fp->Cell(i,j)]);
    for (i=0; i<size; i++)
      for (j=0; j<size; j++)
        if (tp->Fgroup[i+size*j]) change_cell(fp, i, j, 0);
    for (i=0; i<size; i++)
      for (j=0; j<size; j++)
        tp->Fgroup[i+size*j] = 0;
  }

}

/* simulates 'events' events */
void simulate(tube *tp, int events, double tmax, int emax, int smax)
{
  int i,j,n,oldn; double dt; flake *fp=tp->flake_list;
  unsigned char ringi;  double total_rate; long int emaxL;

  if (fp==NULL) return;  /* no flakes! */

   emaxL = (emax==0 || tp->events+events<emax)?(tp->events+events):emax;

   /* make sure seed tiles are in place (is this really needed?) */
   while (fp!=NULL) {  
     change_cell(fp, fp->seed_i, fp->seed_j, fp->seed_n);
     fp=fp->next_flake;
   }

   total_rate = tp->flake_tree->rate+tp->k*tp->conc[0]*tp->flake_tree->empty;
   while (tp->events < emaxL && 
          (tmax==0 || tp->t < tmax) && 
          (smax==0 || tp->stat_a-tp->stat_d < smax) &&
          total_rate > 0) {

     dt = -log(drand48()) / total_rate;

     fp=choose_flake(tp);

     /* let the designated seed site wander around */
     /* must do this very frequently, else treadmilling would get stuck */
     if (wander && (tp->events%10)==0) {  
      int new_i, new_j, new_n, old_i, old_j, old_n, size=(1<<tp->P);
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
            tp->events-=4; /* don't count this bookkeeping! */
            /* ERROR: stats are also modified!!! */
	 }
      }
     }

     choose_cell(fp, &i, &j, &n); oldn = fp->Cell(i,j);
     if (i != fp->seed_i || j != fp->seed_j) { /* can't change seed tile */
       tp->t += dt;
       if (tp->T>0) { /* irreversible Tile Assembly Model */
        if (n>0 && Gse(fp,i,j,n)>=tp->T) change_cell(fp,i,j,n);
        else { tp->stat_a++; tp->stat_d++; tp->events+=2; fp->events+=2; } 
        /* NOTE this is extremely inefficient -- would be much better to
        actually change the rates such that bad on-events and all off-events
        have rate 0.  However, this doesn't fit well into the "empty"
        framework for calculating rates.  Another fast solution would be to
        keep a linked list of active growth sites.  Alas. */
        /* NOTE the above refered to before update_rates() and calc_rates() 
        were modified for the T>0 case.  Now it's faster.  Still, if 
        only one of N tile types can be added at a location, all have
        equal chance to be chosen & rejected -- this is inefficient, but 
        unless fp->empty were to be changed, we can't handle unequal tile
        concentrations. */
       } else {
        /* to add, must have se contact; hydrolysis also happens here */
        if (oldn!=0 || HCONNECTED(fp,i,j,n)) change_cell(fp,i,j,n);
        /* zero-bond tile additions fall off immediately;
           count them or else, if there are lots, the display
           can be super-slow! */
        if (oldn==0 && ~HCONNECTED(fp,i,j,n)) 
           { tp->stat_a++; tp->stat_d++; tp->events+=2; fp->events+=2; } 
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
       }
      } else dprintf("can't move seed!\n");
      d2printf("%d,%d -> %d\n",i,j,n);
     total_rate = tp->flake_tree->rate+tp->k*tp->conc[0]*tp->flake_tree->empty;
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
