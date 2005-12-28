/* grow.c

This code is freely distributable.

by Erik Winfree


*/
  

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <assert.h>
# include <limits.h>
# include <unistd.h>
# include <string.h>

# include "grow.h"
# include "xgrow-tests.h"

/* index is 8 bits N NE E SE S SW W NW where                         */
/*  N  E  S  W    refer to where there is a tile present, and        */
/*   NE SE SW NW  indicate whether the corner tile connects its      */
/*                neighbors (given that all three are present)       */
/* then ring[index]=1 if there is 1 (or fewer) connected groups.     */
unsigned char ring[256];
#define ROTATE(i)          ((((i)&1)<<7) + ((i)>>1))
#define ROTATE_CLEAR(i)    (               ((i)>>1))
#define AVOGADROS_NUMBER 6.022e23

/*************** the fill routine for checking connectedness ***********/

#if 0
/* original (pre-Dec 29, 2003) definitions:                             */
/* A tile is defined to be HCONNECTED to a neighbor if they both have   */
/* non-zero se types  (who cares what the Gse is).                      */
/* This is hypothetical on i,j being tile n != 0.                       */
/* CONNECTED is the non-hypothetical version.                           */
#define HCONNECTED_N(fp,i,j,n) \
    ((fp->tube->tileb)[n][0]!=0 && (fp->tube->tileb)[fp->Cell((i)-1,j)][2]!=0)
#define HCONNECTED_E(fp,i,j,n) \
    ((fp->tube->tileb)[n][1]!=0 && (fp->tube->tileb)[fp->Cell(i,(j)+1)][3]!=0)
#define HCONNECTED_S(fp,i,j,n) \
    ((fp->tube->tileb)[n][2]!=0 && (fp->tube->tileb)[fp->Cell((i)+1,j)][0]!=0)
#define HCONNECTED_W(fp,i,j,n) \
    ((fp->tube->tileb)[n][3]!=0 && (fp->tube->tileb)[fp->Cell(i,(j)-1)][1]!=0)
#else
/* new (post-Dec 29, 2003) definitions:                                 */
/* A tile is defined to be HCONNECTED to a neighbor if the              */
/* Gse of the bond between them is non-zero (regardless of how small).  */
/* This is hypothetical on i,j being tile n != 0.                       */
/* CONNECTED is the non-hypothetical version.                           */
#define HCONNECTED_N(fp,i,j,n) (fp->tube->Gse_NS[fp->Cell((i)-1,j)][n]>0)
#define HCONNECTED_E(fp,i,j,n) (fp->tube->Gse_EW[fp->Cell(i,(j)+1)][n]>0)
#define HCONNECTED_S(fp,i,j,n) (fp->tube->Gse_NS[n][fp->Cell((i)+1,j)]>0)
#define HCONNECTED_W(fp,i,j,n) (fp->tube->Gse_EW[n][fp->Cell(i,(j)-1)]>0)
#endif

#define HCONNECTED(fp,i,j,n) \
    ( HCONNECTED_N(fp,i,j,n) || HCONNECTED_E(fp,i,j,n) || \
      HCONNECTED_S(fp,i,j,n) || HCONNECTED_W(fp,i,j,n) )
#define CONNECTED_N(fp,i,j) HCONNECTED_N(fp,i,j,fp->Cell(i,j))
#define CONNECTED_E(fp,i,j) HCONNECTED_E(fp,i,j,fp->Cell(i,j))
#define CONNECTED_S(fp,i,j) HCONNECTED_S(fp,i,j,fp->Cell(i,j))
#define CONNECTED_W(fp,i,j) HCONNECTED_W(fp,i,j,fp->Cell(i,j))
#define CONNECTED(fp,i,j) HCONNECTED(fp,i,j,fp->Cell(i,j))
/* all of these can only be used for 0 <= i,j < size                     */

int num_flakes=0;
int double_tile_count=0;
flake *blank_flakes = NULL;

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
  if (ftp != NULL) {
    if (ftp->left!=NULL) free_tree(ftp->left);
    if (ftp->right!=NULL) free_tree(ftp->right);
    free(ftp);
  }
}

/* sets up data structures for tube -- tile set, params, scratch, stats  */
tube *init_tube(unsigned char P, unsigned char N, int num_bindings)
{
  int i,j,n,m;
  int size = (1<<P);
  tube *tp = (tube *)malloc(sizeof(tube));

  tp->P = P; tp->N = N; tp->num_bindings = num_bindings;
  tp->hydro=0;  tp->num_flakes=0; tp->total_flakes = 0;
  tp->largest_flake_size = 0;

  tp->tileb = (int**) calloc(sizeof(int*),N+1);
  for (i=0;i<=N;i++) {
    tp->tileb[i] = (int*) calloc(sizeof(int),4);
    for (j=0;j<4;j++) tp->tileb[i][j]=0;
  }
  tp->strength = (double*) calloc(sizeof(double),num_bindings+1);
  tp->glue = (double **) calloc(sizeof(double*),num_bindings+1);
  for (i=0;i<=num_bindings;i++) {
    tp->glue[i] = (double*) calloc(sizeof(double),num_bindings+1);
    for (j=0;j<=num_bindings;j++) { tp->glue[i][j]=0; }  // necessary??  -- EW
  }

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

  tp->events=0; tp->t=0; tp->ewrapped=0;
  tp->stat_a=tp->stat_d=tp->stat_h=tp->stat_f=0;

  tp->rv  = (double *)calloc(sizeof(double),1+N+4);
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
  int n,i; flake *fp;
  
  free(tp->conc);
  free(tp->Gcb);
  for (n=0;n<tp->N+1;n++) free(tp->Gse_EW[n]); free(tp->Gse_EW);
  for (n=0;n<tp->N+1;n++) free(tp->Gse_NS[n]); free(tp->Gse_NS);

  free(tp->rv); free(tp->Fgroup); free(tp->Fnext);
  
  for (n=0;n<tp->N+1;n++) free(tp->tileb[n]); free(tp->tileb);
  for (i=0;i<tp->num_bindings+1;i++) free(tp->glue[i]); free(tp->glue);
  free(tp->strength); 

  free_tree(tp->flake_tree);  

  fp=tp->flake_list;
  while (fp!=NULL) { fp=free_flake(fp); } 

  free(tp);
  /* NOTE because init_flake is not safe to out-of-mem, this could die */
}

void set_Gses(tube *tp, double Gse, double Gseh) {
  int n,m;
  for (n=1; n<=tp->N; n++)
    for (m=1; m<=tp->N; m++) {
      
      tp->Gse_EW[n][m] = ((((tp->tileb)[n][3]==(tp->tileb)[m][1]) *
			   (tp->strength)[(tp->tileb)[m][1]]) +
			  (tp->glue)[(tp->tileb)[n][3]][(tp->tileb)[m][1]])*
	(tp->hydro?(((n>tp->N/2 || m>tp->N/2))?Gseh:Gse):Gse);
      tp->Gse_NS[n][m] = (((tp->tileb)[n][2]==(tp->tileb)[m][0]) *
			  (tp->strength)[(tp->tileb)[m][0]] +
			  (tp->glue)[(tp->tileb)[n][2]][(tp->tileb)[m][0]])* 
	(tp->hydro?(((n>tp->N/2 || m>tp->N/2))?Gseh:Gse):Gse);
    }
}


/* set up info for tile set, in flake data struc  */
/* fp->seed_n should have a defined value before entering set_params */
void set_params(tube *tp, int** tileb, double* strength, double **glue, double* stoic,
		double anneal_g, double anneal_t, int updates_per_RC,
		int *dt_right, int *dt_left, int hydro, double k, double Gmc, double Gse,
		double Gmch, double Gseh, double Ghyd, 
		double Gas, double Gam, double Gae, double Gah, double Gao, double T,
		double tinybox, int seed_i, int seed_j, double Gfc)
{
  int i,j,n;
  /* make our own copy of tile set and strengths, so the calling program
     can do as it pleases with it's tileb & strength information */
  for (i=0;i<=tp->N;i++) for (j=0;j<4;j++) tp->tileb[i][j] = tileb[i][j];
  for (i=0;i<=tp->num_bindings;i++) tp->strength[i] = strength[i]; 
  for (i=0;i<=tp->num_bindings;i++) {
    for (j=0;j<=tp->num_bindings;j++) {
      tp->glue[i][j] = glue[i][j];
    }
  }

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
  if (anneal_t && Gse < anneal_g) {
    fprintf(stderr,"Final Gse must be larger than initial Gse for an anneal.\n");
    exit(-1);
  }
  tp->Gmc = Gmc;
  tp->anneal_g = anneal_g;
  tp->anneal_t = anneal_t;
  tp->Gse_final = Gse;
  if (tp->anneal_t) {
    tp->Gse = anneal_g;
  } else {
    tp->Gse = Gse;
  }
  tp->updates = 1;
  tp->update_freq = updates_per_RC;
  tp->next_update_t = exp(-(((double) tp->updates)*tp->anneal_t*log(2))/tp->update_freq);
  tp->dt_right = dt_right;
  tp->dt_left = dt_left;
  for (n=0; n< tp->N; n++) {
    if (tp->dt_right[n]) {
      double_tile_count++;
    }
  }
  tp->tinybox = tinybox;
  tp->default_seed_i = seed_i;
  tp->default_seed_j = seed_j;
  tp->initial_Gfc = Gfc;

  /* set tp->Gcb from Ghyd */
  for (n=0; n <= tp->N; n++) tp->Gcb[n]=0;
  if (tp->hydro) for (n=tp->N/2+1; n <= tp->N; n++) tp->Gcb[n]=Ghyd;

  /* set Gse_EW Gse_NS from Gse, Gseh rules */
  /* uses (tp->tileb)[] and tp->strength[] and tp->glue[] */
  /* XXX We will want to modify this */  /* See also (change also!) reset_params */
  set_Gses(tp,Gse,Gseh);
  tp->watching_states  = 0;
  tp->tracking_seen_states = 0;
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
} // set_params()

/* recalculate flake energy & rates from scratch                    */
/* assume locations not adjacent to a tile correctly have zero rate */
/* BUG: if conc's go to zero, log returns nan                       */
void recalc_G(flake *fp)
{
  int n,i,j,size=(1<<fp->P);  tube *tp=fp->tube;

  fp->G = 0; fp->mismatches=0; fp->tiles=0;
  /* OLD: don't count the seed tile concentration */
  /* NEW: count it only if 'wander' is on -- but only in the display */

  /* add up all tile's entropy, bond energy, and hydrolysis energy */
  /* while we're at it, make sure 'mismatches' is correct          */
  /* and re-evaluate all off-rate & hydrolysis rates               */
  for (i=0;i<size;i++)
    for(j=0;j<size;j++) {
      if ((n=fp->Cell(i,j))>0) {
	fp->G += -log(tp->conc[n]) - Gse(fp,i,j,n)/2.0 - tp->Gcb[n];
	/* Don't count Gmc for both sides of a double tile */
	/* Don't count Gse between double tile */
	if (tp->dt_right[n]) {
	  fp->G += log(tp->conc[n]);
	  fp->G += fp->tube->Gse_EW[fp->Cell(i,(j+1)%size)][n];
	}
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
} // recalc_G()

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
} // calc_dG_bonds()

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
} // calc_perimeter()


void reset_params(tube *tp, double old_Gmc, double old_Gse, 
		  double new_Gmc, double new_Gse, double Gseh)
{  int n;
 flake *fp;

 if (!(tp->hydro)) {         /* not clear what to do for hydro rules */

   /* reset bond strengths for new Gse  (this code same as in set_params) */

   set_Gses(tp,new_Gse,Gseh);

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
} // reset_params()

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

  /* If the flake_tree is empty, make the root. */
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
  fp->flake_ID=++tp->total_flakes;
  tp->num_flakes++;
} // insert_flake()

void add_flake_to_reserve_list(flake *fp) {
  int p, i, j;
  int size;
  // First clear flake
  size = (1<< (fp->P));
  for (i=0;i<2+size;i++) 
    memset(fp->cell[i],0,2+size*sizeof(char));
  for (p=0;p<=fp->P;p++) {
    size = (1<<p);
    for (i=0;i<size;i++) {
      //memset(fp->empty[p][i],size,sizeof(int));      
      for (j=0;j<size;j++) {
	fp->empty[p][i][j]=0;
	fp->rate[p][i][j]=0;
      }
    }
  }
  fp->G=0; fp->mismatches=0; fp->tiles=1; fp->events=0;
  fp->tree_node = NULL;
  fp->next_flake = blank_flakes;
  blank_flakes = fp;
}

/* Returns a flake from the list of blank flakes, if any are available */
flake * recover_flake (int seed_i, int seed_j, int seed_n, int Gfc) {
  flake *fp;
  if (blank_flakes) {
    fp = blank_flakes;
    fp->seed_i = seed_i;
    fp->seed_j = seed_j;
    fp->seed_n = seed_n;
    fp->flake_conc= (Gfc>0)?exp(-Gfc):0;
    change_cell(fp,seed_i,seed_j,seed_n);  

    blank_flakes = blank_flakes->next_flake;
    return fp;
  }
  else {
    return NULL;
  }
}

/* This function is used with the tinybox option to remove a flake
   that has become just a single tile.  It doesn't try to balance
   anything for now, and since flakes will be entering and leaving
   quickly, the tree should be rebalanced again soon.  */
void remove_flake(flake *fp) {
  tube *tp;
  flake_tree *ftp, *tmp, *u;
  flake *f;

  tp=fp->tube;
  ftp = fp->tree_node;
  assert (ftp != NULL);
  u = ftp->up;
  if (u) {
    /* Get rid of the node in the tree and merge the other child with
     the parent*/
    if (u->right == ftp) {
      tmp = u->left;
    }
    else {
      assert (u->left == ftp);
      tmp = u->right;
    }
    u->fp = tmp->fp;
    if (u->fp) {
      u->fp->tree_node = u;
    }
    u->rate = tmp->rate;
    u->empty = tmp->empty;
    u->right = tmp->right;
    if (u->right) {
      u->right->up = u;
    }
    u->left = tmp->left;
    if (u->left) { 
      u->left->up = u;
    }
    free (tmp);
    /* Now propagate the rate change up the tree */
    u = u->up;
    while (u) {
      //u->rate -= ftp->rate;
      //u->empty -= ftp->empty;
      u->rate = u->left->rate + u->right->rate;
      u->empty = u->left->empty + u->right->empty;
      u = u->up;
    }
  }
  else {
    /* Tree is now empty */
    tp->flake_tree = NULL;
  }
  free (ftp);
  // Remove the flake from the flake list
  if (fp == tp->flake_list) {
    tp->flake_list = fp->next_flake;
  }
  else {
    for (f = tp->flake_list; f != NULL; f = f->next_flake) {
      if (f->next_flake == fp) {
	f->next_flake = fp->next_flake;
	break;
      }
    }
  }
  add_flake_to_reserve_list(fp);
  //free_flake (fp);
  tp->num_flakes--;
}

/* gives concentration-independent rates                                   */
/* for non-empty cells i,j, computes rate rv[n] to convert to type n       */
/*  so rv[0] gives the off-rate  (not the sum)                             */
/*  for chunk_fission, r=rv[0] gives the total off-rate over single tile,  */
/*     EW pair, NS pair, & 2x2 block, individually in [1+N+0]...[1+N+3]    */
/*  n=0 (empty) is the most common case.  the returned value is sum rv[n]  */
/* for empty cells i,j, returns 0                                          */
/* rv must already exist, of size at least fp->N+1 (+4 for chunks)         */
/* 0 <= i,j < 2^P                                                          */
double calc_rates(flake *fp, int i, int j, double *rv)
{
  int n,mi,ei,hi,mo,eo,ho; double r, sumr; tube *tp=fp->tube;
  unsigned char nN,nE,nS,nW; int N=fp->N; int size=(1<<fp->P);
  int seedchunk[4]; 

  if (rv!=NULL) for (n=0;n<=N+4;n++) rv[n]=0;
  if (tp==NULL) return 0;
  if (tp->T>0) return 0;   /* no off-rates: irreversible Tile Assembly Model */
  n = fp->Cell(i,j);
  if (n==0) return 0;                           /* no off-rate for empties   */
  if (tp->dt_left[n]) return 0;                 /* similarly, no off-rate for the  
						   right side of a double tile */
  // NOTE: w/o wander, seed site can't dissociate.  So set rate to zero.
  //       w/  wander, seed site can dissociate if there is a neighbor to
  //           move the seed to.  So set rate to zero if flake is monomer.
  // Similarly, for chunk_fission, we must zero the rates for seedchunks,
  // unless wander is on -- in which case we zero the rates if there are
  // no neighbors to move the seed to.  This is done below.
  seedchunk[0] = (i   == fp->seed_i && j   == fp->seed_j) || 
    (i == fp->seed_i && ((tp->dt_right[fp->seed_n] && j == fp->seed_j + 1)  ||
			 (tp->dt_left[fp->seed_n] && j == fp->seed_j - 1)));
  r = tp->k * exp(-Gse(fp,i,j,n)); 
  if (seedchunk[0] && (!wander || fp->tiles==1 || (fp->tiles==2 && fp->seed_is_double_tile)))
    r=0; 
  if (tp->dt_right[n]) {
    r *= exp(-Gse(fp,i,j+1,tp->dt_right[n])+2*fp->tube->Gse_EW[tp->dt_right[n]][n]);
    //printf("Setting right rate to %e.\n",r);
  }
  sumr=r;
  if (fission_allowed==2) {              // rates for pairs and 2x2 block dissoc
    if (rv!=NULL) rv[1+N+0]=r;
    seedchunk[1] = (i   == fp->seed_i && j+1 == fp->seed_j) || seedchunk[0];
    seedchunk[2] = (i+1 == fp->seed_i && j   == fp->seed_j) || seedchunk[0];
    seedchunk[3] = (i+1 == fp->seed_i && j+1 == fp->seed_j) || seedchunk[1] || seedchunk[2];
    if ( (seedchunk[1] && (!wander || fp->tiles==2)) || (!periodic && j+1==size) ) r=0; 
    else r = tp->k * exp(-chunk_Gse_EW(fp,i,j,n)) * (fp->Cell(i,j+1)!=0); 
    sumr+=r; if (rv!=NULL) rv[1+N+1]=r; 
    if ( (seedchunk[2] && (!wander || fp->tiles==2)) || (!periodic && i+1==size) ) r=0; 
    else r = tp->k * exp(-chunk_Gse_NS(fp,i,j,n)) * (fp->Cell(i+1,j)!=0); 
    sumr+=r; if (rv!=NULL) rv[1+N+2]=r; 
    if ( (seedchunk[3] && (!wander || fp->tiles==4)) || 
         (!periodic && (i+1==size || j+1==size)) ) r=0; 
    else r = tp->k * exp(-chunk_Gse_2x2(fp,i,j,n)) * 
	   (fp->Cell(i+1,j)!=0 && fp->Cell(i,j+1)!=0 && fp->Cell(i+1,j+1)!=0); 
    sumr+=r; if (rv!=NULL) rv[1+N+3]=r; 
  }
  if (rv!=NULL) rv[0]=sumr; 

  /*  hydrolysis model assumes tiles 1...N/2 non hydrolyzed, N/2+1...N hydro */
  if (tp->hydro) {
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
    nS=fp->Cell(i+1,j); nE=fp->Cell(i,j+1);
    ei = ((tp->tileb)[n][2] != 0 && (tp->tileb)[nS][0]==0) + 
      ((tp->tileb)[n][1] != 0 && (tp->tileb)[nE][3]==0);
    mi = ((tp->tileb)[n][2] != 
	  (tp->tileb)[nS][0] && (tp->tileb)[n][2]*(tp->tileb)[nS][0] !=0) +
      ((tp->tileb)[n][1] != 
       (tp->tileb)[nE][3] && (tp->tileb)[n][1]*(tp->tileb)[nE][3] !=0);
    hi = (nE>N/2) + (nS>N/2);

    /* empty and mismatched outputs */
    nN=fp->Cell(i-1,j); nW=fp->Cell(i,j-1);
    eo = ((tp->tileb)[n][0] != 0 && (tp->tileb)[nN][2]==0) + 
      ((tp->tileb)[n][3] != 0 && (tp->tileb)[nW][1]==0);
    mo = ((tp->tileb)[n][0] != 
	  (tp->tileb)[nN][2] && (tp->tileb)[n][0]*(tp->tileb)[nN][2] !=0) +
      ((tp->tileb)[n][3] != 
       (tp->tileb)[nW][1] && (tp->tileb)[n][3]*(tp->tileb)[nW][1] !=0);
    ho = (nN>N/2) + (nW>N/2);

    r = tp->k * (tp->kas + mi * tp->kam + ei * tp->kae + hi * tp->kah +
		 tp->kao * (mo * tp->kam + eo * tp->kae + ho * tp->kah));
    if (n<=N/2) {
      if (rv!=NULL) rv[n+N/2]=r; sumr += r;
    } else {
      r = r * exp(-tp->Gcb[n]+tp->Gcb[n-N/2]
                  -Gse(fp,i,j,n)+Gse(fp,i,j,n-N/2));
      if (rv!=NULL) rv[n-N/2]=r; sumr += r;
    }
  } 

  return sumr;
} // calc_rates()


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
/* Assumes that ii,jj can be anything -- boundary conditions are    */
/*   taken care of here.                                            */
void update_rates(flake *fp, int ii, int jj)
{
  int n,p; int size=(1<<fp->P); tube *tp=fp->tube;
  double oldrate, newrate; int oldempty, newempty;

  // wrap in case ii,jj go beyond the central field of 1-cell protection zone
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
	fp->rate[p][ii][jj] = fp->rate[p+1][ii<<1][jj<<1]   +
	  fp->rate[p+1][ii<<1][(jj<<1)+1] +
	  fp->rate[p+1][(ii<<1)+1][jj<<1] +
	  fp->rate[p+1][(ii<<1)+1][(jj<<1)+1];
      ii = (ii>>1); jj = (jj>>1);
    }
  }
  assert (fp->rate[0][0][0] >= 0);
} // update_rates()

void update_tube_rates(flake *fp)
{
  flake_tree *ftp=fp->tree_node; 
  double oldrate, newrate; int oldempty, newempty;
 
  if (ftp==NULL) return;
 
  oldrate = ftp->rate; oldempty=ftp->empty;
  newrate = fp->rate[0][0][0]; newempty=fp->empty[0][0][0];
  assert (newrate >= 0);
  while (ftp!=NULL) {
    ftp->rate+=newrate-oldrate;
    if (ftp->rate < -0.1) {
      printf("Old rate was %e, new rate is %e.\n",oldrate,ftp->rate);
      printf("Bad news -- negative rate.\n");
    }
    else {
      oldrate = oldrate +0;
    }
    ftp->empty+=newempty-oldempty;
    ftp=ftp->up;
  }
  
} // update_tube_rates()


int between_double_tile (flake *fp, tube *tp, int i, int j, unsigned char n) {
  if (n == 0) {
    return (tp->dt_right[fp->Cell(i,j-1)] || tp->dt_left[fp->Cell(i,j+1)]);
  }
  if (tp->dt_right[n]) {
    return (fp->Cell(i,j+1) != tp->dt_right[n]);
  }
  if (tp->dt_left[n]) {
    return (fp->Cell(i,j-1) != tp->dt_left[n]);
  }
  return 0;
}

/* convert Cell(i,j) to type n.                                       */
/* update all hierarchical rates and empty counts, in flake and tube. */
/* ALL modifications to the cell array go through this interface!     */
/* we only modify cells with 0 <= i,j < (1<<fp->P)                    */
/* but since it might be called out of range in FILL, we fix it up.   */
/* BUG: changes in concentration should change G for every tile, but  */
/* we don't update G automatically; also, if conc[n]==0, nan results. */
void change_cell(flake *fp, int i, int j, unsigned char n)
{
  int size=(1<<fp->P);  tube *tp=fp->tube; 
  if (periodic) { i=(i+size)%size; j=(j+size)%size; }
  else if (i<0 || i>=size || j<0 || j>=size) return; // can't change tiles beyond central field

  if (fp->Cell(i,j)==n) return;  // nothing to change!
  if (tp!=NULL) { /* flake has been added to a tube */
    //printf("Changing %d, %d from %d to %d.\n",i,j,fp->Cell(i,j),n);
    if (fp->Cell(i,j)==0) {                         /* tile addition */
      if (tp->conc[n]<=fp->flake_conc) return; // conc's can't go to zero!
      if (fp->tiles==1 && tp->conc[fp->seed_n]<=fp->flake_conc) return; // ditto
      fp->G += -log(tp->conc[n]) - Gse(fp,i,j,n);   
      tp->conc[n] -= fp->flake_conc; 
      tp->conc[0] -= fp->flake_conc;
      if (fp->tiles==1 || (fp->tiles==2 && fp->seed_is_double_tile)) { 
	// monomer flakes don't deplete []; now no longer monomer!
        tp->conc[fp->seed_n] -= fp->flake_conc; 
        tp->conc[0]          -= fp->flake_conc;
	if (tp->dt_right[fp->seed_n]) {
	  tp->conc[tp->dt_right[fp->seed_n]] -= fp->flake_conc;
	  tp->conc[0]          -= fp->flake_conc;
	}
	if (tp->dt_left[fp->seed_n]) {
	  tp->conc[tp->dt_left[fp->seed_n]] -= fp->flake_conc;
	  tp->conc[0]          -= fp->flake_conc;
	}
      }
      tp->stat_a++; fp->tiles++; 
      if (fp->tiles > tp->largest_flake_size) {
	tp->largest_flake = fp->flake_ID;
	tp->largest_flake_size = fp->tiles;
      }
      fp->mismatches += Mism(fp,i,j,n);
    } else if (n==0) {                              /* tile loss */
      tp->conc[0] += fp->flake_conc; 
      tp->conc[fp->Cell(i,j)] += fp->flake_conc; 
      fp->G += log(tp->conc[fp->Cell(i,j)]) + Gse(fp,i,j,fp->Cell(i,j));

      // monomer flakes don't deplete []; just became monomer!
      // zzz check this
      //if (fp->tiles==2 || (fp->tiles==3 && tp->dt_right[fp->Cell(i,j)])) { 
      if (fp->tiles==2 || (fp->tiles==3 && fp->seed_is_double_tile)) { 
        tp->conc[fp->seed_n] += fp->flake_conc; 
        tp->conc[0]          += fp->flake_conc;
	if (tp->dt_right[fp->seed_n]) {
	  tp->conc[tp->dt_right[fp->seed_n]] += fp->flake_conc;
	  tp->conc[0]          += fp->flake_conc;
	}
	if (tp->dt_left[fp->seed_n]) {
	  tp->conc[tp->dt_left[fp->seed_n]] += fp->flake_conc;
	  tp->conc[0]          += fp->flake_conc;
	}
      }
      tp->stat_d++; fp->tiles--; 
      fp->mismatches -= Mism(fp,i,j,fp->Cell(i,j));
    } else {                               /* tile hydrolysis or replacement */
      fp->G += Gse(fp,i,j,fp->Cell(i,j)) - Gse(fp,i,j,n) +
	log(tp->conc[fp->Cell(i,j)]) - log(tp->conc[n]) + 
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

  // If we've changed to a state we haven't seen before, and we're counting
  // unique visited states, record it.
  if (tp && tp->tracking_seen_states && !between_double_tile (fp,tp,i,j,n) &&
      !assembly_is_a_duplicate(tp->states_seen_hash,
			       fp->cell,size)) {
    add_assembly_to_seen(tp);
  }
  if (tp && tp->watching_states) {
    if (fp->chain_state) {
      update_state_off_indicator(fp);
    }  
    if (!between_double_tile (fp,tp,i,j,n)) {
      update_state_on_indicator(fp, fp->cell, size);
    }
  }
  // note: this recalculates all these rates from scratch, although we know only some can change
  update_rates(fp,i,j);
  update_rates(fp,i+1,j);
  update_rates(fp,i-1,j);
  update_rates(fp,i,j+1);
  update_rates(fp,i,j-1);
  // Also change these in case of double tiles or chunk_fission
  // TODO : just check if we have double tiles or chunk fission before
  // doing these tedious updates
  update_rates(fp,i-1,j+1);
  update_rates(fp,i+1,j+1);
  update_rates(fp,i-1,j-1);
  update_rates(fp,i+1,j-1);
  update_rates(fp,i,j+2);
  update_rates(fp,i-1,j+2);
  update_rates(fp,i+1,j+2);
  update_rates(fp,i,j-2);
  update_rates(fp,i-1,j-2);
  update_rates(fp,i+1,j-2);
  if (tp!=NULL) update_tube_rates(fp);
} // change_cell()

void change_seed(flake *fp, int new_i, int new_j)
{  int old_i=fp->seed_i; int old_j=fp->seed_j;
 fp->seed_n = fp->Cell(new_i,new_j); fp->seed_i=new_i; fp->seed_j=new_j;
 assert (!fp->tube->dt_left[fp->seed_n]);
 //printf("Seed is now tile %d at %d,%d.\n",fp->seed_n,fp->seed_i,fp->seed_j);
 // TODO: When double tiles are present, more rates should be updated.
 update_rates(fp, old_i, old_j);  // no longer being seed may allow dissoc => rates change
 update_rates(fp, old_i-1, old_j); 
 update_rates(fp, old_i+1, old_j); 
 update_rates(fp, old_i, old_j-1); 
 update_rates(fp, old_i, old_j+1); 
 update_rates(fp, new_i, new_j);  // now being seed may prevent dissoc => rates change
 update_rates(fp, new_i-1, new_j); 
 update_rates(fp, new_i+1, new_j); 
 update_rates(fp, new_i+1, new_j-1); 
 update_rates(fp, new_i-1, new_j-1); 
 // all this updating is painfully slow, since it must be done with every seed
 // motion during WANDER.
 if (fp->seed_n==0) 
   printf("seed wandered to empty site %d,%d!\n",new_i,new_j);
 if (fp->tube != NULL) 
   update_tube_rates(fp);
} // change_seed()


int choose_tile_type (tube *tp) {
  double r, cum;
  int oops, n;
  // This is exactly the same code as inside choose_cell,
  // but because that code is used so often, and this will be used for
  // adding a flake, which we imagine doing much less often, the
  // other code was left inline.
  r = drand48();
  do {
    r = r * tp->conc[0];  cum = 0;  oops=0;
    for (n=1; n<=tp->N; n++) if (r < (cum += tp->conc[n])) break; 
    if (n>tp->N) { // apparently conc[0] is not the sum of conc[n], oops
      printf("Concentration sum error!!! %f =!= %f\n",tp->conc[0],cum); 
      r=drand48(); oops=1; 
      tp->conc[0]=0; for (n=1; n <= tp->N; n++) tp->conc[0]+=tp->conc[n];
    }
  } while (oops);
  return n;
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
	printf("Concentration sum error!!! %f =!= %f\n",tp->conc[0],cum); 
	r=drand48(); oops=1; 
	tp->conc[0]=0; for (n=1; n <= tp->N; n++) tp->conc[0]+=tp->conc[n];
      }
    } while (oops);
  } else {                    /* choose off-event 0 or conversion to 1...N */
    if (tp->hydro) {
      sum = calc_rates(fp,i,j,tp->rv);
      if (sum==0) printf("Zero-sum hydro rate was chosen!!!\n");
      r = r * sum;  cum = 0;
      for (n=0; n<=fp->N; n++) if (r < (cum += tp->rv[n])) break; 
      if (n>fp->N) { printf("Hydro failed to choose anyone!!!\n"); n=0; }
    } else {
      n=0;  // always an off-event, unless hydrolysis rules are used.
    }
  }
  *np = n;
} // choose_cell()

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
    assert (ftp->left->rate >= -0.1);
    assert (ftp->right->rate >= -0.1);
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
} // choose_flake()

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


/* here, we lose the part that's unconnected to the seed, rather      */
/* than saving it and creating a new flake.                           */
/* the cell at ii, jj has already been removed (set to 0).            */
/* if fission_allowed==0, then the unconnected tiles are not removed. */
/* returns whether this dissociation does/would break the flake.      */
int flake_fission(flake *fp, int ii, int jj)
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
  if (i>0 || j>0) printf("FILL: %d groupies and %d wannabies on entry\n",i,j);
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

  /* all groups merged to one. clean up; we're fine. or, fission not allowed */
  if (ngroups==1 || fission_allowed==0) { 
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
  if (fission_allowed>0) tp->stat_f += ngroups-1;
  return (ngroups != 1); // would fission occur w/o this tile?
} // flake_fission()


/* Tile type n just dissociated from site i,j.  Was this safe to do?        */
/* If we can verify that what used to be connected still must be connected, */
/* then it is safe.                                                         */
int locally_fission_proof(flake *fp, int i, int j, int oldn)
{ unsigned char ringi; 
 ringi = ((fp->Cell(i-1,j)!=0)<<7) +
   ((CONNECTED_W(fp,i-1,j+1) && CONNECTED_S(fp,i-1,j+1))<<6) +
   ((fp->Cell(i,j+1)!=0)<<5) +
   ((CONNECTED_N(fp,i+1,j+1) && CONNECTED_W(fp,i+1,j+1))<<4) +
   ((fp->Cell(i+1,j)!=0)<<3) +
   ((CONNECTED_E(fp,i+1,j-1) && CONNECTED_N(fp,i+1,j-1))<<2) +
   ((fp->Cell(i,j-1)!=0)<<1) +
   ((CONNECTED_S(fp,i-1,j-1) && CONNECTED_E(fp,i-1,j-1))<<0);

 // safe if neighbors form one fully connected group, w/o central tile:
 if (ring[ringi]==1) return 1; 

 // next check, a little harder: label each neighbor, and merge connected groups.
 // if same labels regardless of whether central tile was used for merging,
 // then we're safe again.

 // This could be done a lot faster with a lookup table, but it's sufficiently
 // rare that I'm not bothering... for now...

 { int Nw=1,Nwo=1,Sw=2,Swo=2,Ew=3,Ewo=3,Ww=4,Wwo=4; 
 int Nc=(fp->Cell(i-1,j)!=0), Sc=(fp->Cell(i+1,j)!=0), 
   Ec=(fp->Cell(i,j+1)!=0), Wc=(fp->Cell(i,j-1)!=0), 
   NEc=(CONNECTED_W(fp,i-1,j+1) && CONNECTED_S(fp,i-1,j+1)),
   SEc=(CONNECTED_N(fp,i+1,j+1) && CONNECTED_W(fp,i+1,j+1)),
   SWc=(CONNECTED_E(fp,i+1,j-1) && CONNECTED_N(fp,i+1,j-1)),
   NWc=(CONNECTED_S(fp,i-1,j-1) && CONNECTED_E(fp,i-1,j-1));
 int Nh=HCONNECTED_N(fp,i,j,oldn), Eh=HCONNECTED_E(fp,i,j,oldn),
   Sh=HCONNECTED_S(fp,i,j,oldn), Wh=HCONNECTED_W(fp,i,j,oldn);
 int changed=1;

 while (changed) { changed=0;
 if (Nc && Ec && NEc && Nwo!=Ewo) {Nwo=Ewo=MIN(Nwo,Ewo); changed=1;}
 if (Nc && Wc && NWc && Nwo!=Wwo) {Nwo=Wwo=MIN(Nwo,Wwo); changed=1;}
 if (Sc && Ec && SEc && Swo!=Ewo) {Swo=Ewo=MIN(Swo,Ewo); changed=1;}
 if (Sc && Wc && SWc && Swo!=Wwo) {Swo=Wwo=MIN(Swo,Wwo); changed=1;}

 if (Nc && Ec && NEc && Nw!=Ew) {Nw=Ew=MIN(Nw,Ew); changed=1;}
 if (Nc && Wc && NWc && Nw!=Ww) {Nw=Ww=MIN(Nw,Ww); changed=1;}
 if (Sc && Ec && SEc && Sw!=Ew) {Sw=Ew=MIN(Sw,Ew); changed=1;}
 if (Sc && Wc && SWc && Sw!=Ww) {Sw=Ww=MIN(Sw,Ww); changed=1;}
 if (Nc && Ec && Nh && Eh && Nw!=Ew) {Nw=Ew=MIN(Nw,Ew); changed=1;}
 if (Nc && Wc && Nh && Wh && Nw!=Ww) {Nw=Ww=MIN(Nw,Ww); changed=1;}
 if (Sc && Ec && Sh && Eh && Sw!=Ew) {Sw=Ew=MIN(Sw,Ew); changed=1;}
 if (Sc && Wc && Sh && Wh && Sw!=Ww) {Sw=Ww=MIN(Sw,Ww); changed=1;}
 if (Sc && Nc && Sh && Nh && Sw!=Nw) {Sw=Nw=MIN(Sw,Nw); changed=1;}
 if (Ec && Wc && Eh && Wh && Ew!=Ww) {Ew=Ww=MIN(Ew,Ww); changed=1;}
 }

 if (Nw==Nwo && Sw==Swo && Ew==Ewo && Ww==Wwo) return 1;  
 }

 return 0;
}

/* remove all tiles whose off-rate more than is 'X' times faster than its on-rate. */
/* (these are calculated for individual tiles only; chunk_fission has no effect.) */
/* repeat 'iters' times. */
void clean_flake(flake *fp, double X, int iters)
{
  int i,j,n; double kc;  tube *tp=fp->tube;
  int size = (1<<fp->P); int it; int *F;

  F = (int *)calloc(size*size, sizeof(int));  /* scratch space */

  kc = tp->k*tp->conc[0]; /* on-rate */

  /* first memorize, then remove, to avoid changing rates during removal */
  for (it=0; it<iters; it++) {
    for (i=0; i<size; i++)
      for (j=0; j<size; j++) {
        n = fp->Cell(i,j);
        F[i+size*j] = (exp(-Gse(fp,i,j,n)) > X * tp->conc[n]);
      }
    for (i=0; i<size; i++)
      for (j=0; j<size; j++) 
        if (F[i+size*j]) {
	  n = fp->Cell(i,j); change_cell(fp, i,j,0);
	  if (!locally_fission_proof(fp,i,j,n)) /* couldn't quickly confirm... */
	    if (flake_fission(fp,i,j) && fission_allowed==0) {
	      change_cell(fp,i,j,n); tp->stat_a--; tp->stat_d--;
	    }
	}
  }
  free(F);
} // clean_flake()

/* add tiles whose off-rate less than 'X' times faster than its on-rate. */
/* (these are calculated for individual tiles only; chunk_fission has no effect.) */
/* repeat 'iters' times. */
void fill_flake(flake *fp, double X, int iters)
{
  int i,j,n; double secure, most_secure; double kc;  tube *tp=fp->tube;
  int size = (1<<fp->P); int it; int *F;

  F = (int *)calloc(size*size, sizeof(int));  /* scratch space */

  kc = tp->k*tp->conc[0]; /* on-rate */

  /* first memorize, then add, to avoid changing rates during removal */
  for (it=0; it<iters; it++) {
    for (i=0; i<size; i++)
      for (j=0; j<size; j++) {
        most_secure=0;
        if (fp->Cell(i,j)==0) for (n=1; n<=fp->N; n++) {
          secure = (exp(-Gse(fp,i,j,n)) - X * tp->conc[n]);
          if (secure<most_secure) { 
            F[i+size*j] = n; most_secure=secure;
	  }
	}
      }
    for (i=0; i<size; i++)
      for (j=0; j<size; j++) {
        if (F[i+size*j]) change_cell(fp, i, j, F[i+size*j]);
        F[i+size*j] = 0;
      }
  }
  free(F);
} // fill_flake()


/* Repair, as well as possible, a flake 
   [constructed by a tile set in which all errorless assemblies have no holes]
   * first, remove all tiles involved with mismatches. 
   * identify holes -- perform an inefficient fill from the edges of empty space
   -- unfilled, empty cells are "interior"
   * fill in interior sites where a unique strength-T tile may be added
   -- repeat until no longer possible
   * fill in interior sites with "the" tile that makes the strongest bond
   -- perform in rounds with decreasing minT, requiring at least minT*Gse to add 
   -- within each round, repeat until no longer possible

   NOTE: because this will always completely fill in any interior holes, it may construct
   a flake that is NOT connected by non-zero Gse -- thus, if simulation is continues,
   flake fission may not reach all relevant tiles when removing tiles in said interior hole, 
   leading to visibly non-connected regions.  This is a BUG.  But seldom relevant.
*/ 
void repair_flake(flake *fp, double T, double Gse)
{
  int i,j,n; tube *tp=fp->tube; int size = (1<<fp->P); 
  int morefill, bestn, numn; double bestGse; int *F; 
  double minT,bigT; int n1,n2;

  F = (int *)calloc(size*size, sizeof(int));  /* scratch space */

  // first, remove all tiles involved with mismatches -- identify, then remove
  for (i=0; i<size; i++)
    for (j=0; j<size; j++)  
      if ((n=fp->Cell(i,j))>0 && Mism(fp,i,j,n)) F[i+size*j] = 1;
  for (i=0; i<size; i++)
    for (j=0; j<size; j++) 
      if (F[i+size*j]) {
	n = fp->Cell(i,j); change_cell(fp, i,j,0); F[i+size*j]=0;
	if (!locally_fission_proof(fp,i,j,n)) /* couldn't quickly confirm... */
	  if (flake_fission(fp,i,j) && fission_allowed==0) {
	    change_cell(fp,i,j,n); tp->stat_a--; tp->stat_d--;
	  }
      }

  // identify holes -- perform an inefficient fill from the edges of empty space 
  //                -- unfilled, empty cells are "interior"                       
  for (i=0; i<size; i++) {
    j=0;      while (j<size && fp->Cell(i,j)==0) { F[i+size*j] = 1; j++; }
    j=size-1; while (j>=0   && fp->Cell(i,j)==0) { F[i+size*j] = 1; j--; }
    j=0;      while (j<size && fp->Cell(j,i)==0) { F[j+size*i] = 1; j++; }
    j=size-1; while (j>=0   && fp->Cell(j,i)==0) { F[j+size*i] = 1; j--; }
  }
  do { morefill=0;
  for (i=1; i<size-1; i++)
    for (j=1; j<size-1; j++)  
      if (F[i+size*j]==0 && fp->Cell(i,j)==0)  
	if (F[i+1+size*j]==1 || F[i+size*j+size]==1 || 
	    F[i-1+size*j]==1 || F[i+size*j-size]==1) 
	  { F[i+size*j] = 1; morefill=1; }
  } while (morefill);
  // now all "exterior" cells have Fgroup set to 1; tiles and interior empties are 0

  // fill in interior sites where a unique strength-T tile may be added
  //                -- repeat until no longer possible
  do { morefill=0;
  for (i=0; i<size; i++)
    for (j=0; j<size; j++)  
      if (fp->Cell(i,j)==0 && F[i+size*j]==0) {
	bestn=0; numn=0;
	for (n=1; n<=fp->N; n++) if (Gse(fp,i,j,n)>=T*Gse) { bestn=n; numn++; }

	// we can change immediately, since if there is a correct fill-in, 
	// it doesn't matter what order the fill-in occurs by unique steps

	if (numn==1) { change_cell(fp,i,j,bestn); morefill=1; }
      }
  } while (morefill);

  // fill in interior sites with "the" tile that makes the strongest bond
  //                -- perform in rounds with decreasing minT, requiring at least minT*Gse to add 
  //                -- within each round, repeat until no longer possible

  bestGse=0;
  for (n1=1; n1<=tp->N; n1++) 
    for (n2=1; n2<=tp->N; n2++) 
      bestGse=MAX(bestGse,MAX(tp->Gse_NS[n1][n2],tp->Gse_EW[n1][n2]));
  bigT=ceil(4*bestGse/Gse); // over-estimate max bond-strength to hold in a single tile

  for (minT=bigT; minT>-1; minT-=1.0)
    do { morefill=0;
    for (i=0; i<size; i++)
      for (j=0; j<size; j++)  
	if (fp->Cell(i,j)==0 && F[i+size*j]==0) {
	  bestn=0; bestGse=minT*Gse;
	  for (n=1; n<=fp->N; n++) if (Gse(fp,i,j,n)>=bestGse) 
	    { bestn=n; bestGse=Gse(fp,i,j,n); }

	  // we can change immediately, since at this point we're sure to make mistakes anyway

	  if (bestn>0) { change_cell(fp,i,j,bestn); morefill=1; }
	}
    } while (morefill);

  free(F);
} // repair_flake()


/* recalculate # mismatches counting only tiles w/o an empty space within rad.       */
/* (also see recalc_G for original #mismatches count, as displayed in window always) */
void error_radius_flake(flake *fp, double rad)
{
  int n,i,j,ii,jj,size=(1<<fp->P),solid; 

  fp->mismatches=0; 
  for (i=0;i<size;i++)
    for(j=0;j<size;j++) {
      if ((n=fp->Cell(i,j))>0 && Mism(fp,i,j,n)) {
	solid=1;
	for (ii=MAX(0,floor(i-rad)); ii<=MIN(size-1,ceil(i+rad)); ii++)
	  for (jj=MAX(0,floor(j-rad)); jj<=MIN(size-1,ceil(j+rad)); jj++)
	    if ( (ii-i)*(ii-i)+(jj-j)*(jj-j) < rad*rad && fp->Cell(ii,jj)==0 ) solid=0;
	fp->mismatches += solid;
      }
    }
  fp->mismatches/=2;  // errors right on boundary of radius may not be counted twice.

} // error_radius_flake()


int double_tile_allowed(tube *tp, flake *fp, int i, int j, int n) {
  int size = (1<<(tp->P));
  int right_side_can_attach, left_side_can_attach;
  int j_norm;
  if (!tp->dt_right[n] && !tp->dt_left[n]) 
    return 1;
  j_norm = (j+size)%size;
  left_side_can_attach =  tp->dt_right[n] && ((periodic || j + 1 < size) && fp->Cell(i,j_norm+1) == 0 && 
					      !HCONNECTED(fp,i,j_norm+1,tp->dt_right[n]));
  //if (tp->dt_right[n] && !left_side_can_attach) 
  //  printf("Rejecting %d,%d for left side tile %d.\n",i,j,n);
  right_side_can_attach =  tp->dt_left[n] && ((periodic || j - 1 >= 0) && fp->Cell(i,j_norm-1) == 0);
  //if (tp->dt_left[n] && !right_side_can_attach) 
  //  printf("Rejecting %d,%d for right side tile %d.\n",i,j,n);
  return (left_side_can_attach || right_side_can_attach);

}

int not_in_block(int i,int j,const int *di,const int *dj,int t,int n, int size) {
  int x, outside=1;
  for (x = 0; x < n; x++) {
    if (x == t) { continue; }
    if ((di[x] == i && dj[x] == j) || 
	(periodic && ((di[x] + size)%size) == i && ((dj[x] + size)%size == j))) { 
      outside = 0;
      break;
    }
  }
  return outside;
}

int adjacent_tiles(int i_one, int j_one, int i_two, int j_two, int size) {
  return ((i_one == i_two && j_one == j_two + 1) ||
	  (i_one == i_two && j_one == j_two - 1) ||
	  (i_one == i_two + 1 && j_one == j_two) ||
	  (i_one == i_two - 1 && j_one == j_two) ||
	  (periodic && i_one == i_two && 
	   ((j_one+1+size)%size == (j_two+size)%size ||
	    (j_one-1+size)%size == (j_two+size)%size)));
}

void order_removals(tube *tp, flake *fp,
		    int n, const int *di,const int *dj,int removals[]) {
  int index,t;
  int already_chosen, outside_neighbors, adjacent_to_already_chosen;
  int x;
  int size;
  // Remove tiles one by one in an order that satisfies the following
  // properties (before any removals occur)
  // 1.  The last tile to be removed is bound to a tile not to be removed.
  // 2.  Each preceding tile is either bound to a tile to be removed after
  // it or satisfies (1).
  // First, we must leave at least one tile:
  assert (n < fp->tiles);
  index = n-1;
  size = (1<<tp->P);
  // index = place in removal order
  // t = place in unordered chunk
  for (index = n-1; index >= 0; index--) {
    for (t = 0; t < n; t++) {
      already_chosen = 0;
      for (x = n-1; x > index; x--) {
	if (removals[x] == t) {
	  already_chosen = 1;
	  break;
	}
      }
      if (!already_chosen) {
	int di_up, di_down, dj_left, dj_right;
	if (periodic) {
	  di_up = (di[t]-1+size)%size;
	  di_down = (di[t]+1+size)%size;
	  dj_left = (dj[t]-1+size)%size;
	  dj_right = (dj[t]+1+size)%size;
	}
	else {
	  di_up = di[t]-1;
	  di_down = di[t]+1;
	  dj_left = dj[t] - 1;
	  dj_right = dj[t] + 1;
	}
	outside_neighbors = 0;

	// Check bottom
	if (fp->Cell(di_down,dj[t]) && 
	    not_in_block(di_down,dj[t],di,dj,t,n,size)) {
	  outside_neighbors = 1;
	}
	// Top 
	if (fp->Cell(di_up,dj[t]) && 
	    not_in_block(di_up,dj[t],di,dj,t,n,size)) {
	  outside_neighbors = 1;
	}
	// Left
	if (fp->Cell(di[t],dj_left) && 
	    not_in_block(di[t],dj_left,di,dj,t,n,size)) {
	  outside_neighbors = 1;
	}
	// Right
	if (fp->Cell(di[t],dj_right) && 
	    not_in_block(di[t],dj_right,di,dj,t,n,size)) {
	  outside_neighbors = 1;
	}
	if (outside_neighbors) {
	  removals[index] = t;
	  break;
	}
	else {
	  adjacent_to_already_chosen = 0;
	  size = (1<<tp->P);
	  for (x = n-1; x > index; x--) {
	    if (adjacent_tiles(di[removals[x]],dj[removals[x]],di[t],dj[t],size)) {
	      adjacent_to_already_chosen = 1;
	    }
	  }
	  if (adjacent_to_already_chosen) {
	    removals[index] = t;
	    break;
	  }   
	}
      }
    }
    if (t ==n) {
     printf("arg!\n");
    }
    assert(t < n);
  }
}

void update_all_rates(tube *tp) {
  flake *fp;
  for (fp = tp->flake_list; fp != NULL; fp=fp->next_flake) {
    recalc_G(fp);
  }
}

void get_random_wander_permutation (int di[6], int dj[6], 
				    int seed_is_double_tile) {
  int perm,x;
  int perm_nums[6];
  if (seed_is_double_tile) {
    int taken[6] = {0,0,0,0,0,0};
    int posx;
    perm=random()%720;
    for (x = 0; x < 6; x++) {
      posx = perm%(6-x);
      perm_nums[x]=0;
      while(taken[perm_nums[x]] || posx) {
	if (!taken[perm_nums[x]]) {
	  posx--;
	}
	perm_nums[x]++;
      }
      taken[perm_nums[x]] = 1;
      perm = perm/(6-x);
    }
    for (x = 0; x < 6; x++) {
      switch (perm_nums[x]) {
      case 0:
	di[x] = -1;
	dj[x] = 0;
	break;
      case 1:
	di[x] = -1;
	dj[x] = 1;
	break;
      case 2:
	di[x] = 0;
	dj[x] = 2;
	break;
      case 3:
	di[x] = 1;
	dj[x] = 1;
	break;
      case 4:
	di[x] = 1;
	dj[x] = 0;
	break;
      case 5:
	di[x] = 0;
	dj[x] = -1;
	break;
      default:
	assert(0);
      }
    }
  }
  else {
    int taken[4] = {0,0,0,0}, posx;
    perm = random()%24;
    for (x = 0; x < 4; x++) {
      posx = perm%(4-x);
      perm_nums[x]=0;
      while(taken[perm_nums[x]] || posx) {
	if (!taken[perm_nums[x]]) {
	  posx--;
	}
	perm_nums[x]++;
      }
      taken[perm_nums[x]] = 1;
      perm = perm/(4-x);
    }
    for (x = 0; x < 4; x++) {
      switch (perm_nums[x]) {
      case 0:
	di[x] = -1;
	dj[x] = 0;
	break;
      case 1:
	di[x] = 0;
	dj[x] = -1;
	break;
      case 2:
	di[x] = 1;
	dj[x] = 0;
	break;
      case 3:
	di[x] = 0;
	dj[x] = 1;
	break;
      default :
	assert(0);
      }
    }
  }
}

/* simulates 'events' events */
void simulate(tube *tp, int events, double tmax, int emax, int smax, int fsmax)
{
  int i,j,n,oldn; double dt; flake *fp; int chunk, seedchunk[4];
  double total_rate, total_blast_rate, new_flake_rate, event_choice; long int emaxL;
  int size=(1<<tp->P), N=tp->N;  

  if (tp->flake_list==NULL && tp->tinybox == 0) return;  /* no flakes! */
  //   if (tp->events + 2*events > INT_MAX) {
  if (tp->events + 2*events > 1000000000) {
    tp->ewrapped=1; tp->events=0; 
    tp->stat_a-=tp->stat_d; tp->stat_d=0; tp->stat_h=0; tp->stat_f=0;
    fp=tp->flake_list; 
    while (fp!=NULL) {  
      fp->events=0;
      fp=fp->next_flake;
    }
  }

  emaxL = (emax==0 || tp->events+events<emax)?(tp->events+events):emax;


  fp=tp->flake_list; 
  while (fp!=NULL) {  
     assert (!tp->tinybox ||
	     ((!fp->seed_is_double_tile && fp->tiles > 1) || fp->tiles > 2));
     
    //printf("Seed is tile %d at %d,%d.\n",fp->seed_n,fp->seed_i,fp->seed_j);
    if (periodic) {
      assert(fp->Cell((fp->seed_i+size)%size,(fp->seed_j+size)%size) == fp->seed_n);
    } else {
      assert(fp->Cell(fp->seed_i,fp->seed_j) == fp->seed_n);
    }


    if (tp->dt_right[fp->seed_n]) {
      assert(fp->Cell(fp->seed_i,fp->seed_j+1) == tp->dt_right[fp->seed_n]);
    }
    if (tp->dt_left[fp->seed_n]) {
      assert(fp->Cell(fp->seed_i,fp->seed_j-1) == tp->dt_left[fp->seed_n]);
    }
    fp=fp->next_flake;
  }
  for (i = 0; i < N; i++) {
    if (tp->dt_right[i]) {
      if (tp->conc[i] != tp->conc[tp->dt_right[i]]) {
	printf("Concentrations are off!\n");
      }
      assert (tp->conc[i] == tp->conc[tp->dt_right[i]]);
    }
  }

  if (tp->flake_tree) {
    total_rate = tp->flake_tree->rate+tp->k*tp->conc[0]*tp->flake_tree->empty;
  }
  else {
    total_rate = 0;
  }
  total_blast_rate = tp->k*tp->conc[0]*blast_rate*size*size*tp->num_flakes;
  new_flake_rate = tp->k*2*tp->conc[0]*tp->conc[0]*tp->tinybox*AVOGADROS_NUMBER ;
  
  assert (total_rate + total_blast_rate + new_flake_rate >= 0); // can be zero in aTAM if finite-sized assembly is done
  
  while (tp->events < emaxL && 
	 (tmax==0 || tp->t < tmax) && 
          (smax==0 || tp->stat_a-tp->stat_d < smax) &&
	 (fsmax==0 || tp->largest_flake_size < fsmax) &&
	 total_blast_rate+total_rate+new_flake_rate > 0) {

    /* First check if time is such that we need to update the temperature */
    if (tp->anneal_t && (tp->t > tp->next_update_t)) {
      tp->Gse = tp->Gse_final- (tp->Gse_final - tp->anneal_g)*exp(-tp->t/tp->anneal_t);
      set_Gses(tp,tp->Gse,0);  // NOT SAFE FOR HYDROLYSIS
      /* Now we have to update all rates */
      tp->updates++;
      tp->next_update_t = tp->updates*tp->anneal_t*log(2)/tp->update_freq;
      update_all_rates (tp);
    }
    dt = -log(drand48()) / (total_rate + total_blast_rate + new_flake_rate);
    event_choice = drand48()*(total_rate+total_blast_rate+new_flake_rate);
    
    if (blast_rate>0 && event_choice < total_blast_rate) {  
      int kb=size,ii,jj,ic,jc,di,dj,seed_here,flake_n;
      
      while(kb==size) { double dr = drand48()*blast_rate;
         for (kb=1; kb<size; kb++)  // choose blast hole size kb= 1...size
           if ( ( dr -= blast_rate_alpha * exp(-blast_rate_gamma*(kb-1)) / pow(kb*1.0,blast_rate_beta) ) < 0 )
	     break; 
      }
      // printf("zap! %d x %d\n",kb,kb);

       // choose a flake
       flake_n = random()%(tp->num_flakes); fp=tp->flake_list;  for (i=0; i<flake_n; i++) fp=fp->next_flake;

       ic=random()%size; jc=random()%size;  // corner coordinates for kb x kb square to be removed
       di=2*(random()%2)-1; dj=2*(random()%2)-1;  // square goes in random direction from ic, jc

       for (seed_here=0, ii=0; ii<kb; ii++) for (jj=0; jj<kb; jj++) { // make sure seed tile is not in square
         if (periodic) { i=(ic+di*ii+size)%size; j=(jc+dj*jj+size)%size; } else { i=ic+di*ii; j=jc+dj*jj; }
         if (i==fp->seed_i && j==fp->seed_j) seed_here=1;  // square wraps or is cropped
       }
       if (!seed_here) { int vorh=random()%2;
         for (ii=0; ii<kb; ii++) for (jj=0; jj<kb; jj++) {
           if (vorh) { if (periodic) { i=(ic+di*ii+size)%size; j=(jc+dj*jj+size)%size; } else { i=ic+di*ii; j=jc+dj*jj; } }
           else      { if (periodic) { i=(ic+di*jj+size)%size; j=(jc+dj*ii+size)%size; } else { i=ic+di*jj; j=jc+dj*ii; } }
           if (i>=0 && j>=0 && i<size && j<size && fp->Cell(i,j)>0) {  
             // might have been removed already by previous fission or was never there; or maybe i j needs to be cropped
             oldn = fp->Cell(i,j); change_cell(fp,i,j,0);  // now it's gone!
             if (!locally_fission_proof(fp,i,j,oldn)) // may remove additional tiles (not seed)
	       if (flake_fission(fp,i,j) && fission_allowed==0) {  // see below under "dissociation" for comments
		change_cell(fp,i,j,oldn); tp->stat_a--; tp->stat_d--;
                ii=kb; jj=kb;  // stop the blast (without this, it also works, but looks weird)
	      }
           }
         }
       }
    } else if (new_flake_rate && event_choice < (total_blast_rate + new_flake_rate)) {
      int m,r,x,d,di,dj,c;
      // Add a new flake
      // Select the seed tile for the new flake
      n = choose_tile_type (tp);

      m = choose_tile_type (tp);
      // Choose a second cell to add to the new tile to
      // Determine an orientation for the two tiles
      r = random();

      x = ((r>>2) % 2) * 2 - 1;
      d = r % 2;
      // Determine whether the two are connected
      if (d) {
	// Connect up and down
	dj = 0;
	di = x;
	if (x > 0) {
	  c = tp->Gse_NS[n][m];
	} else { c = tp->Gse_NS[m][n]; }
      }
      else {
	// Connect left or right
	dj = x;
	di = 0;
	if ((tp->dt_left[n] && dj < 0) ||
	    (tp->dt_right[n] && dj > 0)) {
	  c = 0;
	}
	else {
	  if (x > 0) { 
	    c = tp->Gse_EW[m][n]; 
	  } else { c = tp->Gse_EW[n][m]; }
	}
      }
      if (c) {
	int s_n, s_j;
	//printf("Initting flake with tile %d and tile %d at %d,%d and %d,%d.\n",
	//     n,m,tp->default_seed_i,tp->default_seed_j,tp->default_seed_i+di,tp->default_seed_j+dj); 
	if (tp->dt_left[n]) {
	  s_n = tp->dt_left[n];
	  s_j = tp->default_seed_j - 1;
	}
	else {
	  s_n = n;
	  s_j = tp->default_seed_j;
	}
	if ((fp = recover_flake (tp->default_seed_i,s_j,s_n,tp->initial_Gfc)) == NULL) {
	  fp = init_flake (tp->P,tp->N,tp->default_seed_i, s_j, s_n, tp->initial_Gfc);
	}
	insert_flake (fp, tp);
	fp->tiles = 1;
	fp->seed_is_double_tile = tp->dt_right[fp->seed_n];
	if (tp->dt_right[s_n]) {
	  change_cell(fp, tp->default_seed_i, s_j+1,tp->dt_right[s_n]);
	}
	change_cell(fp,tp->default_seed_i+di,tp->default_seed_j+dj,m);
	if (tp->dt_right[m]) {
	  change_cell(fp, tp->default_seed_i+di, tp->default_seed_j+dj+1,tp->dt_right[m]);
	}
	if (tp->dt_left[m]) {
	  change_cell(fp, tp->default_seed_i+di, tp->default_seed_j+dj-1,tp->dt_left[m]);
	}
	assert (!tp->tinybox ||
	       ((!fp->seed_is_double_tile && fp->tiles > 1) || fp->tiles > 2));
      }
      else {
	tp->stat_a++; tp->stat_d++; tp->events+=2; 
      }       
    }

    else { // blast error case and new flake case above, kTAM / aTAM below
	
     fp=choose_flake(tp);
     assert (!tp->dt_left[fp->Cell(fp->seed_i,fp->seed_j)]);
     assert (!tp->tinybox ||
	     fp->tiles > 2 || (!fp->seed_is_double_tile && fp->tiles > 1));

     /* let the designated seed site wander around */
     /* must do this very frequently, else treadmilling would get stuck */
     if (wander) {  
      int new_i, new_j;
      // Pick a new seed adjacent to the old one
      new_i = fp->seed_i-1+random()%3;
      new_j = fp->seed_j-1+random()%3;
      if (periodic  || (new_i>=0 && new_i<size && new_j>=0 && new_j<size)) {
	//printf("Size is %d, new_i is %d, new_j is %d.\n",size,new_i,new_j);
	if (periodic) { new_i=(new_i+size)%size; new_j=(new_j+size)%size; }
	//printf("After renormalize, size is %d, new_i is %d, new_j is %d.\n",size,new_i,new_j);
	if (fp->Cell(new_i,new_j) != 0 && (new_i != fp->seed_i || new_j != fp->seed_j) && 
	    !tp->dt_left[fp->Cell(new_i,new_j)]) { 
	  change_seed(fp,new_i,new_j);
	  fp->seed_is_double_tile = tp->dt_right[fp->seed_n];
	  //assert((!fp->seed_is_double_tile && fp->tiles > 1) || fp->tiles > 2);
	}
      }
      // pick a new seed if 0 of the old one are left, and we are a single tile
      if ((fp->tiles==1 || (fp->tiles==2 && fp->seed_is_double_tile)) && 
	  tp->conc[fp->seed_n]<=fp->flake_conc) {
	int e,i;
	// Assuming there is enough of any free tile left:
	e = 0;
	for (i = 1; i < N; i++) {
	  if (tp->conc[i] > fp->flake_conc) { e=1; }
	}
	if (e) {
	  // Find a random new tile
	  while (tp->conc[fp->seed_n] < fp->flake_conc || tp->dt_left[fp->seed_n]) {
	    fp->seed_n=(random()%N)+1; 
	  }
	  change_cell(fp,fp->seed_i,fp->seed_j,0);
	  // In case the old seed was a double tile
	  change_cell(fp,fp->seed_i,fp->seed_j+1,0);
	  change_cell(fp,fp->seed_i,fp->seed_j,fp->seed_n); 
	  if (tp->dt_right[fp->seed_n]) {
	    change_cell(fp,fp->seed_i,fp->seed_j+1,tp->dt_right[fp->seed_n]); 
	  }
	}
      }
    }

    choose_cell(fp, &i, &j, &n); chunk=0;
    //printf("Chose cell %d,%d tile %d.\n",i,j,n);
    if (fission_allowed==2 && n==0) { // must choose either single tile, EW/NS pairs, or block
      double sum=0, rsum; 
      sum = calc_rates(fp,i,j,tp->rv); 
      rsum=sum*drand48();
      if (sum == 0) {
	// In this case, which chunk to choose is arbitrary, but we'll select 0
	chunk = 0;
      }
      else {
	for (chunk=0; chunk<4; chunk++) 
	  if (rsum<tp->rv[1+N+chunk]) break; else rsum-=tp->rv[1+N+chunk];
      }
      if (chunk<0 || chunk>3) printf("impossible chunk chosen!!!\n");
    }
    seedchunk[0] = ((i   == fp->seed_i && j   == fp->seed_j) ||
		    (tp->dt_right[fp->seed_n] && i == fp->seed_i && j == fp->seed_j + 1) ||
		    (tp->dt_left[fp->seed_n] && i == fp->seed_i && j == fp->seed_j - 1));
    seedchunk[1] = (n==0 && i   == fp->seed_i && j+1 == fp->seed_j) || seedchunk[0];
    seedchunk[2] = (n==0 && i+1 == fp->seed_i && j   == fp->seed_j) || seedchunk[0];
    seedchunk[3] = (n==0 && i+1 == fp->seed_i && j+1 == fp->seed_j) || seedchunk[1] || seedchunk[2];
    if (wander && n==0) {
      int new_i=fp->seed_i, new_j=fp->seed_j;
      if (chunk==0 && seedchunk[0]) {
	// looks like we're trying to dissociate the seed tile.
	// if 'wander' is set, this can happen -- if there is a neighboring
	// tile to move the seed to.
	int mi[6],mj[6], x, limit;
	assert (!tp->dt_left[fp->seed_n]);
	get_random_wander_permutation (mi, mj, fp->seed_is_double_tile);
	if (fp->seed_is_double_tile) {
	  limit = 6;
	}
	else { limit = 4; }
	for (x = 0; x < limit; x++) {
	  if ((fp->Cell(fp->seed_i+mi[x],fp->seed_j+mj[x]) &&
	      !tp->dt_left[fp->Cell(fp->seed_i+mi[x],fp->seed_j+mj[x])]) ||
	      (periodic && fp->Cell((fp->seed_i+mi[x]+size)%size,
				    (fp->seed_j+mj[x]+size)%size) &&
	      !tp->dt_left[fp->Cell((fp->seed_i+mi[x]+size)%size,
				(fp->seed_j+mj[x]+size)%size)])) {
	    new_i = new_i + mi[x];
	    new_j = new_j + mj[x];
	    break;
	  }
	}
	if ((!fp->seed_is_double_tile && x == 4) || x == 6) {
	  // Only place to move to was the right side of a double tile.
	  // Repeat and go one to the left to avoid the double tile
	  for (x = 0; x < limit; x++) {
	    if (fp->Cell(fp->seed_i+mi[x],fp->seed_j+mj[x]) ||
		(periodic && fp->Cell((fp->seed_i+mi[x]+size)%size,
				      (fp->seed_j+mj[x]+size)%size))) {
	      assert (tp->dt_left[fp->Cell(fp->seed_i+mi[x],fp->seed_j+mj[x])]);
	      new_i = new_i + mi[x];
	      new_j = new_j + mj[x]-1;
	      break;
	    }
	  }
	}
	if (periodic) { new_i = (new_i+size)%size; new_j = (new_j+size)%size; }
      } else if (chunk==1 && seedchunk[1]) {
	int mi,mj,mk; mi=((random()/17)%2)*2-1; mj=((random()/17)%2)*3-1; mk=(random()/17)%2;
	if      (fp->Cell(i-mi,j+mk)!=0)   { new_i=i-mi; new_j=j+mk; }
	else if (fp->Cell(i+mi,j+mk)!=0)   { new_i=i+mi; new_j=j+mk; }
	else if (fp->Cell(i-mi,j+1-mk)!=0) { new_i=i-mi; new_j=j+1-mk; }
	else if (fp->Cell(i+mi,j+1-mk)!=0) { new_i=i+mi; new_j=j+1-mk; }
	else if (fp->CellM(i,j-mj+1)!=0)   { new_i=i;    new_j=j-mj+1; }
	else if (fp->CellM(i,j+mj)!=0)     { new_i=i;    new_j=j+mj; }
      } else if (chunk==2 && seedchunk[2]) {
	int mi,mj,mk; mi=((random()/17)%2)*3-1; mj=((random()/17)%2)*2-1; mk=(random()/17)%2;
	if      (fp->Cell(i+mk,j+mj)!=0)   { new_i=i+mk;   new_j=j+mj; }
	else if (fp->Cell(i+mk,j-mj)!=0)   { new_i=i+mk;   new_j=j-mj; }
	else if (fp->Cell(i+1-mk,j+mj)!=0) { new_i=i+1-mk; new_j=j+mj; }
	else if (fp->Cell(i+1-mk,j-mj)!=0) { new_i=i+1-mk; new_j=j-mj; }
	else if (fp->CellM(i-mi+1,j)!=0)   { new_i=i-mi+1; new_j=j; }
	else if (fp->CellM(i+mi,j)!=0)     { new_i=i+mi;   new_j=j; }
      } else if (chunk==3 && seedchunk[3]) {
	int mi,mj,mk; mi=((random()/17)%2)*3-1; mj=((random()/17)%2)*3-1; mk=(random()/17)%2;
	if      (fp->CellM(i-mi+1,j+mk)!=0)   { new_i=i-mi+1; new_j=j+mk; }
	else if (fp->CellM(i+mi,j+mk)!=0)     { new_i=i+mi; new_j=j+mk; }
	else if (fp->CellM(i-mi+1,j+1-mk)!=0) { new_i=i-mi+1; new_j=j+1-mk; }
	else if (fp->CellM(i+mi,j+1-mk)!=0)   { new_i=i+mi; new_j=j+1-mk; }
	else if (fp->CellM(i+mk,j+mj)!=0)     { new_i=i+mk; new_j=j+mj; }
	else if (fp->CellM(i+mk,j-mj+1)!=0)   { new_i=i+mk; new_j=j-mj+1; }
	else if (fp->CellM(i+1-mk,j+mj)!=0)   { new_i=i+1-mk; new_j=j+mj; }
	else if (fp->CellM(i+1-mk,j-mj+1)!=0) { new_i=i+1-mk; new_j=j-mj+1; }
      }
      if (new_i != fp->seed_i || new_j != fp->seed_j) { 
	change_seed(fp,new_i,new_j);
	assert (!tp->dt_left[fp->seed_n]);
	fp->seed_is_double_tile = tp->dt_right[fp->seed_n];
      }
      else {
	assert(!seedchunk[0] ||
	       fp->tiles==1 || (fp->seed_is_double_tile && fp->tiles==2));
      }
      seedchunk[0] = ((i   == fp->seed_i && j   == fp->seed_j) ||
		      (tp->dt_right[fp->seed_n] && i == fp->seed_i && j == fp->seed_j + 1) ||
		      (tp->dt_left[fp->seed_n] && i == fp->seed_i && j == fp->seed_j - 1) ||
		      (periodic && i == fp->seed_i && 
		       ((tp->dt_right[fp->seed_n] && j == (fp->seed_j+1+size)%size) ||
		       (tp->dt_left[fp->seed_n] && j == (fp->seed_j-1+size)%size))));
      seedchunk[1] = (i   == fp->seed_i && j+1 == fp->seed_j) || seedchunk[0];
      seedchunk[2] = (i+1 == fp->seed_i && j   == fp->seed_j) || seedchunk[0];
      seedchunk[3] = (i+1 == fp->seed_i && j+1 == fp->seed_j) || seedchunk[1] || seedchunk[2];
    }
    oldn = fp->Cell(i,j);
    tp->t += dt;  // increment time whether or not seed tile event is effective
    //     if (seedchunk[chunk]) 
    //       printf("seedchunk triggered by %d,%d chunk %d !\n",i,j,chunk);
    if (!seedchunk[chunk]) { /* can't change seed tile */
      if (tp->T>0) { /* irreversible Tile Assembly Model */
        if (n>0 && Gse(fp,i,j,n)>=tp->T) { 
	  change_cell(fp,i,j,n);
	  if (tp->dt_right[n])
	      change_cell(fp,i,j+1,tp->dt_right[n]);
	  if (tp->dt_left[n])
	    change_cell(fp,i,j-1,tp->dt_left[n]);
	}
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
      } else { /* reversible kTAM model */
	if (zero_bonds_allowed==0) {
	  /* to add, must have se contact */
	  if (oldn==0 && HCONNECTED(fp,i,j,n) && 
	      double_tile_allowed(tp,fp,i,j,n)) {
	    change_cell(fp,i,j,n);
	    if (tp->dt_right[n])
	      change_cell(fp,i,j+1,tp->dt_right[n]);
	    if (tp->dt_left[n])
	      change_cell(fp,i,j-1,tp->dt_left[n]);
	  }
	  /* zero-bond tile additions fall off immediately;
	     count them or else, if there are lots, the display
	     can be super-slow! */
	  if (oldn==0 && !HCONNECTED(fp,i,j,n)) 
	    { tp->stat_a++; tp->stat_d++; tp->events+=2; fp->events+=2; } 
	}
	else if (oldn==0 && double_tile_allowed(tp,fp,i,j,n)) {
	  change_cell(fp,i,j,n);
	  if (tp->dt_right[n])
	    change_cell(fp,i,j+1,n);
	  if (tp->dt_left[n])
	    change_cell(fp,i,j-1,n);
	}
	
	/* hydrolysis happens here */
	if (oldn>0 && n>0) change_cell(fp,i,j,n);
	if (n==0) {  /* dissociation: check connectedness */
          int d, k, dn, di[4], dj[4], oldns[4], removals[4];
          if      (chunk==0 && !tp->dt_right[oldn]) { 
	    dn=1; di[0]=i; dj[0]=j; 
	  }
	  else if (chunk==1 || (chunk==0 && tp->dt_right[oldn])) { 
	    dn=2; 
	    di[0]=i; 
	    dj[0]=j; 
	    di[1]=i; 
	    dj[1]=j+1; 
	  }
          else if (chunk==2) { dn=2; di[0]=i; dj[0]=j; di[1]=i+1; dj[1]=j; }
          else { dn=4; di[0]=i; dj[0]=j;   di[1]=i; dj[1]=j+1; 
	  di[2]=i; dj[2]=j+1; di[3]=i+1; dj[3]=j+1; }
	
	  if (dn > 1) {
	    // determine the ordering of tiles to be removed
	    order_removals(tp,fp,dn,di,dj,removals);
	  }
	  else {
	    removals[0] = 0;
	  }
          for (d=0; d<dn; d++) { // delete each tile to be removed
	    i=di[removals[d]]; j=dj[removals[d]]; if (periodic) { i=(i+size)%size; j=(j+size)%size; }
	    oldns[removals[d]]=fp->Cell(i,j); // must make sure oldn is correct for each tile in chunk
	    if (i==fp->seed_i && j==fp->seed_j)
	      printf("removing seed at %d, %d! chunk=%d from %d,%d\n",i,j,chunk,di[0],dj[0]);
	    change_cell(fp,i,j,0);
	    if (!locally_fission_proof(fp,i,j,oldns[removals[d]])) { /* couldn't quickly confirm... */
	      if (flake_fission(fp,i,j)) {
		if (fission_allowed==0) {
		  // If we are collecting seen states, remove this
		  // state from the list of seen states, since it was
		  // never really "seen".
		  if (tp->tracking_seen_states && 
		      !between_double_tile(fp,tp,i,j,0)) {
		    remove_assembly_from_seen (tp);
		  }
		  
		  // re-attach cell:
		  // dissociation was chosen, but rejected because it would
		  // cause fission. (note that flake_fission calculates but
		  // doesn't remove cells if fission_allowed==0.)

		  for (k=0; k<=d; k++) {
		    change_cell(fp,di[removals[k]],dj[removals[k]],oldns[removals[k]]); tp->stat_a--; tp->stat_d--;
		  }
		  // If we are watching states to count how often they are entered, we 
		  // didn't actually leave the state we thought we left.
		  if (tp->watching_states && fp->chain_state) {
		    undo_state_off_indicator(fp);
		  }
		  break;
		}
	      } else if (0) { // should be fission_proof, according to local test
		// this is time consuming!  here only for debugging!
		if (flake_fission(fp,i,j)) {
		  printf("Fission_proof locale fissioned at %d,%d [chunk cell %d/%d]: "
			 "\n %3d %3d %3d\n %3d %3d %3d\n %3d %3d %3d\n", i,j, d, dn,
			 fp->Cell(i-1,j-1), fp->Cell(i-1,j  ), fp->Cell(i-1,j+1), 
			 fp->Cell(i  ,j-1),     oldn         , fp->Cell(i  ,j+1), 
			 fp->Cell(i+1,j-1), fp->Cell(i+1,j  ), fp->Cell(i+1,j+1) );
		}
	      }
	    }
	  } i=di[0]; j=dj[0];
	}
      } 
    } // else dprintf("can't move seed!\n"); NEW: no error, since this event exists
    else {
      assert (0);
    }
    // If we are using tinybox and only a single tile is left, remove the flake
    if (tp->tinybox && 
	(fp->tiles == 1 || (fp->tiles == 2 && fp->seed_is_double_tile))) {
      remove_flake(fp);
      //printf("Removed flake.  There are now %d flakes remaining.\n",tp->num_flakes);
    }
    d2printf("%d,%d -> %d\n",i,j,n);
    } // end of kTAM / aTAM section
    if (tp->flake_tree) 
      total_rate = tp->flake_tree->rate+tp->k*tp->conc[0]*tp->flake_tree->empty;
    else
      total_rate = 0;
    assert (total_rate >= 0);
    total_blast_rate = tp->k*tp->conc[0]*blast_rate*size*size*tp->num_flakes;
   } // end while
} // simulate


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
} // linear_simulate
 
    
