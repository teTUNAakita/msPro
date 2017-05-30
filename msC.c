/***** ms.c     ************************************************
*
*       Generates samples of gametes ( theta given or fixed number
*						of segregating sites.)
*	Usage is shown by typing ms without arguments.
usage: ms nsam howmany  -t  theta  [options]
or
ms nsam howmany -s segsites  [options]

nsam is the number of gametes per sample.
howmany is the number of samples to produce.
With -t the numbers of segregating sites will randomly vary
from one sample to the next.
with -s segsites,  the number of segregating sites will be
segsites in each sample.

Other options: See msdoc.pdf or after downloading and compiling, type ms<CR>.


*	  Arguments of the options are explained here:

npop:  Number of subpopulations which make up the total population
ni:  the sample size from the i th subpopulation (all must be
specified.) The output will have the gametes in order such that
the first n1 gametes are from the first island, the next n2 are
from the second island, etc.
nsites: number of sites between which recombination can occur.
theta: 4No times the neutral mutation rate
rho: recombination rate between ends of segment times 4No
f: ratio of conversion rate to recombination rate. (Wiuf and Hein model.)
track_len:  mean length of conversion track in units of sites.  The
total number of sites is nsites, specified with the -r option.
mig_rate: migration rate: the fraction of each subpop made up of
migrants times 4No.
howmany: howmany samples to generate.

Note:  In the above definition, No is the total diploid population if
npop is one, otherwise, No is the diploid population size of each
subpopulation.
A seed file called "seedms" will be created  if it doesn't exist. The
seed(s) in this file will be modified by the program.
So subsequent runs
will produce new output.  The initial contents of seedms will be
printed on the second line of the output.
Output consists of one line with the command line arguments and one
line with the seed(s).
The samples appear sequentially following that line.
Each sample begins with "//", then the number of segregating sites, the positions
of the segregating sites (on a scale of 0.0 - 1.0). On the following
lines are the sampled gametes, with mutants alleles represented as
ones and ancestral alleles as zeros.
To compile:  cc -o ms  ms.c  streec.c  rand1.c -lm
or:  cc -o ms ms.c streec.c rand2.c -lm
(Of course, gcc would be used instead of cc on some machines.  And -O3 or
some other optimization switches might be usefully employed with some
compilers.) ( rand1.c uses drand48(), whereas rand2.c uses rand() ).

*
*   Modifications made to combine ms and mss on 25 Feb 2001
*	Modifications to command line options to use switches  25 Feb 2001
*	Modifications to add // before each sample  25 Feb 2001
Modifications to add gene conversion 5 Mar 2001
Added demographic options -d  13 Mar 2001
Changed ran1() to use rand(). Changed seed i/o to accomodate this change. 20 April.
Changed cleftr() to check for zero rand() .13 June 2001
Move seed stuff to subroutine seedit()  11 July 2001
Modified streec.c to handle zero length demographic intervals 9 Aug 2001
Corrected problem with negative growth rates (Thanks to D. Posada and C. Wiuf) 13 May 2002
Changed sample_stats.c to output thetah - pi rather than pi - thetah.  March 8 2003.
Changed many command line options, allowing arbitrary migration matrix, and subpopulation
sizes.  Also allows parameters to come from a file. Option to output trees.  Option to
split and join subpopulations.   March 8, 2003. (Old versions saved in msold.tar ).
!!! Fixed bug in -en option.  Earlier versions may have produced garbage when -en ... used. 9 Dec 2003
Fixed bug which resulted in incorrect results for the case where
rho = 0.0 and gene conversion rate > 0.0. This case was not handled
correctly in early versions of the program. 5 Apr 2004.  (Thanks to
Vincent Plagnol for pointing out this problem.)
Fixed bug in prtree().  Earlier versions may have produced garbage when the -T option was used.
1 Jul 2004.
Fixed bug in -e. options that caused problems with -f option  13 Aug 2004.
Fixed bug in -es option, which was a problem when used with -eG. (Thanks again to V. Plagnol.) 6 Nov. 2004
Added -F option:  -F minfreq  produces output with sites with minor allele freq < minfreq filtered out.  11 Nov. 2004.
Fixed bug in streec.c (isseg() ).  Bug caused segmentation fault, crash on some machines. (Thanks
to Melissa Jane Todd-Hubisz for finding and reporting this bug.)
Added -seeds option 4 Nov 2006
Added "tbs" arguments feature 4 Nov 2006
Added -L option.  10 May 2007
Changed -ej option to set Mki = 0 pastward of the ej event.  See msdoc.pdf.  May 19 2007.
fixed bug with memory allocation when using -I option. This caused problems expecially on Windows
machines.  Thanks to several people, including Vitor Sousa and Stephane De Mita for help on this one.
Oct. 17, 2007.
Modified pickb() and pickbmf() to eliminate rare occurrence of fixed mutations Thanks to J. Liechty and K. Thornton. 4 Feb 2010.
***************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "ms.h"

#define SITESINC 100000

unsigned maxsites = SITESINC ;

// N E W ========
unsigned max_out = 100000 ; //later, must manage this value only when >100000
int **tmp = NULL ;
double *time_out = NULL ;
double *time_mut = NULL ;
int *node_mut = NULL ;
int **full_list = NULL ;
int *nch = NULL ;
int *dive = NULL ;
// N E W =====end


// N E W 2 ========
unsigned max_len = 100000 ;
double **two_dim = NULL ;
double *divergence = NULL ;
// N E W 2 ========end





struct node {
  int abv ;
  int ndes ;
  float time ;
} ;

struct segl {
  int beg ;
  struct node *ptree ;
  int next ;
} ;

double *posit ;
double segfac ;
int count, ntbs, nseeds ;
struct params pars ;

main ( argc,argv )
int argc ;
char *argv[] ;
{
  int i, j, k, howmany, segsites ;
  char **list, **cmatrix(), **tbsparamstrs ;
  FILE *pf, *fopen () ;
  double probss, tmrca, ttot ;
  void seedit ( const char * ) ;
  void getpars ( int argc, char *argv[], int *howmany, char* )  ; //st
  int gensam ( char **list, double *probss, double *ptmrca, double *pttot, int **tmp, double *time_out,
    int **full_list, double *time_mut, int *node_mut, int *nch, int *dive, double **two_dim,
    double *divergence, double **tractl ) ;

    char tract_div_file[100] ; //st	filename of tract and divergence
    int st_tl, en_tl ; //st for cleft2
    double qq ; //st for cleft2

    ntbs = 0 ;   /* these next few lines are for reading in parameters from a file (for each sample) */
    tbsparamstrs = (char **) malloc ( argc*sizeof(char *) ) ;

    for ( i=0; i<argc; i++ ){ printf("%s ",argv[i]) ; }
    for ( i=0; i<argc; i++ ){ tbsparamstrs[i] = (char *)malloc(30*sizeof(char) ) ; }
    for ( i=1; i<argc; i++ ){
      if ( strcmp( argv[i],"tbs") == 0 ){argv[i] = tbsparamstrs[ntbs++] ; }
    }
    count = 0 ;

    if ( ntbs > 0 )  for( k=0; k<ntbs; k++)  scanf(" %s", tbsparamstrs[k] ) ;

    // read and set parameters
    getpars ( argc, argv, &howmany, tract_div_file ) ;   /* results are stored in global variable, pars */ //st add filename

    //printf ("\n") ;
    //printf ("1 %s\n", tract_div_file) ; exit (0) ;

    //seeding
    if ( !pars.commandlineseedflag ) seedit ("s") ;
    pf = stdout ; // output

    // allocate memory for sequence list
    if ( pars.mp.segsitesin ==  0 ) { // use -t option
      list = cmatrix (pars.cp.nsam, maxsites+1) ;
      posit = (double *)malloc( (unsigned)( maxsites*sizeof( double)) ) ;
    }
    else { // use -s option
      list = cmatrix (pars.cp.nsam, pars.mp.segsitesin+1 ) ;
      posit = (double *)malloc( (unsigned)( pars.mp.segsitesin*sizeof( double)) ) ;
      if ( pars.mp.theta > 0.0 ){ // for calculating "prob:" by -t -s option
      segfac = 1.0 ;
      for ( i=pars.mp.segsitesin; i > 1; i-- ) segfac *= i ;
    }
  }

  // N E W 2 ========
  // Definition of matrix of tract and divergence that is given by user
  //st NOTE: two_dim[][2] is used for a different probability later
  two_dim = (double **) malloc ( max_len * sizeof(double *) ) ;
  //st comment out
  //for ( i = 0; i < max_len ; i++ ) {
  //    two_dim[i] = (double *) malloc ( 5 * sizeof(double) ) ; // len, div, freq, prob, accum
  //}
  if ( two_dim == NULL ) perror ( "malloc error" ) ;



  FILE *fp ;
  fp = fopen ( tract_div_file, "r" ) ;  //st
  double c1, c2, c3 ;
  double sum_c, tmp_c ;

  //st species info
  double **species ; // including divergence, max_tract, relative proportion
  species = (double **) malloc ( 11 * sizeof(double *) ) ;//ta 10=>11
  for ( i=0; i < 11; i++ ){//ta 10=>11
    species[i] = (double *) malloc ( 3 * sizeof(double) ) ;
    species[i][0] = 0.0 ;
    species[i][1] = -777.0 ;
    species[i][2] = 0.0 ;
  }
  int n_species = 0 ;
  int f_species = 0 ;

  // Read and strore the data
  if ( fp == NULL ){
    printf ( "cannot open %s\n", tract_div_file ) ; //st
    exit (1) ;
  }
  else { //st
    //printf ( "two_dim: input date with scaling\ntract\tdiv\tfreq\tprob\tcumm_prob\n" ) ;
    j = 0 ;
    sum_c = 0 ;
    while ( fscanf( fp, "%lf%lf%lf", &c1, &c2, &c3 ) == 3 ) {
      //printf ( "%.0lf\t%.2lf\t%lf\n", c1, c2, c3 ) ;
      two_dim[j] = (double *) malloc ( 5 * sizeof(double) ) ; //st
      two_dim[j][0] = c1 ; // tract
      two_dim[j][1] = c2 ; // diverfence
      two_dim[j][2] = c3 ; // proportion => changing prob of tract lenght for cleft2 later
      two_dim[j][3] = 0.0 ; //st probability
      two_dim[j][4] = 0.0 ; //st cumulative probability
      sum_c += c3 ;
      j++ ;
      //printf("j=%d\n",j);
      //st Error check
      if ( c1 == 0.0 ){
        printf ("tract length must be > 0\n") ;
        exit (1) ;
      }
      if ( (c2 == 0.0) || (1.0 <= c2) ){
        printf ("Divergence value, d, must be 0.0 < d < 1.0\n") ;
      }

      //st species info
      if ( n_species == 0 ){
        species[0][0] = c2 ;
        species[0][1] = c1 ;
        species[0][2] = c3 ;
        n_species++ ;
      }
      else {
        f_species = 0 ;
        for ( i=0; i < n_species; i++ ){
          if ( species[i][0] == c2 ){ // species already exists
            species[i][2] += c3 ;
            if ( species[i][1] < c1 ){
              species[i][1] = c1 ;
            }
            f_species = 1 ;
            break ;
          }
        }
        if ( f_species == 0 ){ // species does not exist
          species[n_species][0] = c2 ;
          species[n_species][1] = c1 ;
          species[n_species][2] = c3 ;
          n_species++ ;
          if ( n_species == 11 ){
            printf ("No. species must be <= 10\n") ;
            exit (1) ;
          }
        }

      }

      // st memory allocation
      if ( j == max_len ){ //st
        max_len += 100 ;
        two_dim = (double **) realloc ( two_dim, max_len * sizeof(double *) ) ;
      }
    }
    fclose ( fp ) ;
  }
  pars.cp.n_two_dim = j ;

  // Manipulate the matrix
  // sum of prob (two_dim[][3]) must be 1
  // cumulative prob in two_dim[][4]
  tmp_c = 0.0 ;
  for ( i=0; i < pars.cp.n_two_dim; i++ ) {
    two_dim[i][3] = two_dim[i][2] / sum_c ;
    tmp_c += two_dim[i][3] ;
    //printf ("sum_c = %f, tmp_c = %f\n",sum_c, tmp_c ) ;
    two_dim[i][4] = tmp_c ;
    //printf ( "%.0f\t%.2f\t%.3f\t%.3f\t%.3f\n", two_dim[i][0], two_dim[i][1], two_dim[i][2], two_dim[i][3], two_dim[i][4] ) ;
  }
  //st
  tmp_c = 0.0 ;
  int n_tractl = 0 ; // # of **tractl
  for ( i=0; i < n_species; i++ ){
    tmp_c += species[i][2] ;
    n_tractl += (int) species[i][1] ;
  }
  for ( i=0; i < n_species; i++ ){
    species[i][2] /= tmp_c ;
  }
  pars.cp.n_tractl = n_tractl ;

  //printf ( "length of given matrix = %d\n", j ) ;
if(0){
  for ( i=0; i < pars.cp.n_two_dim; i++ ) {
    printf ("%d\t%f\t%f\t%f\t%f\n", (int) two_dim[i][0], two_dim[i][1], two_dim[i][2], two_dim[i][3], two_dim[i][4]) ;
  }
  printf ("no. two_dim = %d\n", pars.cp.n_two_dim) ;
  printf ("no. species = %d\n", n_species) ;
  for ( i=0; i < n_species; i++ ){
    printf ("  %dth, d = %.2f, maxq=%d, ratio=%f\n", i, species[i][0], (int) species[i][1], species[i][2]) ;
  }
  printf ("no. tractl = %d\n", n_tractl) ;
}
  //st calculate the distribution of tract length for cleft2
  double **tractl ; // including tract length, divergence, probability, cumulative prob
  tractl = (double **) malloc ( n_tractl * sizeof(double *) ) ;

  //st
  // initial setting
  // two_dim[][2] is used as a conditional probability
  // Probability is normalized for tract length
  // The longer tract length is, the higher probability of cleft2 is
  tmp_c = 0.0 ;
  for ( k=0; k < pars.cp.n_two_dim; k++ ){
    tmp_c += two_dim[k][0] ;
  }
  for ( k=0; k < pars.cp.n_two_dim; k++ ){
    two_dim[k][2] = two_dim[k][3]*two_dim[k][0]/tmp_c ;
  }

  st_tl = 0 ;
  en_tl = 0 ;
  sum_c = 0.0 ;
  for ( i=0; i < n_species; i++ ){
    en_tl = st_tl + species[i][1] - 1 ;
    //printf ("%d %d\n", st_tl, en_tl) ;

    // calculate probability of cleft2
    qq = 1.0 ;
    for ( j=st_tl; j <= en_tl; j++ ){
      tractl[j] = (double *) malloc ( 4 * sizeof(double) ) ; // tract, divergence, prob, accum
      tractl[j][0] = qq ;
      tractl[j][1] = species[i][0] ;
      tractl[j][2] = 0.0 ;
      // proportion
      for ( k=0; k < pars.cp.n_two_dim; k++ ){
        if ( (two_dim[k][1] == species[i][0]) && (tractl[j][0] <= two_dim[k][0]) ){
          tractl[j][2] += ( two_dim[k][2] * (1.0/two_dim[k][0]) ) ;
        }
      }
      sum_c += tractl[j][2] ;
      tractl[j][3] = 0.0 ;
      qq += 1.0 ;
    }
    st_tl = en_tl + 1 ;
  }

  //st normalize and calculate cumulative probability
  tmp_c = 0.0 ;
  for ( i=0; i < n_tractl; i++ ){
    tractl[i][2] /= sum_c ;
    tmp_c += tractl[i][2] ;
    tractl[i][3] = tmp_c ;
  }

  FILE *hoge ;
  if ( (hoge = fopen ("hoge.out", "w") ) == NULL ){printf ("Error 1 for file handle\n"); exit(1); }
  for ( i=0; i < n_tractl; i++ ){
    fprintf (hoge, "%d\t%.1f\t%f\t%f\n", (int) tractl[i][0], tractl[i][1], tractl[i][2], tractl[i][3]) ;
  }
  fclose (hoge) ;

  //printf ("OK\n") ;
  //exit (0) ;


  while ( howmany-count++ ) {
    // N E W ========
    if ( nch == NULL ){
      nch = (int *) malloc( 1 * sizeof(int) ) ;
      if ( nch == NULL ) perror("malloc error. segtre_mig") ;
    }
    if ( time_out == NULL ){
      time_out = (double *) malloc( max_out * sizeof(double) ) ;
      if ( time_out == NULL ) perror("malloc error. segtre_mig") ;
    }
    if ( divergence == NULL ){
      divergence = (double *) malloc( max_out * sizeof(double) ) ;
      if ( divergence == NULL ) perror("malloc error. segtre_mig") ;
    }

    if ( tmp == NULL ){
      tmp = (int **) malloc( max_out * sizeof(int *) ) ;
      for ( i = 0; i < max_out; i++ ){
        tmp[i] = (int *) malloc ( (pars.cp.nsam + 4) * sizeof(int) );// [start] + [end] + [nsam] + [in] + [out]
      }
      if ( tmp == NULL ) perror("malloc error. segtre_mig") ;
    }

    if ( full_list == NULL ){
      full_list = (int **) calloc( pars.cp.nsam, sizeof(int *) ) ;
      for ( i = 0; i < pars.cp.nsam; i++ ){
        full_list[i] = (int *) calloc ( pars.cp.nsites, sizeof(int) );
      }
      if ( full_list == NULL ) perror("malloc error. segtre_mig") ;
    }
    for( i = 0; i < pars.cp.nsam; i++ ){
      for( j = 0; j < pars.cp.nsites; j++ ){
        full_list[i][j] = 0 ;
        //printf("%d", full_list[i][j]);
      }
      //printf("\n") ;
    }

    if ( time_mut == NULL ){
      time_mut = (double *) calloc( maxsites+1, sizeof(double) ) ;
      if ( time_mut == NULL ) perror("malloc error. segtre_mig") ;
    }
    if ( node_mut == NULL ){
      node_mut = (int *) calloc( maxsites+1, sizeof(int) ) ;
      if ( node_mut == NULL ) perror("malloc error. segtre_mig") ;
    }

    if ( dive == NULL ){
      dive = (int *) calloc( pars.cp.nsites, sizeof(int) ) ;
      if ( dive == NULL ) perror("malloc error. segtre_mig") ;
    }
    // N E W =====end

    // not used
    if ( (ntbs > 0) && (count >1 ) ){
      for ( k=0; k < ntbs; k++){
        if ( scanf (" %s", tbsparamstrs[k]) == EOF ){if ( !pars.commandlineseedflag ) seedit( "end" ); exit(0) ; }
      }
      getpars ( argc, argv, &howmany, tract_div_file) ;
    }
    //following is needed!!

    //fprintf (pf,"\n//") ;
    if ( ntbs > 0 ){for ( k=0; k< ntbs; k++ ) printf("\t%s", tbsparamstrs[k] ) ; } // not used
    //printf("\n");

    // generate samples
    segsites = gensam ( list, &probss, &tmrca, &ttot, tmp, time_out, full_list, time_mut, node_mut, nch, dive,
      two_dim, divergence, tractl ) ;//st add tractl

      // output res
      if ( pars.mp.timeflag ){ // -L total length of tree
        fprintf(pf,"time:\t%lf\t%lf\n",tmrca, ttot ) ;
      }
      // list[sample][site]
      if ( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) { // there is seg site or 4Nu > 0
        if ( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 )){ // -t & -s
          fprintf(pf,"prob: %g\n", probss ) ;
        }
        // # of seg sites
        //following is needed!!
        fprintf (pf,"segsites: %d, nch: %d\n", segsites, *nch) ;
        // potition
        /*
        if ( segsites > 0 ){
        fprintf(pf,"positions: ") ;
      }

      for ( i=0; i < segsites; i++){
      fprintf(pf,"%6.4lf ",posit[i] );
    }
    fprintf(pf,"\n");
    if ( segsites > 0 ){
    for ( i=0;i < pars.cp.nsam; i++ ){
    fprintf (pf,"%s\n", list[i] );
  }
}
*/
//fprintf (pf,"\nhoge\n");
//fprintf (pf,"%s\n", list[0]);

}
// N E W ========
/*
for (i = 0; i < segsites; i++){
printf("\t%d", node_mut[i]);
}
printf("\n");
for (i = 0; i < segsites; i++){
printf("\t%.2f", time_mut[i]);
}
printf("\n\n");
*/

//following is needed!!

for (i = 0; i < pars.cp.nsam; i++){
  for (j = 0; j < pars.cp.nsites; j++) {
    printf("%d", full_list[i][j]);
  }
  printf("\n");
}
printf("\n");

// N E W =====end



}


// THIS LINE IS NEEDED IF RNADOM NUMBER IS CHANGED FOR RECURREATION
if( !pars.commandlineseedflag ) seedit( "end" );
}


// segsites = gensam ( list, &probss, &tmrca, &ttot ) ;
int gensam ( char **list, double *pprobss, double *ptmrca, double *pttot, int **tmp, double *time_out, int **full_list,
  double *time_mut, int *node_mut, int *nch, int *dive, double **two_dim, double *divergence,
  double **tractl ){//st add tractl
    int nsegs, h, i, k, j, seg, ns, start, end, len, segsit ;
    struct segl *seglst, *segtre_mig (struct c_params *p, int *nsegs, int **tmp, double *time_out, int *nch, double **two_dim,
      double *divergence, double **tractl ) ;
      /* used to be: [MAXSEG];  */ //ADD int **tmp, int max_out
      double nsinv,  tseg, tt, ttime (struct node *, int nsam), ttimemf (struct node *, int nsam, int mfreq) ;
      double *pk;
      int *ss;
      int segsitesin,nsites;
      double theta, es ;
      int nsam, mfreq ;
      void prtree( struct node *ptree, int nsam);
      //int make_gametes(int nsam, int mfreq,  struct node *ptree, double tt, int newsites, int ns, char **list );
      void ndes_setup( struct node *, int nsam );

      // N E W ========
      int make_gametes2(int nsam, int mfreq,  struct node *ptree, double tt, int newsites, int ns, char **list, int **full_list,
        double *time_mut, int *node_mut, int seg, int start, int end, int *dive );
        // N E W =====end



        nsites = pars.cp.nsites ; // length
        nsinv = 1./nsites;        // ??

        /* generate tree */
        seglst = segtre_mig ( &(pars.cp), &nsegs, tmp, time_out, nch, two_dim, divergence, tractl ) ; // generate tree in streec.c, //ADD int **tmp, int max_out

        nsam = pars.cp.nsam; // # of smaples
        segsitesin = pars.mp.segsitesin ; // S
        theta = pars.mp.theta ; // theta per region
        mfreq = pars.mp.mfreq ; // minor allele freq??

        if ( pars.mp.treeflag ) { // output tree phy
          ns = 0 ;
          for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
            if( (pars.cp.r > 0.0 ) || (pars.cp.f > 0.0) ){
              end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
              start = seglst[seg].beg ;
              len = end - start + 1 ;
              fprintf (stdout,"[%d]", len);
            }
            prtree( seglst[seg].ptree, nsam ) ;
            if( (segsitesin == 0) && ( theta == 0.0 ) && ( pars.mp.timeflag == 0 ) )
            free(seglst[seg].ptree) ;
          }
        }

        if ( pars.mp.timeflag ) { // calculate length of tree
          tt = 0.0 ;
          for ( seg=0, k=0; k < nsegs; seg=seglst[seg].next, k++ ) {
            if ( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam ) ; // default is 1
            end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 ) ;
            start = seglst[seg].beg ;
            if ( (nsegs == 1) || ( ( start <= nsites/2) && ( end >= nsites/2 ) ) ){
              *ptmrca = (seglst[seg].ptree + 2*nsam-2) -> time ;
            }
            len = end - start + 1 ;
            tseg = len/(double)nsites ;
            if ( mfreq == 1 ){
              tt += ttime(seglst[seg].ptree,nsam)*tseg ;
            }
            else {
              tt += ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
            }
            if ( (segsitesin == 0) && ( theta == 0.0 )  ){
              free (seglst[seg].ptree) ;
            }
          }
          *pttot = tt ;
        }

        /*
        ===========
        */
        if ( (segsitesin == 0) && (theta > 0.0)  ) { // use only -t option
          ns = 0 ;
          for ( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++ ) {
            // printf (" seg = %d, k = %d\n", seg,k) ; only once
            if ( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam ) ; // minor allele freq
            end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 ) ;
            start = seglst[seg].beg ;
            len = end - start + 1 ; // len = 2 without recombination
            //printf (" len = %d\n", len) ;
            tseg = len*(theta/nsites) ;
            if ( mfreq == 1 ){ // default
              tt = ttime (seglst[seg].ptree, nsam) ; // total tree length
            }
            else {
              tt = ttimemf (seglst[seg].ptree, nsam, mfreq ) ;
            }
            segsit = poisso ( tseg*tt ) ; // number of mutations
            if ( (segsit + ns) >= maxsites ) {
              maxsites = segsit + ns + SITESINC ;
              posit = (double *)realloc(posit, maxsites*sizeof(double) ) ;
              biggerlist(nsam, list) ;
            }

            // assign mutations on the genealogy and
            // generate sequence data
            //make_gametes (nsam, mfreq, seglst[seg].ptree, tt, segsit, ns, list ) ;
            make_gametes2 (nsam, mfreq, seglst[seg].ptree, tt, segsit, ns, list, full_list, time_mut, node_mut, seg, start, end, dive ) ;
            free (seglst[seg].ptree) ;

            // determine the position of the mutation
            //locate (segsit,start*nsinv, len*nsinv,posit+ns) ;

            // number of segregating sites incremented
            ns += segsit;
          }
        }
        else if ( segsitesin > 0 ) { // use -s option

          pk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
          ss = (int *)malloc((unsigned)(nsegs*sizeof(int)));
          if( (pk==NULL) || (ss==NULL) ) perror("malloc error. gensam.2");


          tt = 0.0 ;
          for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
            if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
            end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
            start = seglst[seg].beg ;
            len = end - start + 1 ;
            tseg = len/(double)nsites ;
            if( mfreq == 1 ) pk[k] = ttime(seglst[seg].ptree,nsam)*tseg ;
            else pk[k] = ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
            tt += pk[k] ;
          }
          if( theta > 0.0 ) {
            es = theta * tt ;
            *pprobss = exp( -es )*pow( es, (double) segsitesin) / segfac ;
          }
          if( tt > 0.0 ) {
            for (k=0;k<nsegs;k++) pk[k] /= tt ;
            mnmial(segsitesin,nsegs,pk,ss);
          }
          else
          for( k=0; k<nsegs; k++) ss[k] = 0 ;
          ns = 0 ;
          for ( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
            end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 ) ;
            start = seglst[seg].beg ;
            len = end - start + 1 ;
            tseg = len/(double)nsites;
            make_gametes (nsam,mfreq,seglst[seg].ptree,tt*pk[k]/tseg, ss[k], ns, list);

            free(seglst[seg].ptree) ;
            locate(ss[k],start*nsinv, len*nsinv,posit+ns);
            ns += ss[k] ;
          }
          free (pk) ;
          free (ss) ;

        }

        // last of strings
        /*
        for (i=0;i<nsam;i++) list[i][ns] = '\0' ;
        */
        return( ns ) ;
      }

      void ndes_setup(struct node *ptree, int nsam ){
        int i ;

        for( i=0; i<nsam; i++) (ptree+i)->ndes = 1 ;
        for( i = nsam; i< 2*nsam -1; i++) (ptree+i)->ndes = 0 ;
        for( i= 0; i< 2*nsam -2 ; i++)  (ptree+((ptree+i)->abv))->ndes += (ptree+i)->ndes ;

      }

      int biggerlist ( nsam,  list )
      int nsam ;
      char ** list ;
      {
        int i;

        /*  fprintf(stderr,"maxsites: %d\n",maxsites);  */
        for( i=0; i<nsam; i++){
          list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
          if( list[i] == NULL ) perror( "realloc error. bigger");
        }
      }



      /* allocates space for gametes (character strings) */
      // list = cmatrix (pars.cp.nsam, maxsites+1) ;
      char **cmatrix( nsam, len )
      int nsam, len ;
      {
        int i ;
        char **m ; // output

        // list[sample][site]
        if ( ! ( m = (char **) malloc ( (unsigned) nsam*sizeof( char* ) ) ) ){
          perror ("alloc error in cmatrix") ;
        }
        for( i=0; i < nsam; i++ ) {
          if( ! ( m[i] = (char *) malloc( (unsigned) len*sizeof( char ) ))){
            perror("alloc error in cmatric. 2");
          }
        }
        return ( m ) ;
      }



      int
      locate(n,beg,len,ptr)
      int n;
      double beg, len, *ptr;
      {
        int i;

        ordran(n,ptr);
        for(i=0; i<n; i++)
        ptr[i] = beg + ptr[i]*len ;

      }

      int NSEEDS = 3 ;

      void getpars ( int argc, char *argv[], int *phowmany, char *tdfile ){

        int arg, i, j, sum , pop , argstart, npop , npop2, pop2 ;
        double migr, mij, psize, palpha ;
        void addtoelist( struct devent *pt, struct devent *elist );
        void argcheck( int arg, int argc, char ** ) ;
        int commandlineseed( char ** ) ;
        void free_eventlist( struct devent *pt, int npop );
        struct devent *ptemp , *pt ;
        FILE *pf ;
        char ch3 ;


        if ( count == 0 ) { // initial setting or if tbs is not used
          if ( argc < 4 ){fprintf(stderr,"Too few command line arguments\n") ; usage() ; }
          pars.cp.nsam = atoi( argv[1] ) ; // # of samples
          if ( pars.cp.nsam <= 0 ) {fprintf(stderr,"First argument error. nsam <= 0. \n") ; usage() ; }
          *phowmany = atoi( argv[2] ) ; // # of rep
          if ( *phowmany  <= 0 ) {fprintf(stderr,"Second argument error. howmany <= 0. \n") ; usage() ; }
          pars.commandlineseedflag = 0 ; //
          pars.cp.r = pars.mp.theta =  pars.cp.f = 0.0 ; // g/f
          pars.cp.track_len = 0. ;  // 1/q

          // N E W ========
          pars.cp.f2 = 0.0 ; // g'
          pars.cp.track_len2 = 0. ;  // 1/q'
          pars.cp.divergence = 0. ;  // 1/q2
          // N E W =====end

          // N E W 2 ======
          pars.cp.size_tract = 0 ;
          pars.cp.size_div = 0 ;
          // N E W 2 ======


          pars.cp.npop = npop = 1 ;
          pars.cp.mig_mat = (double **)malloc( (unsigned) sizeof(double*));
          pars.cp.mig_mat[0] = (double *)malloc( (unsigned)sizeof(double));
          pars.cp.mig_mat[0][0] =  0.0 ;
          pars.mp.segsitesin = 0 ;  // -s
          pars.mp.treeflag = 0 ; // -T tree output
          pars.mp.timeflag = 0 ; // -L length of tree
          pars.mp.mfreq = 1 ;    // ?minor allele?
          pars.cp.config = (int *) malloc((unsigned)(( pars.cp.npop +1 ) *sizeof(int)));
          (pars.cp.config)[0] = pars.cp.nsam ;
          pars.cp.size = (double *) malloc((unsigned)( pars.cp.npop *sizeof(double)));
          (pars.cp.size)[0] = 1.0  ;
          pars.cp.alphag = (double *) malloc((unsigned)(( pars.cp.npop ) *sizeof(double)));
          (pars.cp.alphag)[0] = 0.0  ;
          pars.cp.nsites = 2 ; // Length
        }
        else{ // if tbs is used
          npop = pars.cp.npop ;
          free_eventlist( pars.cp.deventlist, npop ) ;
        }
        pars.cp.deventlist = NULL ;

        arg = 3 ;

        int lenf ;

        while ( arg < argc ){
          if ( argv[arg][0] != '-' ){fprintf (stderr," argument should be -%s ?\n", argv[arg]) ; usage() ; }
          switch ( argv[arg][1] ){
            case 'f' : // arguments from FILE
            if ( ntbs > 0 ) { fprintf(stderr," can't use tbs args and -f option.\n"); exit(1); }
            arg++ ;
            argcheck ( arg, argc, argv) ;
            pf = fopen( argv[arg], "r" ) ;
            if( pf == NULL ) {fprintf(stderr," no parameter file %s\n", argv[arg] ); exit(0);}
            arg++;
            argc++ ;
            argv = (char **)malloc(  (unsigned)(argc+1)*sizeof( char *) ) ;
            argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
            argstart = arg ;
            while( fscanf(pf," %s", argv[arg]) != EOF ) {
              arg++;
              argc++;
              argv = (char **)realloc( argv, (unsigned)argc*sizeof( char*) ) ;
              argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
            }
            fclose(pf);
            argc--;
            arg = argstart ;
            break;
            case 'r' : // 4Nr*(len-1)
            arg++;
            argcheck ( arg, argc, argv) ;
            pars.cp.r = atof ( argv[arg++] ) ;
            argcheck( arg, argc, argv) ;
            pars.cp.nsites = atoi( argv[arg++]) ;
            if ( pars.cp.nsites < 2 ){
              fprintf (stderr,"with -r option must specify both rec_rate and nsites>1\n") ;
              usage ();
            }
            break;
            case 'c' : // g/r and 1/q
            arg++;
            argcheck( arg, argc, argv);
            pars.cp.f = atof(  argv[arg++] );
            argcheck( arg, argc, argv);
            pars.cp.track_len = atof( argv[arg++]);
            if( pars.cp.track_len <1. ){
              fprintf(stderr,"with -c option must specify both f and track_len>0\n");
              usage();
            }
            break;
            case 'b': //st
            arg++ ;
            argcheck( arg, argc, argv);
            pars.cp.f2 = atof(  argv[arg++] );
            strcpy (tdfile, argv[arg]) ;
            arg++ ;
            break;

            case 't' : // 4Nu
            arg++;
            argcheck( arg, argc, argv);
            pars.mp.theta = atof(  argv[arg++] );
            break;
            case 's' : // S
            arg++;
            argcheck (arg, argc, argv) ;
            if ( argv[arg-1][2] == 'e' ){  // command line seeds, not care
              pars.commandlineseedflag = 1 ;
              if( count == 0 ) nseeds = commandlineseed(argv+arg );
              arg += nseeds ;
            }
            else { // S
              pars.mp.segsitesin = atoi (argv[arg++]) ;
            }
            break;
            case 'F' :
            arg++;
            argcheck( arg, argc, argv);
            pars.mp.mfreq = atoi(  argv[arg++] );
            if( (pars.mp.mfreq < 2 ) || (pars.mp.mfreq > pars.cp.nsam/2 ) ){
              fprintf(stderr," mfreq must be >= 2 and <= nsam/2.\n");
              usage();
            }
            break;
            case 'T' : // -T tree output
            pars.mp.treeflag = 1 ;
            arg++;
            break;
            case 'L' : // -L length of tree
            pars.mp.timeflag = 1 ;
            arg++;
            break;
            case 'I' : // subpopulation
            arg++;
            if ( count == 0 ) {
              argcheck( arg, argc, argv);
              pars.cp.npop = atoi( argv[arg]);
              pars.cp.config = (int *) realloc ( pars.cp.config, (unsigned)( pars.cp.npop*sizeof( int)));
              npop = pars.cp.npop ;
            }
            arg++ ;
            for ( i=0; i< pars.cp.npop; i++) {
              argcheck( arg, argc, argv);
              pars.cp.config[i] = atoi( argv[arg++]);
            }
            if ( count == 0 ){
              pars.cp.mig_mat =
              (double **)realloc(pars.cp.mig_mat, (unsigned)(pars.cp.npop*sizeof(double *) )) ;
              pars.cp.mig_mat[0] =
              (double *)realloc(pars.cp.mig_mat[0], (unsigned)( pars.cp.npop*sizeof(double)));
              for(i=1; i<pars.cp.npop; i++) pars.cp.mig_mat[i] =
              (double *)malloc( (unsigned)( pars.cp.npop*sizeof(double)));
              pars.cp.size = (double *)realloc( pars.cp.size, (unsigned)( pars.cp.npop*sizeof( double )));
              pars.cp.alphag =
              (double *) realloc( pars.cp.alphag, (unsigned)( pars.cp.npop*sizeof( double )));
              for( i=1; i< pars.cp.npop ; i++) {
                (pars.cp.size)[i] = (pars.cp.size)[0]  ;
                (pars.cp.alphag)[i] = (pars.cp.alphag)[0] ;
              }
            }
            if ( (arg <argc) && ( argv[arg][0] != '-' ) ) {
              argcheck( arg, argc, argv);
              migr = atof(  argv[arg++] );
            }
            else migr = 0.0 ;
            for ( i=0; i<pars.cp.npop; i++)
            for ( j=0; j<pars.cp.npop; j++) pars.cp.mig_mat[i][j] = migr/(pars.cp.npop-1) ;
            for ( i=0; i< pars.cp.npop; i++) pars.cp.mig_mat[i][i] = migr ;
            break;
            case 'm' :
            if ( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
            if ( argv[arg][2] == 'a' ) {
              arg++;
              for( pop = 0; pop <npop; pop++)
              for( pop2 = 0; pop2 <npop; pop2++){
                argcheck( arg, argc, argv);
                pars.cp.mig_mat[pop][pop2]= atof( argv[arg++] ) ;
              }
              for( pop = 0; pop < npop; pop++) {
                pars.cp.mig_mat[pop][pop] = 0.0 ;
                for( pop2 = 0; pop2 < npop; pop2++){
                  if( pop2 != pop ) pars.cp.mig_mat[pop][pop] += pars.cp.mig_mat[pop][pop2] ;
                }
              }
            }
            else {
              arg++;
              argcheck( arg, argc, argv);
              i = atoi( argv[arg++] ) -1;
              argcheck( arg, argc, argv);
              j = atoi( argv[arg++] ) -1;
              argcheck( arg, argc, argv);
              mij = atof( argv[arg++] );
              pars.cp.mig_mat[i][i] += mij -  pars.cp.mig_mat[i][j]  ;
              pars.cp.mig_mat[i][j] = mij;
            }
            break;
            case 'n' : // size of subpop in t = 0, -n popID, x (x*N0)
            if ( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
            arg++;
            argcheck ( arg, argc, argv);
            pop = atoi( argv[arg++] ) -1;
            argcheck( arg, argc, argv);
            psize = atof( argv[arg++] );
            pars.cp.size[pop] = psize ;
            break;
            case 'g' :
            if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
            arg++;
            argcheck( arg, argc, argv);
            pop = atoi( argv[arg++] ) -1;
            if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); usage(); }
            palpha = atof( argv[arg++] );
            pars.cp.alphag[pop] = palpha ;
            break;
            case 'G' :
            arg++;
            if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); usage(); }
            palpha = atof( argv[arg++] );
            for( i=0; i<pars.cp.npop; i++)
            pars.cp.alphag[i] = palpha ;
            break;
            case 'e' :
            pt = (struct devent *)malloc( sizeof( struct devent) ) ;
            pt->detype = argv[arg][2] ;
            ch3 = argv[arg][3] ; // for migration
            arg++;
            argcheck (arg, argc, argv) ;
            pt->time = atof( argv[arg++] ) ;
            pt->nextde = NULL ;
            if ( pars.cp.deventlist == NULL )
            pars.cp.deventlist = pt ;
            else if ( pt->time < pars.cp.deventlist->time ) {
              ptemp = pars.cp.deventlist ;
              pars.cp.deventlist = pt ;
              pt->nextde = ptemp ;
            }
            else {
              addtoelist( pt, pars.cp.deventlist ) ;
            }
            switch ( pt->detype ) {
              case 'N' :
              argcheck( arg, argc, argv);
              pt->paramv = atof( argv[arg++] ) ;
              break;
              case 'G' :
              if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eG.\n"); usage(); }
              pt->paramv = atof( argv[arg++] ) ;
              break;
              case 'M' :
              argcheck( arg, argc, argv);
              pt->paramv = atof( argv[arg++] ) ;
              break;
              case 'n' : // -en, change the size of subpop
              argcheck( arg, argc, argv);
              pt->popi = atoi( argv[arg++] ) -1 ;
              argcheck( arg, argc, argv);
              pt->paramv = atof( argv[arg++] ) ;
              break;
              case 'g' :
              argcheck( arg, argc, argv);
              pt->popi = atoi( argv[arg++] ) -1 ;
              if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eg.\n"); usage(); }
              pt->paramv = atof( argv[arg++] ) ;
              break;
              case 's' :
              argcheck( arg, argc, argv);
              pt->popi = atoi( argv[arg++] ) -1 ;
              argcheck( arg, argc, argv);
              pt->paramv = atof( argv[arg++] ) ;
              break;
              case 'm' :
              if( ch3 == 'a' ) {
                pt->detype = 'a' ;
                argcheck( arg, argc, argv);
                npop2 = atoi( argv[arg++] ) ;
                pt->mat = (double **)malloc( (unsigned)npop2*sizeof( double *) ) ;
                for( pop =0; pop <npop2; pop++){
                  (pt->mat)[pop] = (double *)malloc( (unsigned)npop2*sizeof( double) );
                  for( i=0; i<npop2; i++){
                    if( i == pop ) arg++;
                    else {
                      argcheck( arg, argc, argv);
                      (pt->mat)[pop][i] = atof( argv[arg++] ) ;
                    }
                  }
                }
                for( pop = 0; pop < npop2; pop++) {
                  (pt->mat)[pop][pop] = 0.0 ;
                  for( pop2 = 0; pop2 < npop2; pop2++){
                    if( pop2 != pop ) (pt->mat)[pop][pop] += (pt->mat)[pop][pop2] ;
                  }
                }
              }
              else {
                argcheck( arg, argc, argv);
                pt->popi = atoi( argv[arg++] ) -1 ;
                argcheck( arg, argc, argv);
                pt->popj = atoi( argv[arg++] ) -1 ;
                argcheck( arg, argc, argv);
                pt->paramv = atof( argv[arg++] ) ;
              }
              break;
              case 'j' : // merge subpop
              argcheck( arg, argc, argv);
              pt->popi = atoi( argv[arg++] ) -1 ;
              argcheck( arg, argc, argv);
              pt->popj = atoi( argv[arg++] ) -1 ;
              break;
              default: fprintf(stderr,"e event\n");  usage();
            }
            break;
            default: fprintf(stderr," option default\n");  usage() ;
          }
        }
        if( (pars.mp.theta == 0.0) && ( pars.mp.segsitesin == 0 ) && ( pars.mp.treeflag == 0 ) && (pars.mp.timeflag == 0) ) {
          fprintf(stderr," either -s or -t or -T option must be used. \n");
          usage();
          exit(1);
        }

        sum = 0 ;
        for( i=0; i< pars.cp.npop; i++) sum += (pars.cp.config)[i] ;
        if( sum != pars.cp.nsam ) {
          fprintf(stderr," sum sample sizes != nsam\n");
          usage();
          exit(1);
        }
      }


      void
      argcheck( int arg, int argc, char *argv[] )
      {
        if( (arg >= argc ) || ( argv[arg][0] == '-') ) {
          fprintf(stderr,"not enough arguments after %s\n", argv[arg-1] ) ;
          fprintf(stderr,"For usage type: ms<return>\n");
          exit(0);
        }
      }

      int
      usage()
      {
        fprintf(stderr,"usage: ms nsam howmany \n");
        fprintf(stderr,"  Options: \n");
        fprintf(stderr,"\t -t theta   (this option and/or the next must be used. Theta = 4*N0*u )\n");
        fprintf(stderr,"\t -s segsites   ( fixed number of segregating sites)\n");
        fprintf(stderr,"\t -T          (Output gene tree.)\n");
        fprintf(stderr,"\t -F minfreq     Output only sites with freq of minor allele >= minfreq.\n");
        fprintf(stderr,"\t -r rho nsites     (rho here is 4Nc)\n");
        fprintf(stderr,"\t\t -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) \n");
        fprintf(stderr,"\t\t\t if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.\n");
        fprintf(stderr,"\t -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t\n");
        fprintf(stderr,"\t -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) \n");
        fprintf(stderr,"\t\t -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)\n");
        fprintf(stderr,"\t\t -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");
        fprintf(stderr,"\t\t -n i size_i   (popi has size set to size_i*N0 \n");
        fprintf(stderr,"\t\t -g i alpha_i  (If used must appear after -M option.)\n");
        fprintf(stderr,"\t   The following options modify parameters at the time 't' specified as the first argument:\n");
        fprintf(stderr,"\t -eG t alpha  (Modify growth rate of all pop's.)\n");
        fprintf(stderr,"\t -eg t i alpha_i  (Modify growth rate of pop i.) \n");
        fprintf(stderr,"\t -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n");
        fprintf(stderr,"\t -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )\n");
        fprintf(stderr,"\t -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");
        fprintf(stderr,"\t -eN t size  (Modify pop sizes. New sizes = size*N0 ) \n");
        fprintf(stderr,"\t -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
        fprintf(stderr,"\t -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");
        fprintf(stderr,"\t\t proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.\n");
        fprintf(stderr,"\t\t Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.\n");
        fprintf(stderr,"\t -ej t i j   ( Join lineages in pop i and pop j into pop j\n");
        fprintf(stderr,"\t\t  size, alpha and M are unchanged.\n");
        fprintf(stderr,"\t  -f filename     ( Read command line arguments from file filename.)\n");
        fprintf(stderr," See msdoc.pdf for explanation of these parameters.\n");

        exit(1) ;
      }

      // Sort demographic event??
      void addtoelist( struct devent *pt, struct devent *elist ) {
        struct devent *plast, *pevent, *ptemp  ;

        pevent = elist ;
        while (  (pevent != NULL ) && ( pevent->time <= pt->time ) )  {
          plast = pevent ;
          pevent = pevent->nextde ;
        }
        ptemp = plast->nextde ;
        plast->nextde = pt ;
        pt->nextde = ptemp ;
      }

      void
      free_eventlist( struct devent *pt, int npop )
      {
        struct devent *next ;
        int pop ;

        while( pt != NULL){
          next = pt->nextde ;
          if( pt->detype == 'a' ) {
            for( pop = 0; pop < npop; pop++) free( (pt->mat)[pop] );
            free( pt->mat );
          }
          free(pt);
          pt = next ;
        }
      }


      /************ make_gametes.c  *******************************************
      *
      *
      *****************************************************************************/

      #define STATE1 '1'
      #define STATE2 '0'


      // tt : total tree length
      // newsites : number of mutations
      int make_gametes ( int nsam, int mfreq, struct node *ptree, double tt, int newsites, int ns, char **list ){
        int  tip, j,  node ;
        int pickb(int nsam, struct node *ptree, double tt),
        pickbmf(int nsam, int mfreq, struct node *ptree, double tt) ;

        for (  j=ns; j < ns+newsites;  j++ ) {
          if ( mfreq == 1 ){ // Default
            node = pickb ( nsam, ptree, tt ) ;
          }
          else {
            node = pickbmf ( nsam, mfreq, ptree, tt) ;
          }

          for ( tip=0; tip < nsam; tip++ ) {
            if ( tdesn (ptree, tip, node) ){
              list[tip][j] = STATE1 ;
            }
            else {
              list[tip][j] = STATE2 ;
            }
          }
        }
      }

      // N E W ========
      int make_gametes2 ( int nsam, int mfreq, struct node *ptree, double tt, int newsites, int ns, char **list, int **full_list, double *time_mut, int *node_mut, int seg, int start, int end, int *dive ){
        int  tip, j,  node ;
        int pickb2(int nsam, struct node *ptree, double tt, int j, double *time_mut ) ;

        //
        int i, len, t, it, m, count, im, flag ;
        double ran1(), tlim ;
        int len2, start2, end2, flag2 ;
        //


        // Following, I don't care "mfreq" at all
        // Divergence is expressed by non-multiple mutation in a site, while a normal mutation can accept it


        for (  j=ns; j < ns+newsites;  j++ ) {// j is # mutation along the segment
          node_mut[j] = pickb2 ( nsam, ptree, tt, j, time_mut ) ;
        }

        len = end - start + 1 ;

        count = 0;

        //printf("\n------------------------------------\n");
        for ( t = 0; t < *nch; t++ ) {
          //printf("in=%d, out=%d\n", tmp[t][nsam+2], tmp[t][nsam+3]) ;
          if (
            ( (start <= tmp[t][nsam+2])&&(tmp[t][nsam+2] <= end) )
            || ( (start <= tmp[t][nsam+3])&&(tmp[t][nsam+3] <= end) )
            || ( (start >= tmp[t][nsam+2])&&(tmp[t][nsam+3] >= end) )
          ) {
            count++ ;
            //printf ("t = %d\n",t) ;
          }
        }


        //printf("start = %d, end = %d, count = %d\n", start, end, count);
        //printf("newsites = %d (pre + last = newsites, when count > 0)\n", newsites);



        if (newsites>0){
          /*
          printf("node_mut[i]: ");
          for (i = ns; i < ns+newsites; i++){
          printf("\t%d", node_mut[i]);
        }
        printf("\n");
        printf("time_mut[i]: ");
        for (i = ns; i < ns+newsites; i++){
        printf("\t%.2f", time_mut[i]);
      }
      printf("\n\n");
      */
    }

    int k = 0 ;

    double t_count = 0.0 ;

    if (count > 0){
      //printf("%d\n", newsites) ;
      for ( t = *nch - 1; t >= 0; t-- ) {
        flag2 = 0 ;
        if ( (tmp[t][nsam+2] <= start) && (end <= tmp[t][nsam+3]) )  {
          start2 = start ;
          end2 = end ;
          len2 = len ;
          flag2 = 1 ;
        } else if ( (tmp[t][nsam+2] <= start) && (end > tmp[t][nsam+3]) && ( start <= tmp[t][nsam+3] ) ) { // BUG fixed !!
          start2 = start ;
          end2 = tmp[t][nsam+3] ;
          len2 = end2 - start2 + 1 ;
          flag2 = 2 ;
        } else if ( (tmp[t][nsam+2] > start) && (end <= tmp[t][nsam+3]) && ( end >= tmp[t][nsam+2] )  ) { // BUG fixed !!
          start2 = tmp[t][nsam+2] ;
          end2 = end ;
          len2 = end2 - start2 + 1 ;
          flag2 = 3 ;
        } else if ( (tmp[t][nsam+2] > start) && (end > tmp[t][nsam+3]) ) {
          start2 = tmp[t][nsam+2] ;
          end2 = tmp[t][nsam+3] ;
          len2 = end2 - start2 + 1 ;
          flag2 = 4 ;
        } else {
          flag2 = 0 ;
          len2 = 0 ;
          //break ;
        }
        //if ( flag2 > 0 ) printf( "t = %d, in = %d, out = %d, len2 = %d, flag2 = %d\n", t, tmp[t][nsam+2], tmp[t][nsam+3], len2, flag2 ) ;

        if ( len2 > 0 ){

          //printf( "t = %d, in = %d, out = %d, len2 = %d, flag2 = %d\n", t, tmp[t][nsam+2], tmp[t][nsam+3], len2, flag2 ) ;
          //printf("start=%d, end=%d, start1=%d, end1=%d\n", start, end, tmp[t][0], tmp[t][1]);


          if ( t_count == 0 ) {
            tlim = 1000000.0 ;
          } else {
            tlim = t_count ;
          }
          t_count = time_out[t] ;
          //printf ( "t_count = %.2f, tlim = %.2f\n", t_count, tlim ) ;


          for ( j=ns; j < ns+newsites;  j++ ){ // mutation occurs "newsites" times
          // i: location of mutation along sites
          //printf("time_out[%d] = %.2f, time_mut[%d] = %.2f\n", t, time_out[t], j, time_mut[j]);
          if ( (time_out[t] < time_mut[j]) && (time_mut[j] < tlim) ) {
            i = start + (int)( len * ran1() ) ; // if a same pair of "i" is drawn, reverse mutation occurs
            //printf("pre_i = %d\n",i);
            k++ ;
            for ( tip=0; tip < nsam; tip++ ) {
              if ( tdesn (ptree, tip, node_mut[j]) ){
                full_list[tip][i] = (full_list[tip][i] == 0) ? 1 : 0 ;//printf("mut\n");
              }
              else {
                full_list[tip][i] = (full_list[tip][i] == 0) ? 0 : 1 ;
              }
            }
          }//if
        }//for

        ////////// NEW way start (2014 10/16)
        for ( i = start2; i < start2 + len2; i++ ) {
          if ( ran1() < divergence[t] ) {
            for ( tip = 0; tip < nsam; tip++ ) {
              if ( tmp[t][tip+2] ){
                full_list[tip][i] = (full_list[tip][i] == 0) ? 1 : 0 ;//printf("between\n");
              }
            }
          }
        }
        ////////// NEW way end

        //printf( "t = %d, start2 = %d, len2 = %d\n", t, start2, len2 ) ;

        /////////////////////////////////////// OLD way start
        /*
        // confering a divergence to the segment
        for( i = 0; i < pars.cp.nsites; i++ ){
        dive[i] = -1 ;
      }
      m = len2 + 1;
      while (m > len2) m = poisso ( len2 * divergence[t] ) ; //

      //while (m > len2) m = poisso ( len2 * 0.3 ) ;//example, d=0.3
      //printf("m = %d\n", m);
      //printf("m = %d, start2 = %d, end2= %d, len2 = %d\n", m, start2, end2, len2);
      for ( j = 0; j < m; j++ ){ //howmany??
      flag = 1 ;
      while ( flag == 1 ) {
      flag = 0 ;
      i = start2 + (int)( len2 * ran1() ) ;// later, multiple mutations are not accepted!!!!!!!
      //printf("while\n");





      for ( im = 0; im < m; im++ ) {
      if ( dive[im] == i ) flag = 1 ;
    }

  }
  dive[j] = i;
  //printf("div_i = %d\n", i);

  for ( tip = 0; tip < nsam; tip++ ) {
  if ( tmp[t][tip+2] ){
  full_list[tip][i] = (full_list[tip][i] == 0) ? 1 : 0 ;//printf("between\n");
}
}
}
*/
/////////////////////////////////////// OLD way end















//for ( j = 0; j < m; j++ ){
//printf("dive[%d]=%d\n", j, dive[j]);
//}

count-- ;
//printf("count = %d\n", count);
if (count == 0) {
  // final normal mutations
  //printf("last_mut\n");
  for ( j=ns; j < ns+newsites;  j++ ){ // mutation occurs "newsites" times
  // i: location of mutation along sites
  //printf("time_out[%d] = %.2f, time_mut[%d] = %.2f\n", t, time_out[t], j, time_mut[j]);
  if ( time_mut[j] < time_out[t] ) {//printf(">_<\n");
  i = start + (int)( len * ran1() ) ; // if a same pair of "i" is drawn, reverse mutation occurs
  //printf("las_i = %d\n",i);
  k++ ;
  for ( tip=0; tip < nsam; tip++ ) {
    if ( tdesn (ptree, tip, node_mut[j]) ){
      full_list[tip][i] = (full_list[tip][i] == 0) ? 1 : 0 ;
    }
    else {
      full_list[tip][i] = (full_list[tip][i] == 0) ? 0 : 1 ;
    }
  }
}//if
}//for
break ;
}// if


}//if
}//for ( t = *nch - 1; t >= 0; t-- )
}//if (count > 0)
else{
  //printf("%d\n", newsites) ;
  for ( j=ns; j < ns+newsites;  j++ ){ // mutation occurs "newsites" times
  // i: location of mutation along sites
  i = start + (int)( len * ran1() ) ; // if a same pair of "i" is drawn, reverse mutation occurs
  //printf("i = %d\n",i);
  //printf("%d\t",i);
  k++ ;
  for ( tip=0; tip < nsam; tip++ ) {
    if ( tdesn (ptree, tip, node_mut[j]) ){
      full_list[tip][i] = (full_list[tip][i] == 0) ? 1 : 0 ;
    }
    else {
      full_list[tip][i] = (full_list[tip][i] == 0) ? 0 : 1 ;
    }
  }
}
}
if (k!=newsites) {
  printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~") ;
  printf("k=%d,newsites=%d\n",k,newsites);
  exit(0) ;
}




}
// N E W =====end




/***  ttime.c : Returns the total time in the tree, *ptree, with nsam tips. **/

double ttime ( ptree, nsam )
struct node *ptree ;
int nsam ;
{
  double t ;
  int i ;

  t = (ptree + 2*nsam-2) -> time ;
  //printf ("%f\n", t) ;
  for ( i=nsam; i < 2*nsam-1 ; i++ ){
    t += (ptree + i)-> time ;
    //printf ("  time = %f\n", (ptree + i)-> time) ;
  }

  //for ( i=0; i < 2*nsam-1; i++ ){
  //	printf ("%d, t = %f, abv = %d, ndes = %d\n", i, (ptree + i)-> time, (ptree + i)-> abv, (ptree + i)-> ndes) ;
  //}
  //printf ("\n") ;

  return (t) ;
}


double
ttimemf( ptree, nsam, mfreq)
struct node *ptree;
int nsam, mfreq;
{
  double t;
  int i;

  t = 0. ;
  for( i=0;  i< 2*nsam-2  ; i++)
  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
  t += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
  return(t);
}


void
prtree( ptree, nsam)
struct node *ptree;
int nsam;
{
  double t;
  int i, *descl, *descr ;
  void parens( struct node *ptree, int *descl, int *descr, int noden );

  descl = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
  descr = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
  for( i=0; i<2*nsam-1; i++) descl[i] = descr[i] = -1 ;
  for( i = 0; i< 2*nsam-2; i++){
    if( descl[ (ptree+i)->abv ] == -1 ) descl[(ptree+i)->abv] = i ;
    else descr[ (ptree+i)->abv] = i ;
  }
  parens( ptree, descl, descr, 2*nsam-2);
  free( descl ) ;
  free( descr ) ;
}

void
parens( struct node *ptree, int *descl, int *descr,  int noden)
{
  double time ;

  if( descl[noden] == -1 ) {
    printf("%d:%5.3lf", noden+1, (ptree+ ((ptree+noden)->abv))->time );
  }
  else{
    printf("(");
    parens( ptree, descl,descr, descl[noden] ) ;
    printf(",");
    parens(ptree, descl, descr, descr[noden] ) ;
    if( (ptree+noden)->abv == 0 ) printf(");\n");
    else {
      time = (ptree + (ptree+noden)->abv )->time - (ptree+noden)->time ;
      printf("):%5.3lf", time );
    }
  }
}

/***  pickb : returns a random branch from the tree. The probability of picking
a particular branch is proportional to its duration. tt is total
time in tree.   ****/
int pickb ( nsam, ptree, tt )
int nsam ;
struct node *ptree ;
double tt ;
{
  double x, y, ran1() ;
  int i ;

  x = ran1()*tt ;
  for( i=0, y=0; i < 2*nsam-2 ; i++) {
    y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
    if ( y >= x ) return ( i ) ;
  }
  return ( 2*nsam - 3  );  /* changed 4 Feb 2010 */
}

// N E W ========
int pickb2 ( nsam, ptree, tt, j, time_mut )
int nsam ;
struct node *ptree ;
double tt ;
///
int j ;
double *time_mut ;
///
{
  double x, y, ran1() ;
  int i ;

  x = ran1()*tt ;
  for( i=0, y=0; i < 2*nsam-2 ; i++) {
    y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
    if ( y >= x ) {
      ///
      time_mut[j] = (ptree+i)->time + (y - x) ;
      ///
      return ( i ) ;
    }
  }
  time_mut[j] = 0;
  return ( 2*nsam - 3  );  /* changed 4 Feb 2010 */
}
// N E W =====end




int
pickbmf(nsam, mfreq, ptree, tt )
int nsam, mfreq;
struct node *ptree;
double tt;
{
  double x, y, ran1();
  int i, lastbranch = 0 ;

  x = ran1()*tt;
  for( i=0, y=0; i < 2*nsam-2 ; i++) {
    if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ){
      y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
      lastbranch = i ;    /* changed 4 Feb 2010 */
    }
    if( y >= x ) return( i ) ;
  }
  return( lastbranch );   /*  changed 4 Feb 2010 */
}

/****  tdesn : returns 1 if tip is a descendant of node in *ptree, otherwise 0. **/
int tdesn (ptree, tip, node )
struct node *ptree ;
int tip, node ;
{
  int k ;
  for ( k=tip ; k < node ; k = (ptree+k)->abv ) ;
  if ( k == node ) return (1) ;
  else return (0) ;
}

// N E W ========
int tdesn2 (ptree, tip, node )
struct node *ptree ;
int tip, node ;
{
  int k ;
  ///
  //printf ( "call tdesn2, tip = %d, node = %d\n", tip, node );
  for ( k=tip ; k < node ; k = (ptree+k)->abv )
  {
    //printf ( "1: (ptree+%d)->abv = %d\n", k, (ptree+k)->abv );
    if((ptree+k)->abv == 0){
      return (0) ;
    }
  }
  ///
  k = 0 ;
  for ( k=tip ; k < node ; k = (ptree+k)->abv ) {
    //printf ( "2: (ptree+%d)->abv = %d\n", k, (ptree+k)->abv );
  }

  if ( k == node ) {
    return (1) ;
  } else {
    return (0) ;
  }

}
// N E W =====end


/* pick2()  */

int
pick2(n,i,j)
int n, *i, *j;
{
  double ran1();

  *i = n * ran1() ;
  while( ( *j = n * ran1() ) == *i )
  ;
  return(0) ;
}

/**** ordran.c  ***/

void ordran(int n, double pbuf[])// modified 4/21,2014
{
  ranvec(n,pbuf);
  order(n,pbuf);
  return;
}


int
mnmial(n,nclass,p,rv)
int n, nclass, rv[];
double p[];
{
  double ran1();
  double x, s;
  int i, j;

  for(i=0; i<nclass; i++) rv[i]=0;
  for(i=0; i<n ; i++) {
    x = ran1();
    j=0;
    s = p[0];
    while( (x > s) && ( j<(nclass-1) ) )  s += p[++j];
    rv[j]++;
  }
  return(j);
}

int
order(n,pbuf)
int n;
double pbuf[];
{
  int gap, i, j;
  double temp;

  for( gap= n/2; gap>0; gap /= 2)
  for( i=gap; i<n; i++)
  for( j=i-gap; j>=0 && pbuf[j]>pbuf[j+gap]; j -=gap) {
    temp = pbuf[j];
    pbuf[j] = pbuf[j+gap];
    pbuf[j+gap] = temp;
  }
}


void ranvec(int n, double pbuf[])// modified 4/21,2014
{
  int i;
  double ran1();

  for(i=0; i<n; i++)
  pbuf[i] = ran1();

  return;
}



int poisso (u)
double u ;
{
  double  cump, ru, ran1(), p, gasdev(double, double) ;
  int i=1;

  if( u > 30. ){
    i =  (int)(0.5 + gasdev(u,u)) ;
    if( i < 0 ) return( 0 ) ;
    else return( i ) ;
  }

  ru = ran1();
  p = exp(-u);
  if( ru < p) return(0);
  cump = p;

  while( ru > ( cump += (p *= u/i ) ) )
  i++;
  return(i);
}


/* a slight modification of crecipes version */

double gasdev (m,v)
double m, v;
{
  static int iset=0;
  static float gset;
  float fac,r,v1,v2;
  double ran1();

  if  (iset == 0) {
    do {
      v1=2.0*ran1()-1.0;
      v2=2.0*ran1()-1.0;
      r=v1*v1+v2*v2;
    } while (r >= 1.0);
    fac=sqrt(-2.0*log(r)/r);
    gset= v1*fac;
    iset=1;
    return( m + sqrt(v)*v2*fac);
  } else {
    iset=0;
    return( m + sqrt(v)*gset ) ;
  }
}
