/**********  segtre_mig.c **********************************
*
*	This subroutine uses a Monte Carlo algorithm described in
*	Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201, to produce
*	a history of a random sample of gametes under a neutral
*	Wright-Fisher model with recombination and geographic structure.
*	Input parameters
*	are the sample size (nsam), the number of sites between
*	which recombination can occur (nsites), and the recombination
*	rate between the ends of the gametes (r). The function returns
*	nsegs, the number of segments the gametes were broken into
*	in tracing back the history of the gametes.  The histories of
*	these segments are passed back to the calling function in the
*	array of structures seglst[]. An element of this array,  seglst[i],
* 	consists of three parts: (1) beg, the starting point of
*	of segment i, (2) ptree, which points to the first node of the
*	tree representing the history of the segment, (3) next, which
*	is the index number of the next segment.
*	     A tree is a contiguous set of 2*nsam nodes. The first nsam
*	nodes are the tips of the tree, the sampled gametes.  The other
*	nodes are the nodes ancestral to the sampled gametes. Each node
*	consists of an "abv" which is the number of the node ancestral to
*	that node, an "ndes", which is not used or assigned to in this routine,
*	and a "time", which is the time (in units of 4N generations) of the
*	node. For the tips, time equals zero.
*	Returns a pointer to an array of segments, seglst.

**************************************************************************/

/*--------------------------------------------------------
141007

--------------------------------------------------------*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ms.h"
#define NL putchar('\n')
#define size_t unsigned

#define MIN(x, y) ( (x)<(y) ? (x) : (y) )

#define ERROR(message) fprintf(stderr,message),NL,exit(1)

#define SEGINC 80

extern int flag ;

int nchrom, begs, nsegs ;
long nlinks ;
static int *nnodes = NULL ;
double t, cleft , pc, lnpc ;

// N E W ========
double cleft2, pc2, lnpc2, call ;
int between ; // if between rec occurs, flag = 1 // not used
long nlinks2 ;
// N E W =====end

static unsigned seglimit = SEGINC ;
static unsigned maxchr ;

struct seg {
	int beg ;
	int end ;
	int desc ;
} ;

struct chromo {
	int nseg ;
	int pop ;
	struct seg *pseg ;
} ;

static struct chromo *chrom = NULL ;

struct node {
	int abv ;
	int ndes ;
	float time ;
} *ptree1, *ptree2 ;

struct segl {
	int beg ;
	struct node *ptree ;
	int next ;
}  ;

static struct segl *seglst = NULL ;

void printon () ; //st

struct segl *segtre_mig ( struct c_params *cp, int *pnsegs, int **tmp, double *time_out, int *nch, double **two_dim,
	double *divergence, **tractl ) { //st add tractl

		int i, j, k, seg, dec, pop, pop2, c1, c2, ind, rchrom  ;
		int migrant, source_pop, *config ;
		double  ran1(), sum, x, tcoal, ttemp, rft, clefta,  tmin, p  ;
		double prec, cin,  prect, nnm1, nnm0, mig, ran, coal_prob, prob, rdum , arg ;
		char c, event ;
		int re(), xover(), cinr(), cleftr(), eflag, cpop, ic  ;
		int nsam, npop, nsites, *inconfig ;
		double r,  f, rf,  track_len, *nrec, *npast, *tpast, **migm ;
		double *size, *alphag, *tlast ;
		struct devent *nextevent ;

		// N E W ========
		double f2, rf2, rft2 ;
		double prect2, clefta2, cin2, calla ;
		int cinr2(), cleftr2(), xover2() ;
		void callar() ;
		// N E W =====end

		// N E W 2 ========
		double sum_w2 ;
		// N E W 2 ========end
		int sum_bag = 0 ;

		nsam = cp->nsam ;
		npop = cp->npop ;
		nsites = cp->nsites ;
		inconfig = cp->config ;
		r = cp->r ;
		f = cp->f ;
		track_len = cp->track_len ;

		//printf("r=%f, nsites=%d\n", r, nsites);
		// N E W ========
		f2 = cp->f2 ;
		// N E W =====end

		migm = (double **)malloc( (unsigned)npop*sizeof(double *) ) ;
		for ( i=0; i < npop; i++ ) {
			migm[i] = (double *)malloc( (unsigned)npop*sizeof( double) ) ;
			for ( j=0; j < npop; j++ ) migm[i][j] = (cp->mig_mat)[i][j] ;
		}
		nextevent = cp->deventlist ;

		/* Initialization */
		if ( chrom == NULL ) {
			maxchr = nsam + 20 ;
			chrom = (struct chromo *)malloc( (unsigned)( maxchr*sizeof( struct chromo) )) ;
			if ( chrom == NULL ) perror( "malloc error. segtre") ;
		}
		if ( nnodes == NULL ){
			nnodes = (int*) malloc((unsigned)(seglimit*sizeof(int))) ;
			if ( nnodes == NULL ) perror("malloc error. segtre_mig") ;
		}
		if ( seglst == NULL ) {
			seglst = (struct segl *)malloc((unsigned)(seglimit*sizeof(struct segl)) ) ;
			if ( seglst == NULL ) perror("malloc error. segtre_mig.c 2") ;
		}

		config = (int *)malloc( (unsigned) ((npop+1)*sizeof(int) )) ;
		if ( config == NULL ) perror("malloc error. segtre.");
		size = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
		if ( size == NULL ) perror("malloc error. segtre.");
		alphag = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
		if ( alphag == NULL ) perror("malloc error. segtre.");
		tlast = (double *)malloc( (unsigned) ((npop)*sizeof(double) )) ;
		if( alphag == NULL ) perror("malloc error. segtre.");

		for ( pop=0; pop < npop; pop++ ) {
			config[pop] = inconfig[pop] ; // sample size
			size[pop] = (cp->size)[pop] ; // pop size
			alphag[pop] = (cp->alphag)[pop] ;
			tlast[pop] = 0.0 ;
		}
		for ( pop=ind=0; pop < npop; pop++ ){
			for ( j=0; j < inconfig[pop]; j++,ind++ ) {

				chrom[ind].nseg = 1 ;
				if ( !(chrom[ind].pseg = (struct seg*)malloc((unsigned)sizeof(struct seg)) ))
				ERROR ("calloc error. se1");

				(chrom[ind].pseg)->beg = 0 ;
				(chrom[ind].pseg)->end = nsites-1 ;
				(chrom[ind].pseg)->desc = ind ;
				chrom[ind].pop = pop ;
			}
		}
		seglst[0].beg = 0 ;
		if( !(seglst[0].ptree = (struct node *)calloc((unsigned)(2*nsam),sizeof(struct node)) ))
		perror ("calloc error. se2") ;

		nnodes[0] = nsam - 1 ;
		nchrom = nsam ;
		nlinks = ((long)(nsam))*(nsites-1) ;
		nsegs = 1 ;
		t = 0.;
		r /= (nsites-1) ;
		if ( f > 0.0 ) 	pc = (track_len -1.0)/track_len ;
		else pc = 1.0 ;
		lnpc = log ( pc ) ;
		cleft = nsam* ( 1.0 - pow( pc, (double)(nsites-1) ) ) ;
		if ( r > 0.0 ) rf = r*f ;
		else rf = f /(nsites-1) ;
		rft = rf*track_len ;

		nlinks2 = nlinks ;
		*nch = 0 ;

		if ( r > 0.0 ) rf2 = r * f2 ;
		else rf2 = f2 /(nsites-1) ; // per site


		// Define the probability of cleft conditional on the occurence of cleft2 in a focsing region
		//st
		sum_w2 = 0.0 ;
		for ( j=0; j < cp->n_two_dim; j++ ){
			sum_w2 += (two_dim[j][0])*two_dim[j][3] ;
		}

		/*printf ("1 = %f\n", sum_w2) ;
		sum_w2 = 0.0 ;
		//for ( i=1; i <= nsites-1; i++ ){
		for ( i=1; i <= 100000; i++ ){
		for ( j=0; j < cp->n_two_dim; j++ ){
		if ( two_dim[j][0] >= i ){
		sum_w2 += two_dim[j][3] ;
	}
}
}
printf ("2 = %f\n", sum_w2) ;
exit (0) ;*/

while ( nchrom > 1 ) {
	//printf("------------- start subroutine -------------\n");
	//printf("\n<start subroutine>,  nchrom = %d, nlinks = %d, nlinks2 = %d\n", nchrom, nlinks, nlinks2);

	if ( nlinks < 0 ) printf ( "nlinks2 is negative !!\n" ) ;;

	prec = nlinks*r;
	cin = nlinks*rf ;
	clefta = cleft*rft ;
	prect = prec + cin + clefta ;

	// N E W 2 ========

	// --- sum_w2 does not deppend on the effective regions, just on the number of chrom, because "call" (overall region is covered by foregin DNA) is likely to occur in shorter region, thus, due to this canceling out, the cleft2/nchrom (including "call") is constant ---

	cin2 = nlinks2 * rf2 ; //cin2=0

	cleft2 = nchrom * sum_w2 ;
	clefta2 = cleft2 * rf2 ;

	prect2 = cin2 + clefta2 ;

	//printf("cin = %.2f, clefta = %.2f\n", cin, clefta);
	//printf("cin2 = %.2f, clefta2 = %.2f\n", cin2, clefta2);
	// N E W 2 =====end

	mig = 0.0;
	for ( i=0; i<npop; i++ ) mig += config[i]*migm[i][i] ;
	if ( (npop > 1) && ( mig == 0.0) && (nextevent == NULL) ) {
		i = 0 ;
		for ( j=0; j < npop; j++ )
		if ( config[j] > 0 ) i++;
		if ( i > 1 ) {
			fprintf(stderr," Infinite coalescent time. No migration.\n");
			exit(1);
		}
	}
	eflag = 0 ;

	if ( prect > 0.0 ) {      /* cross-over or gene conversion */
		while( (rdum = ran1() )  == 0.0 ) ;
		ttemp = -log( rdum)/prect ;
		if( (eflag == 0) || (ttemp < tmin ) ){
			tmin = ttemp;
			event = 'r' ;
			eflag = 1;
		}
	}

	// N E W ========
	if ( prect2 > 0.0 ) { // recombination BETWEEN species
		while( (rdum = ran1() )  == 0.0 ) ;
		ttemp = -log( rdum)/prect2 ;
		if( (eflag == 0) || (ttemp < tmin ) ){
			tmin = ttemp;
			event = 'b' ;
			eflag = 1;
		}

	}
	// N E W =====end

	if ( mig > 0.0 ) {         /* migration   */
		while( (rdum = ran1() ) == 0.0 ) ;
		ttemp = -log( rdum)/mig ;
		if( (eflag == 0) || (ttemp < tmin ) ){
			tmin = ttemp;
			event = 'm' ;
			eflag = 1 ;
		}
	}
	/*if ( theta > 0.0 ){  /* mutation
	while( (rdum = ran1() )  == 0.0 ) ;
	ttemp = -log( rdum)/theta/2.0 ;
	if ( (eflag == 0) || (ttemp < tmin ) ){
	tmin = ttemp;
	event = 'u' ;
	eflag = 1 ;
}
}*/

for ( pop=0; pop < npop ; pop++ ) {     /* coalescent */
	coal_prob = ((double)config[pop])*(config[pop]-1.) ;
	if ( coal_prob > 0.0 ) {
		while ( ( rdum = ran1() ) == .0 ) ;
		if ( alphag[pop] == 0 ){ // constant population size
			ttemp = -log ( rdum )*size[pop] /coal_prob ;
			if ( (eflag == 0) || (ttemp < tmin) ){
				tmin = ttemp ;
				event = 'c' ;
				eflag = 1 ;
				cpop = pop ;
			}
		}
		else { // population growth
			arg  = 1. - alphag[pop]*size[pop]*exp(-alphag[pop]*(t - tlast[pop] ) )* log(rdum) / coal_prob ;
			if ( arg > 0.0 ) {                          /*if arg <= 0,  no coalescent within interval */
				ttemp = log( arg ) / alphag[pop]  ;
				if ( (eflag == 0) || (ttemp < tmin) ){
					tmin = ttemp ;
					event = 'c' ;
					eflag = 1 ;
					cpop = pop ;
				}
			}
		}
	}
}

if ( (eflag == 0) && ( nextevent == NULL) ) {
	fprintf(stderr,
		" infinite time to next event. Negative growth rate in last time interval or non-communicating subpops.\n");
		exit( 0);
	}

	if( ( ( eflag == 0) && (nextevent != NULL))|| ( (nextevent != NULL) &&  ( (t+tmin) >=  nextevent->time)) ) { // demography
		t = nextevent->time ;
		switch ( nextevent->detype ) {
			case 'N' :
			for(pop =0; pop <npop; pop++){
				size[pop]= nextevent->paramv ;
				alphag[pop] = 0.0 ;
			}
			nextevent = nextevent->nextde ;
			break;
			case 'n' :
			size[nextevent->popi]= nextevent->paramv ;
			alphag[nextevent->popi] = 0.0 ;
			nextevent = nextevent->nextde ;
			break;
			case 'G' :
			for(pop =0; pop <npop; pop++){
				size[pop] = size[pop]*exp( -alphag[pop]*(t - tlast[pop]) ) ;
				alphag[pop]= nextevent->paramv ;
				tlast[pop] = t ;
			}
			nextevent = nextevent->nextde ;
			break;
			case 'g' :
			pop = nextevent->popi ;
			size[pop] = size[pop]*exp( - alphag[pop]*(t-tlast[pop]) ) ;
			alphag[pop]= nextevent->paramv ;
			tlast[pop] = t ;
			nextevent = nextevent->nextde ;
			break;
			case 'M' :
			for(pop =0; pop <npop; pop++)
			for( pop2 = 0; pop2 <npop; pop2++) migm[pop][pop2] = (nextevent->paramv) /(npop-1.0) ;
			for( pop = 0; pop <npop; pop++)
			migm[pop][pop]= nextevent->paramv ;
			nextevent = nextevent->nextde ;
			break;
			case 'a' :
			for(pop =0; pop <npop; pop++)
			for( pop2 = 0; pop2 <npop; pop2++) migm[pop][pop2] = (nextevent->mat)[pop][pop2]  ;
			nextevent = nextevent->nextde ;
			break;
			case 'm' :
			i = nextevent->popi ;
			j = nextevent->popj ;
			migm[i][i] += nextevent->paramv - migm[i][j];
			migm[i][j]= nextevent->paramv ;
			nextevent = nextevent->nextde ;
			break;
			case 'j' :         /* merge pop i into pop j  (join) */
			i = nextevent->popi ;
			j = nextevent->popj ;
			config[j] += config[i] ;
			config[i] = 0 ;
			for( ic = 0; ic<nchrom; ic++) if( chrom[ic].pop == i ) chrom[ic].pop = j ;
			/*  the following was added 19 May 2007 */
			for( k=0; k < npop; k++){
				if( k != i) {
					migm[k][k] -= migm[k][i] ;
					migm[k][i] = 0. ;
				}
			}
			/* end addition */
			nextevent = nextevent->nextde ;
			break;
			case 's' :         /*split  pop i into two;p is the proportion from pop i, and 1-p from pop n+1  */
			i = nextevent->popi ;
			p = nextevent->paramv ;
			npop++;
			config = (int *)realloc( config, (unsigned)(npop*sizeof( int) ));
			size = (double *)realloc(size, (unsigned)(npop*sizeof(double) ));
			alphag = (double *)realloc(alphag, (unsigned)(npop*sizeof(double) ));
			tlast = (double *)realloc(tlast,(unsigned)(npop*sizeof(double) ) ) ;
			tlast[npop-1] = t ;
			size[npop-1] = 1.0 ;
			alphag[npop-1] = 0.0 ;
			migm = (double **)realloc(migm, (unsigned)(npop*sizeof( double *)));
			for( j=0; j< npop-1; j++)
			migm[j] = (double *)realloc(migm[j],(unsigned)(npop*sizeof(double)));
			migm[npop-1] = (double *)malloc( (unsigned)(npop*sizeof( double) ) ) ;
			for( j=0; j<npop; j++) migm[npop-1][j] = migm[j][npop-1] = 0.0 ;
			config[npop-1] = 0 ;
			config[i] = 0 ;
			for( ic = 0; ic<nchrom; ic++){
				if( chrom[ic].pop == i ) {
					if( ran1() < p ) config[i]++;
					else {
						chrom[ic].pop = npop-1 ;
						config[npop-1]++;
					}
				}
			}
			nextevent = nextevent->nextde ;
			break;
		}
	}
	else { // mutational mechanisms
		t += tmin ;
		if ( event == 'r' ) {
			//printf("event = within\n");
			between = 0 ;
			if( (ran = ran1()) < ( prec / prect ) ){ /*recombination*/
				rchrom = re(nsam );
				config[ chrom[rchrom].pop ] += 1 ;
			}
			else if( ran < (prec + clefta)/(prect) ){    /*  cleft event */
				//printf("===========  cleftr WITHIN species  ===========\n");
				//printf( "t = %.3f\n", t ) ;
				rchrom = cleftr(nsam );
				config[ chrom[rchrom].pop ] += 1 ;
				//printf("===========  cleftr WITHIN species  ===========\n");
				//printf("nchrom = %d\n",nchrom);
				//printf("\n");
			}
			else  {         /* cin event */
				//printf("===========  cinr WITHIN species  ===========\n");
				//printf( "t = %.3f\n", t ) ;
				rchrom = cinr(nsam,nsites );
				if ( rchrom >= 0 ) config[ chrom[rchrom].pop ] += 1 ;
				//printf("===========  cinr WITHIN species  ===========\n");
				//printf("nchrom = %d\n",nchrom);
				//printf("\n");

			}
		}

		// N E W ========
		else if ( event == 'b' ) {

			// we do not need to distinguish call event from cleft because they are simply managed in cleftr()

			//printf("event = between\n");
			between = 1 ;
			//printf("*nch=%d\n",*nch);
			if ( (ran = ran1()) < (clefta2/prect2) ) {
				//printf("===========  cleftr2 BETWEEN species  ===========~~~~~~~~~~~~\n");
				//printf( "t = %.3f\n", t ) ;
				rchrom = cleftr2(nsam, nsites, tmp, time_out, nch, divergence, tractl, cp->n_tractl );//st
				//printf("===========  cleftr2 BETWEEN species  ===========\n");
			}
			else  {
				//printf("===========  cinr2 BETWEEN species  ===========~~~~~~~~~~~~\n");
				rchrom = cinr2(nsam, nsites, tmp, time_out, nch, two_dim, cp->n_two_dim, divergence );//st
			}
		}
		// N E W =====end

		else if ( event == 'm' ) {  /* migration event */
			x = mig*ran1();
			sum = 0.0 ;
			for( i=0; i<nchrom; i++) {
				sum += migm[chrom[i].pop][chrom[i].pop] ;
				if( x <sum ) break;
			}
			migrant = i ;
			x = ran1()*migm[chrom[i].pop][chrom[i].pop];
			sum = 0.0;
			for(i=0; i<npop; i++){
				if ( i != chrom[migrant].pop ){
					sum += migm[chrom[migrant].pop][i];
					if( x < sum ) break;
				}
			}
			source_pop = i;
			config[chrom[migrant].pop] -= 1;
			config[source_pop] += 1;
			chrom[migrant].pop = source_pop ;
		}

		else { 								 /* coalescent event */

			//printf("===========  coalescent  ===========\n");
			//printf( "t = %.3f\n", t ) ;
			/* pick the two, c1, c2  */
			pick2_chrom ( cpop, config, &c1, &c2 ) ;  /* c1 and c2 are chrom's to coalesce */
			//printf ("c1 = %d, c2 = %d\n", c1, c2) ; //===//
			//printf( "lnks(c1) = %d\n", links(c1) ) ;
			//coalescent event
			dec = ca (nsam,nsites,c1,c2 ) ;
			config[cpop] -= dec ;
			//printf ("  nnodes = %d\n", nnodes[0]) ;
			//printf("===========  coalescent  ===========\n");
			//printf("nchrom = %d\n",nchrom);
			//printf("\n");
		}
		//printf("tmp[%d][nsam+2]=%d\n\n",*nch,tmp[*nch][nsam+2]) ;
	}
	//Following is for debug mode

	// N E W ========


	/*
	for(i = 0; i < nchrom; i++){
	printf ( "links(%d) = %d, links2(%d) = %d\n", i, links(i), i, links2(i) ) ;
	printf("chrom[%d].nseg = %d\n", i, chrom[i].nseg);
	for(j = 0; j < chrom[i].nseg; j++){
	printf("(chrom[%d].pseg + %d)->beg = %d\n", i, j, (chrom[i].pseg+j)->beg);
	printf("(chrom[%d].pseg + %d)->end = %d\n", i, j, (chrom[i].pseg+j)->end);
}
}
printf("\n");




for(i = 0; i < nsegs; i++){
printf("nnodes[%d] = %d\n", i, nnodes[i]);
}
printf("\n");
struct node *ptree ;
for ( seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++ ) {
ptree = seglst[seg].ptree ;
printf("seg = %d, seglst[%d].beg = %d\n", seg, seg, seglst[seg].beg);
for (i = 0; i < nnodes[seg]+1; i++){
printf ("(ptree + %d)->abv = %d\n", i, (ptree+i)->abv) ;
}
}
printf("nch=%d\n",*nch);
printf("nch, the time, start, end, whether each sample is involved or not\n") ;
for (i = 0; i < *nch; i++){
printf("%d\t%.2f\t", i, time_out[i]);
for (j = 0; j < 2 + nsam + 2; j++){
printf("%d\t", tmp[i][j]);
}
printf("\n");
}
*/


// N E W =====end



}  // while

/*
for(i = 0; i < nchrom; i++){
printf ( "links(%d) = %d, links2(%d) = %d\n", i, links(i), i, links2(i) ) ;
printf("chrom[%d].nseg = %d\n", i, chrom[i].nseg);
for(j = 0; j < chrom[i].nseg; j++){
printf("(chrom[%d].pseg + %d)->beg = %d\n", i, j, (chrom[i].pseg+j)->beg);
printf("(chrom[%d].pseg + %d)->end = %d\n", i, j, (chrom[i].pseg+j)->end);
}
}
printf("\n");
*/


// N E W ========
/*
printf("nsegs = %d, nch = %d\n",nsegs, *nch);

if (*nch>0) {
printf("nch\ttime\tdiv\tstart\tend\tnsam\t...\tin\tout\n") ;
}

for (i = 0; i < *nch; i++){
printf("%d\t%.3f\t%.2f\t", i, time_out[i], divergence[i]);
for (j = 0; j < 4 + nsam; j++){
printf("%d\t", tmp[i][j]);
}
printf("\n");
}*/
// N E W =====end

// N E W 3=======
/*
FILE *fp2 ;
fp2 = fopen( "summary.txt", "a" ) ;

//
fprintf( fp2, "nsegs = %d, nch = %d\n", nsegs, *nch);
if (*nch>0) {
fprintf( fp2, "nch\ttime\tdiv\tstart\tend\tnsam\t...\tin\tout\n" ) ;
}

for (i = 0; i < *nch; i++){
fprintf( fp2, "%d\t%.3f\t%.2f\t", i, time_out[i], divergence[i]);
for (j = 0; j < 4 + nsam; j++){
fprintf( fp2, "%d\t", tmp[i][j]);
if ( (j > 1) && (j < 2 + nsam) ) {
sum_bag += tmp[i][j] ;
}
}
fprintf( fp2, "\n");

if (sum_bag == 0) {
printf ("sum_bag == 0!!!!!!!!!!!!!!!!!!!!!\n") ;
}
sum_bag = 0 ;
}
//
fclose( fp2 ) ;
*/
// N E W 3====end


*pnsegs = nsegs ;
//printf ("nsegs= %d\n", nsegs) ;

free (config);
free (size) ;
free (alphag);
free (tlast);
for ( i=0; i<npop; i++) free ( migm[i] ) ;
free (migm) ;
return ( seglst ) ;

}


void printon (){ //st
	int i, j ;
	for ( i=0; i < nchrom; i++ ){
		printf ("  %dth %d seg\n", i, chrom[i].nseg) ;

		for ( j=0; j < chrom[i].nseg; j++ ){
			printf ("    %dth seg %d - %d\n", j, chrom[i].pseg[j].beg, chrom[i].pseg[j].end) ;
		}
	}


}

/******  recombination subroutine ***************************
Picks a chromosome and splits it in two parts. If the x-over point
is in a new spot, a new segment is added to seglst and a tree set up
for it.   ****/


int re (nsam)
int nsam ;
{
	struct seg *pseg ;
	int  el, lsg, lsgm1,  ic,  is, in, spot;
	double ran1();

	/* First generate a random x-over spot, then locate it as to chrom and seg. */

	spot = nlinks*ran1() + 1.;

	/* get chromosome # (ic)  */

	for( ic=0; ic<nchrom ; ic++) {
		lsg = chrom[ic].nseg ;
		lsgm1 = lsg - 1;
		pseg = chrom[ic].pseg;
		el = ( (pseg+lsgm1)->end ) - (pseg->beg);
		if( spot <= el ) break;
		spot -= el ;
	}
	is = pseg->beg + spot -1;
	xover(nsam, ic, is); // ADD
	return(ic);
}

int
cleftr( int nsam )
{
	struct seg *pseg ;
	int   lsg, lsgm1,  ic,  is, in, spot;
	double ran1(), x, sum, len  ;
	int xover(int, int, int);
	while( (x = cleft*ran1() )== 0.0 )  ;
	sum = 0.0 ;
	ic = -1 ;
	while ( sum < x ) {
		sum +=  1.0 - pow( pc, links(++ic) )  ;
	}
	pseg = chrom[ic].pseg;
	len = links(ic) ;
	is = pseg->beg + floor( 1.0 + log( 1.0 - (1.0- pow( pc, len))*ran1() )/lnpc  ) -1  ;
	xover( nsam, ic, is );
	return( ic) ;
}


// N E W ========
int
cleftr2( int nsam, int nsites, int **tmp, double *time_out, int *nch,
	double *divergence, double **tractl, int n_tractl )// BUG is HERE!!!!!!!!!!!!!!!!!!!
	{

		struct seg *pseg ;
		int   lsg, lsgm1,  ic,  is, in, spot;
		double ran1(), x, sum, len  ;
		int xover(int, int, int);
		///
		int seg, tip, node, i, k, z ;
		double div ;

		ic = floor( nchrom * ran1() ) ;
		pseg = chrom[ic].pseg ;

		//printf( "ic = %d\n", ic ) ;

		x = ran1() ;
		//st
		for ( i=0; i < n_tractl; i++ ){ //st get tract and divergence
			if ( x < tractl[i][3] ){
				break ;
			}
		}
		is = (int)tractl[i][0] ;
		div = tractl[i][1] ;



		//printf( "div = %.2f, ic = %d, z = %d, is = %d, pseg->beg = %d, pseg->end = %d\n", div, ic, z, is, pseg->beg, (pseg+(chrom[ic].nseg-1))->end) ;
		//printf( "is = %d, pseg->beg = %d, pseg->end = %d\n", is, pseg->beg, (pseg+(chrom[ic].nseg-1))->end) ;

		if ( is >= pseg->beg ) { // OK 2014 11/11
			//printf ( "is > pseg->beg\n" ) ;
			//printf( "cleftr2, t = %.2f, ic = %d, is = %d, pseg->beg = %d\n", t, ic, is, pseg->beg ) ;
			if ( is >= ( pseg + (chrom[ic].nseg-1) )->end ) { // OK 2014 11/11
				// finished at out of efficient segments
				//printf ( "finished outside of region\n" ) ;
				//for ( seg = 0, i = 0; i < chrom[ic].nseg; seg = seglst[seg].next, i++ ) {
				for ( seg = 0, i = 0; i < chrom[ic].nseg; i++ ) {

					while ( (pseg + i)->beg != seglst[seg].beg ) {
						seg = seglst[seg].next ;
					}
					//printf ( "seg=%d\n", seg ) ;

					tmp[*nch][0] = (pseg + i)->beg ;
					tmp[*nch][1] = (pseg + i)->end ;
					time_out[*nch] = t ;
					divergence[*nch] = div ;

					//printf("beg=%d, end=%d, div=%.2f\n", tmp[*nch][0], tmp[*nch][1], div) ;
					node = (pseg + i)->desc ;
					for ( tip=0; tip < nsam; tip++ ) {
						//printf("i = %d, seg=%d, node=%d\n", i, seg, node);
						//if ( tdesn2 (seglst[seg].ptree, tip, node) ){
						if ( tdesn2 (seglst[seg].ptree, tip, node) ){
							tmp[*nch][2+tip] = 1 ;
						}
						else {
							tmp[*nch][2+tip] = 0 ;
						}
					}
					tmp[*nch][nsam+2] = (pseg + i)->beg ;
					tmp[*nch][nsam+3] = (pseg + i)->end ;
					//printf("t=%f, tmp[%d][nsam+2]=%d\n",t, *nch,tmp[*nch][nsam+2]) ;
					(*nch)++ ;
				}
			}
			else {
				//printf ( "%d,", is ) ;
				// finished within efficient segments
				//printf ( "finished inside of region\n" ) ;
				//printf ( " ( pseg + (chrom[ic].nseg-1) )->end = %d\n", ( pseg + (chrom[ic].nseg-1) )->end ) ;
				//for ( seg = 0, i = 0; i < chrom[ic].nseg; seg = seglst[seg].next, i++ ) {
				for ( seg = 0, i = 0; i < chrom[ic].nseg; i++ ) {

					while ( (pseg + i)->beg != seglst[seg].beg ) {
						seg = seglst[seg].next ;
					}

					//printf ( "seg=%d\n", seg ) ;

					if ( is < (pseg + i)->beg ) {
						//printf ( "break!\n" ) ;
						break ; //remove the newly arisen tmp due to the steps over the finished-point
					}

					tmp[*nch][0] = (pseg + i)->beg ;
					tmp[*nch][1] = (pseg + i)->end ;
					time_out[*nch] = t ;
					divergence[*nch] = div ;

					//printf("beg=%d, end=%d, div=%.2f\n", tmp[*nch][0], tmp[*nch][1], div) ;
					node = (pseg + i)->desc ;
					for ( tip=0; tip < nsam; tip++ ) {
						//printf("i = %d, seg=%d, node=%d\n", i, seg, node);
						if ( tdesn2 (seglst[seg].ptree, tip, node) ){
							tmp[*nch][2+tip] = 1 ;
						}
						else {
							tmp[*nch][2+tip] = 0 ;
						}
					}
					// finished within focussing segment
					if ( (pseg + i)->beg <= is && is < (pseg + i)->end ) { // OK 2014 11/11
						//printf("t=%f break!!!\n",t);
						//printf("%d %d\n",tmp[*nch][nsam+2], is) ;
						tmp[*nch][nsam+2] = (pseg + i)->beg ;
						tmp[*nch][nsam+3] = is ;
						//printf("%d,",tmp[*nch][nsam+3]-tmp[*nch][nsam+2]+1) ;
						(*nch)++ ;
						break ;
					}
					tmp[*nch][nsam+2] = (pseg + i)->beg ;
					tmp[*nch][nsam+3] = (pseg + i)->end ;
					(*nch)++ ;
				}
			}
		}
		//printf("cleftr2_end ");
		return( ic ) ;
	}
	// N E W =====end

	int
	cinr( int nsam, int nsites )
	{
		struct seg *pseg ;
		int len,  el, lsg, lsgm1,  ic,  is, in, spot, endic ;
		double ran1();
		int xover(), ca() ;

		/* First generate a random x-over spot, then locate it as to chrom and seg. */

		spot = nlinks*ran1() + 1.;

		/* get chromosome # (ic)  */

		for( ic=0; ic<nchrom ; ic++) {
			lsg = chrom[ic].nseg ;
			lsgm1 = lsg - 1;
			pseg = chrom[ic].pseg;
			el = ( (pseg+lsgm1)->end ) - (pseg->beg);
			if( spot <= el ) break;
			spot -= el ;
		}
		is = pseg->beg + spot -1;
		endic = (pseg+lsgm1)->end ;
		xover(nsam, ic, is);

		len = floor( 1.0 + log( ran1() )/lnpc ) ;
		if( is+len >= endic ) return(ic) ;
		if( is+len < (chrom[nchrom-1].pseg)->beg ){
			ca( nsam, nsites, ic, nchrom-1);
			return(-1) ;
		}
		xover( nsam, nchrom-1, is+len ) ;
		ca( nsam,nsites, ic,  nchrom-1);
		return(ic);

	}


	//cleftr2( int nsam, int **tmp, double *time_out, int *nch, double **two_dim )
	// N E W ========
	int cinr2( int nsam, int nsites, int **tmp, double *time_out, int *nch, double **two_dim, int n_two_dim, double *divergence ){ //st
		struct seg *pseg ;
		int len,  el, lsg, lsgm1,  ic,  is, in, spot, endic ;
		double ran1();
		int xover(), ca() ;
		///
		int seg, tip, node, i, j ;
		struct seg *pseg2 ;
		double x, div;
		///
		//printf( "cinr2\n");
		//printf ("%d\n", n_two_dim) ; exit (0) ;

		/* First generate a random x-over spot, then locate it as to chrom and seg. */

		spot = nlinks2 * ran1() + 1.;

		/* get chromosome # (ic)  */

		for( ic=0; ic<nchrom ; ic++) {
			lsg = chrom[ic].nseg ;
			lsgm1 = lsg - 1;
			pseg = chrom[ic].pseg;
			el = (pseg+lsgm1)->end ; // nlinks2 includes empty zone leftside of effective region
			if( spot <= el ) break;
			spot -= el ;
		}
		//printf("spot = %d \n", spot) ;
		// spot >= 1, and is >= 0
		is = spot - 1 ; // "is" is the distance occuring initiation from position zero
		endic = (pseg+lsgm1)->end ;

		//st
		x = ran1() ;
		for ( i=0; i < n_two_dim; i++ ){ //st get tract and divergence
			if ( x < two_dim[i][4] ){
				break ;
			}
		}
		len = (int)two_dim[i][0] ;
		div = two_dim[i][1] ;
		//printf ("x = %f, len = %d, div = %f\n", x, len, div) ; //exit (0) ;

		if ( len == 0 ){
			printf ("Error: tract length must be > 0\n") ;
			exit (1) ;
		}

		//printf( "ic = %d, is = %d, is + len = %d, pseg->beg = %d\n", ic, is, is + len, pseg->beg ) ;
		//printf( "div = %.2f, ic = %d, len = %d, is = %d, endic = %d\n", div, ic, len, is, endic ) ;
		//printf( "ic = %d, len = %d, is = %d, endic = %d\n", ic, len, is, endic ) ;
		//printf ( "is+len = %d, pseg->beg = %d\n", is+len, pseg->beg ) ;

		//printon () ;
		//printf ("x = %f, len = %d, div = %f\n", x, len, div) ;
		//printf ("ic = %d, is = %d\n", ic, is) ;

		if ( is+len >= pseg->beg ) { // ">=" is needed !
		// end of conversion is before the effective region

		// finished outside of efficient segments
		if( is+len >= endic ) { // OK 2014 11/11

			//printf ("***\n") ;

			//printf ( "finished outside of efficient segments! \n" ) ;

			//for ( seg = 0, i = 0; i < chrom[ic].nseg; seg = seglst[seg].next, i++ ) {
			for ( seg = 0, i = 0; i < chrom[ic].nseg; i++ ) { // OK 2014 11/11

				while ( (pseg + i)->beg != seglst[seg].beg ) {
					seg = seglst[seg].next ;
				}
				//printf("(pseg2 + %d)->beg = %d, (pseg2 + %d)->end=%d\n", i, (pseg2 + i)->beg, i, (pseg2 + i)->end);

				if ( is < (pseg + i)->end ) { // OK 2014 11/11
					//printf ( ">_< YES \n" ) ;
					//printf ( "(pseg + i)->end = %d\n", (pseg + i)->end ) ;
					//printf( "cinr2_out, t = %.2f, ic = %d, len = %d, is = %d, endic = %d\n", t, ic, len, is, endic ) ;
					tmp[*nch][0] = (pseg + i)->beg ;
					tmp[*nch][1] = (pseg + i)->end ;
					time_out[*nch] = t ;
					divergence[*nch] = div ;

					//printf("beg=%d, end=%d, div=%.2f\n", tmp[*nch][0], tmp[*nch][1], div) ;
					node = (pseg + i)->desc ;
					for ( tip=0; tip < nsam; tip++ ) {
						//printf("i = %d, seg=%d, node=%d\n", i, seg, node);
						if ( tdesn2 (seglst[seg].ptree, tip, node) ){
							tmp[*nch][2+tip] = 1 ;
						}
						else {
							tmp[*nch][2+tip] = 0 ;
						}
					}
					if ( (pseg + i)->beg <= is && is < (pseg + i)->end  ) { // OK 2014 11/11
						tmp[*nch][nsam+2] = is + 1 ;
					} else {
						tmp[*nch][nsam+2] = (pseg + i)->beg ;
					}
					tmp[*nch][nsam+3] = (pseg + i)->end ;
					//printf("t=%f, tmp[%d][nsam+2]=%d\n",t, *nch,tmp[*nch][nsam+2]) ;
					(*nch)++ ;
				} else {
					//printf ( ">_< NO \n" ) ;
				}
			}
			return(ic) ;
		}

		//printf ( "finished inside of efficient segments! \n" ) ;

		// finished inside of efficient segments
		//for ( seg = 0, i = 0; i < chrom[ic].nseg; seg = seglst[seg].next, i++ ) {
		//printf ("---\n") ;
		for ( seg = 0, i = 0; i < chrom[ic].nseg; i++ ) {

			while ( (pseg + i)->beg != seglst[seg].beg ) {
				seg = seglst[seg].next ;
			}
			if ( (is < (pseg + i)->end) && (is + len >= (pseg + i)->beg) ) { // OK 2014 11/11
				//printf( "cinr2_in, t = %.2f, ic = %d, len = %d, is = %d, endic = %d\n", t, ic, len, is, endic ) ;
				tmp[*nch][0] = (pseg + i)->beg ;
				tmp[*nch][1] = (pseg + i)->end ;
				time_out[*nch] = t ;
				divergence[*nch] = div ;

				//printf("beg=%d, end=%d, div=%.2f\n", tmp[*nch][0], tmp[*nch][1], div) ;
				node = (pseg + i)->desc ;
				for ( tip=0; tip < nsam; tip++ ) {
					//printf("i = %d, seg=%d, node=%d\n", i, seg, node);
					if ( tdesn2 (seglst[seg].ptree, tip, node) ){
						tmp[*nch][2+tip] = 1 ;
					}
					else {
						tmp[*nch][2+tip] = 0 ;
					}
				}
				if ( (pseg + i)->beg <= is && is < (pseg + i)->end  ) { // maybe OK 2014 11/11
					tmp[*nch][nsam+2] = is + 1 ;
				} else {
					tmp[*nch][nsam+2] = (pseg + i)->beg ;
				}

				if ( (pseg + i)->beg <= is+len && is+len < (pseg + i)->end  ) { // maybe OK 2014 11/11
					tmp[*nch][nsam+3] = is + len ;
					//(*nch)++ ;
					//break ;
				} else {
					tmp[*nch][nsam+3] = (pseg + i)->end ;
				}

				(*nch)++ ;
			}
		}



	}
	else {
		//printf ( "karaburi!!\n" ) ;
	}

	return(ic);

}
// N E W =====end

// N E W ========
void
callar( int nsam, int **tmp, double *time_out, int *nch )
{
	struct seg *pseg ;
	int len,  el, lsg, lsgm1,  ic,  is, in, spot, endic ;
	double ran1(), x, sum;
	int xover(), ca() ;
	///
	int seg, tip, node, i, j ;
	struct seg *pseg2 ;
	///
	//printf("call callar\n");
	while( (x = call * ran1() )== 0.0 )  ;
	sum = 0.0 ;
	ic = -1 ;
	//printf("call = %f, x = %f\n", call, x);
	while ( sum < x ) {
		//printf(">_<\n");
		sum +=  pow( pc2, links(++ic) )  ;
		//printf("ic = %d, sum = %f, x = %f\n", ic, sum, x );
	}

	pseg = chrom[ic].pseg;
	for ( seg = 0, i = 0; i < chrom[ic].nseg; seg = seglst[seg].next, i++ ) {
		//for (i=0; i < chrom[nchrom-1].nseg; i++){
		//printf("(pseg2 + %d)->beg = %d, (pseg2 + %d)->end=%d\n", i, (pseg2 + i)->beg, i, (pseg2 + i)->end);
		tmp[*nch][0] = (pseg + i)->beg ;
		tmp[*nch][1] = (pseg + i)->end ;
		time_out[*nch] = t ;

		node = (pseg + i)->desc ;
		for ( tip=0; tip < nsam; tip++ ) {
			//printf("seg=%d, node=%d\n",seg, node);
			if ( tdesn2 (seglst[seg].ptree, tip, node) ){
				tmp[*nch][2+tip] = 1 ;
			}
			else {
				tmp[*nch][2+tip] = 0 ;
			}
		}
		(*nch)++ ;
	}
}
// N E W =====end


int
xover(int nsam,int ic, int is)
{
	//printf("between = %d\n", between);
	struct seg *pseg, *pseg2;
	int i,  lsg, lsgm1, newsg,  jseg, k,  in, spot;
	double ran1(), len ;

	// N E W 3=======
	/*
	FILE *fp1 ;
	fp1 = fopen( "cut.txt", "a" ) ;
	fprintf( fp1, "%d\t",is ) ;
	fclose( fp1 ) ;
	*/
	// N E W 3====end

	pseg = chrom[ic].pseg ;
	lsg = chrom[ic].nseg ;
	len = (pseg + lsg -1)->end - pseg->beg ;

	cleft -= 1 - pow(pc,len) ;

	//
	nlinks2 -= links2( ic ) ;
	//

	// N E W ========
	//cleft2 -= 1 - pow(pc2,len) ;
	//call -= pow(pc2,len) ;
	// N E W =====end
	/* get seg # (jseg)  */

	for( jseg=0; is >= (pseg+jseg)->end ; jseg++) ;
	if( is >= (pseg+jseg)->beg ) in=1;
	else in=0;
	newsg = lsg - jseg ;

	/* copy last part of chrom to nchrom  */

	nchrom++;
	if( nchrom >= maxchr ) {
		maxchr += 20 ;
		chrom = (struct chromo *)realloc( chrom, (unsigned)(maxchr*sizeof(struct chromo))) ;
		if( chrom == NULL ) perror( "malloc error. segtre2");
	}
	if( !( pseg2 = chrom[nchrom-1].pseg = (struct seg *)calloc((unsigned)newsg,sizeof(struct seg)) ) )
	ERROR(" alloc error. re1");
	chrom[nchrom-1].nseg = newsg;
	chrom[nchrom-1].pop = chrom[ic].pop ;
	pseg2->end = (pseg+jseg)->end ;
	if( in ) {
		pseg2->beg = is + 1 ;
		(pseg+jseg)->end = is;
	}
	else pseg2->beg = (pseg+jseg)->beg ;
	pseg2->desc = (pseg+jseg)->desc ;
	for( k=1; k < newsg; k++ ) {
		(pseg2+k)->beg = (pseg+jseg+k)->beg;
		(pseg2+k)->end = (pseg+jseg+k)->end;
		(pseg2+k)->desc = (pseg+jseg+k)->desc;
	}

	lsg = chrom[ic].nseg = lsg-newsg + in ;
	lsgm1 = lsg - 1 ;
	nlinks -= pseg2->beg - (pseg+lsgm1)->end ;

	len = (pseg+lsgm1)->end - (pseg->beg) ;
	cleft += 1.0 - pow( pc, len) ;
	len = (pseg2 + newsg-1)->end - pseg2->beg ;
	cleft += 1.0 - pow(pc, len) ;

	// N E W ========
	//nlinks2 += pseg2->beg - (pseg+lsgm1)->end ;
	//nlinks2 -= pseg2->beg - (pseg+lsgm1)->end ;
	/*len = (pseg+lsgm1)->end - (pseg->beg) ;
	cleft2 += 1.0 - pow( pc2, len) ;
	len = (pseg2 + newsg-1)->end - pseg2->beg ;
	cleft2 += 1.0 - pow(pc2, len) ;

	len = (pseg+lsgm1)->end - (pseg->beg) ;
	call += pow( pc2, len) ;
	len = (pseg2 + newsg-1)->end - pseg2->beg ;
	call += pow(pc2, len) ;*/
	// N E W =====end

	//printf("ic = %d, lsg = %d\n", ic, lsg);
	if( !(chrom[ic].pseg =
		(struct seg *)realloc(chrom[ic].pseg,(unsigned)(lsg*sizeof(struct seg)) )) ) perror( " realloc error. re2");
		if( in ) {
			begs = pseg2->beg;
			for( i=0,k=0; (k<nsegs-1)&&(begs > seglst[seglst[i].next].beg-1);
			i=seglst[i].next, k++) ;
			if( begs != seglst[i].beg ) {
				/* new tree  */

				if( nsegs >= seglimit ) {
					seglimit += SEGINC ;
					nnodes = (int *)realloc( nnodes,(unsigned)(sizeof(int)*seglimit)) ;
					if( nnodes == NULL) perror("realloc error. 1. segtre_mig.c");
					seglst =
					(struct segl *)realloc( seglst,(unsigned)(sizeof(struct segl)*seglimit));
					if(seglst == NULL ) perror("realloc error. 2. segtre_mig.c");
					/*  printf("seglimit: %d\n",seglimit);  */
				}
				seglst[nsegs].next = seglst[i].next;
				seglst[i].next = nsegs;
				seglst[nsegs].beg = begs ;
				if( !(seglst[nsegs].ptree = (struct node *)calloc((unsigned)(2*nsam), sizeof(struct
					node)) )) perror("calloc error. re3.");
					nnodes[nsegs] = nnodes[i];
					ptree1 = seglst[i].ptree;
					ptree2 = seglst[nsegs].ptree;
					nsegs++ ;
					for( k=0; k<=nnodes[i]; k++) {
						(ptree2+k)->abv = (ptree1+k)->abv ;
						(ptree2+k)->time = (ptree1+k)->time;
					}
				}
			}

			// N E W ========
			nlinks2 += links2( ic ) ;
			nlinks2 += links2( nchrom-1 ) ;

			return(ic) ;
		}

		/***** common ancestor subroutine **********************
		Pick two chromosomes and merge them. Update trees if necessary. **/
		int ca ( nsam, nsites, c1, c2 )
		int nsam,c1,c2;
		int  nsites;
		{
			int yes1, yes2, seg1, seg2, seg ;
			int tseg, start, end, desc, k;
			struct seg *pseg ;
			struct node *ptree ;

			seg1 = 0 ;
			seg2 = 0 ;

			// generate a new list of segments for the ancestral chromosome
			if ( !(pseg = (struct seg *)calloc((unsigned)nsegs,sizeof(struct seg) )))
			perror("alloc error.ca1");

			tseg = -1 ;

			// for each segment
			for ( seg=0, k=0; k < nsegs; seg=seglst[seg].next, k++ ) {

				//===//
				// If no recombination,
				// start = 0
				// yes1 = yes2 = 1
				start = seglst[seg].beg ;
				yes1 = isseg (start, c1, &seg1) ;
				yes2 = isseg (start, c2, &seg2) ;

				//===//
				//printf ("    sta = %d\n", start) ;
				//printf ("    yes = %d, %d\n", yes1, yes2) ;

				// if either one has the segment
				if ( yes1 || yes2 ) {
					tseg++ ; // total number of segments
					(pseg+tseg)->beg=seglst[seg].beg ;
					end = ( k < nsegs-1 ? seglst[seglst[seg].next].beg-1 : nsites-1 ) ;
					(pseg+tseg)->end = end ;
					//printf ("  beg = %d, end = %d\n", (pseg+tseg)->beg, (pseg+tseg)->end) ;
					// beg = 0 , end = 1

					if ( yes1 && yes2 ) {
						nnodes[seg]++ ;
						if ( nnodes[seg] >= (2*nsam-2) ){ // last coalesent only
							tseg-- ; //printf ("A\n") ;
						}
						else {
							// tseg = 0 always. printf ("tseg = %d", tseg) ;
							(pseg+tseg)->desc = nnodes[seg] ; // printf ("B\n") ;
						}
						ptree=seglst[seg].ptree ;
						desc = (chrom[c1].pseg + seg1) -> desc ;
						(ptree+desc)->abv = nnodes[seg] ;
						//printf ("  abv = %d\n", (ptree+desc)->abv) ;
						desc = (chrom[c2].pseg + seg2) -> desc ;
						(ptree+desc)->abv = nnodes[seg] ;
						(ptree+nnodes[seg])->time = t ;
						//printf ("  abv = %d\n", (ptree+desc)->abv) ;

						if ( (c1 == 0) && (c2 == 1) ){
							//printf ("point = %d\n", nnodes[seg]) ;
						}
					}
					else { // not used
						(pseg+tseg)->desc = ( yes1 ?
							(chrom[c1].pseg + seg1)->desc :
							(chrom[c2].pseg + seg2)->desc);
						}
					}
				}
				//printf ("  tseg = %d\n", tseg) ;
				// tseg = 0 when sample > 2, = -1 when sample = -1 (last coalescent)



				nlinks -= links (c1) ;
				//printf ( "links(c1) = %d \n", links(c1) ) ;

				cleft -= 1.0 - pow(pc, (double)links(c1)) ;
				// N E W ========
				nlinks2 -= links2 (c1) ;
				/*cleft2 -= 1.0 - pow(pc2, (double)links(c1)) ;
				call -= pow(pc2, (double)links(c1)) ;*/
				// N E W ========

				free (chrom[c1].pseg) ;
				if ( tseg < 0 ) { // last coalescent
					free (pseg) ;
					chrom[c1].pseg = chrom[nchrom-1].pseg ;
					chrom[c1].nseg = chrom[nchrom-1].nseg ;
					chrom[c1].pop = chrom[nchrom-1].pop ;
					if ( c2 == nchrom-1 ) c2 = c1 ;
					nchrom--;
				}
				else {

					if( !(pseg = (struct seg *)realloc(pseg,(unsigned)((tseg+1)*sizeof(struct seg)))))
					perror(" realloc error. ca1") ;

					chrom[c1].pseg = pseg ;
					chrom[c1].nseg = tseg + 1 ;
					nlinks += links(c1) ;
					cleft += 1.0 - pow(pc, (double)links(c1));
					// N E W ========
					//printf ( "links(c1) = %d in if{} \n", links(c1) ) ;
					nlinks2 += links2(c1) ;
					/*cleft2 += 1.0 - pow(pc2, (double)links(c1)) ;
					call += pow(pc2, (double)links(c1)) ;*/
					// N E W ========
				}


				nlinks -= links(c2) ;
				//printf ( "links(c2) = %d \n", links(c2) ) ;

				cleft -= 1.0 - pow(pc, (double)links(c2)) ;
				// N E W ========
				nlinks2 -= links2(c2) ;
				/*cleft2 -= 1.0 - pow(pc2, (double)links(c2)) ;
				call -= pow(pc2, (double)links(c2)) ;*/
				// N E W ========
				free (chrom[c2].pseg) ;

				chrom[c2].pseg = chrom[nchrom-1].pseg;
				chrom[c2].nseg = chrom[nchrom-1].nseg;
				chrom[c2].pop = chrom[nchrom-1].pop ;
				nchrom--;

				if (tseg<0) return( 2 );  /* decrease of nchrom is two */ // last coalescent
				else return( 1 ) ;
			}

			/*** Isseg: Does chromosome c contain the segment on seglst which starts at
			start? *psg is the segment of chrom[c] at which one is to begin
			looking.  **/

			int isseg ( start, c, psg )
			int start, c, *psg;
			{
				int ns ;
				struct seg *pseg ;

				ns = chrom[c].nseg ;
				pseg = chrom[c].pseg ;

				/*  changed order of test conditions in following line on 6 Dec 2004 */
				for (  ; ((*psg) < ns ) && ( (pseg+(*psg))->beg <= start ) ; ++(*psg) )
				if ( (pseg+(*psg))->end >= start ) return(1);
				return (0);
			}



			int
			pick2_chrom(pop,config,pc1,pc2)
			int pop, *pc1, *pc2, config[];
			{
				int c1, c2, cs,cb,i, count;

				pick2(config[pop],&c1,&c2);
				cs = (c1>c2) ? c2 : c1;
				cb = (c1>c2) ? c1 : c2 ;
				i=count=0;
				for(;;){
					while( chrom[i].pop != pop ) i++;
					if( count == cs ) break;
					count++;
					i++;
				}
				*pc1 = i;
				i++;
				count++;
				for(;;){
					while( chrom[i].pop != pop ) i++;
					if( count == cb ) break;
					count++;
					i++;
				}
				*pc2 = i ;
			}



			/****  links(c): returns the number of links between beginning and end of chrom **/

			int
			links ( c )
			int c ;
			{
				int ns ;

				ns = chrom[c].nseg - 1 ;
				//printf ( "ns = %d, (chrom[c].pseg + ns)->end = %d, (chrom[c].pseg)->beg = %d\n", ns, (chrom[c].pseg + ns)->end, (chrom[c].pseg)->beg ) ;
				return( (chrom[c].pseg + ns)->end - (chrom[c].pseg)->beg);
			}

			// NEW
			/****  links2(c): returns the number of links between position 0 and end of effective regions in chrom for cin2 process**/

			int
			links2 ( c )
			int c ;
			{
				int ns2 ;

				ns2 = chrom[c].nseg - 1 ;
				//printf ( "ns = %d, (chrom[c].pseg + ns)->end = %d, (chrom[c].pseg)->beg = %d\n", ns, (chrom[c].pseg + ns)->end, (chrom[c].pseg)->beg ) ;
				return ( (chrom[c].pseg + ns2)->end ) ;
			}
