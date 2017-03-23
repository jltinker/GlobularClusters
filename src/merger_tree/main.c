/***************************************************************************
                                 Merger tree algorithm,  version 1.0
                             ------------------------------------------
    written                 : Monday, July 23, 2007
    last updated            : Wednesday, August 1, 2007
    by                      : Eyal Neistein, Hebrew University, Israel
    email                   : eyal_n@phys.huji.ac.il
    reference               : Neistein & Dekel, 2007, astro-ph 0708.1599
    web-page                : http://www.phys.huji.ac.il/~eyal_n/merger_tree/
    File list               : main.c, main.h, user_prm.txt
 ***************************************************************************/


#include "main.h"

float gasdev(long *idum);
double erfinv(double x);
double get_redshift(double x);
    
int main(int argc, char *argv[])
{
	// 'prm' structure include all parameters needed
	prm_type prm ;	
	// 'tree' holds the full tree
	tree_type *tree ;
	int i ;

	// parameters defined by the user
	user_parameters( &prm ) ;	
	// other model parameters, do not change!
	init_parameters( &prm ) ;
	
	// loop over all trees
	for ( i=1; i<=prm.ntrees ; i++ ) {			
		
		// if ( prm.main_flag==0 )
		//	printf( "Planting tree #%d out of %d ...", i, prm.ntrees ) ;

		// initialize tree array & allocate memory
		init_tree( &tree, &prm ) ;		
		// generate the tree
		construct_tree( tree, &prm, 0 ) ;
		// save the tree, or add to an histogram
		gather_save_output( tree, &prm, i ) ;	
		// free memory
		free( tree ) ;
		
		// if ( prm.main_flag==0 )
		//	printf( " finished\n" ) ;
	}

	// save histogram into a file
	save_hist( &prm ) ;
	printf( "\n\nDone. \n\n" ) ;
	

}



// ---------------------------------------------------------------------------------
void user_parameters( prm_type *prm )
{
	// Parameters to be changed by the user, see README file
	// Can either use input file: "user_prm.txt", or default prm encoded below

	FILE *f_ptr ;
	int i ;
	float s[20] ;
	char st[100] ;

	// default parameters
	prm->M0 =				2e13 ;
	prm->ntrees =			10 ;
	prm->nsteps =			30 ;
	prm->mmin =				1.72e10 ;
	prm->Omega_M =			0.25 ;
	prm->Gamma_shape =		0.169 ;
	prm->sigma8 =			0.9 ;	
	prm->main_flag =		0 ;	
	prm->save_mode =		0 ;		
	prm->random_with_time = 0 ;		
	prm->hist_var =			1 ;			
	prm->histogram_bins =	100 ;

	// read from file if exist
	if( (f_ptr  = fopen( "user_prm.txt", "r" )) == NULL )
		printf( "The file 'user_prm.txt' was not opened, using default parameters.\n\n" );
	else {
		printf( "The file 'user_prm.txt' was opened.\n\n" );
		for (i=1; i<=12 ; i++ ) {
			fscanf(f_ptr,"%s %f", &st, &s[i] ) ;
		}
		fclose(f_ptr) ;
		prm->M0 = s[1] ; 
		prm->ntrees = s[2] ; 
		prm->nsteps = s[3] ;
		prm->mmin =	s[4] ;
		prm->Omega_M = s[5] ;
		prm->Gamma_shape = s[6] ;
		prm->sigma8 = s[7] ;	
		prm->main_flag = s[8] ;	
		prm->save_mode = s[9] ;		
		prm->random_with_time = s[10] ;		
		prm->hist_var =	s[11] ;			
		prm->histogram_bins = s[12] ;
	}

	// print out which parameters are used
	printf("Parameters used:\n") ;
	printf( "M0 = %e\n",prm->M0 );
	printf( "ntrees = %d\n",prm->ntrees );
	printf( "nsteps = %d\n",prm->nsteps );
	printf( "mmin = %e\n",prm->mmin );
	printf( "Omega_M = %f\n",prm->Omega_M );
	printf( "Gamma_shape = %f\n",prm->Gamma_shape );
	printf( "sigma8 = %f\n",prm->sigma8 );
	printf( "main_flag = %d\n",prm->main_flag );
	printf( "save_mode = %d\n",prm->save_mode );
	printf( "random_with_time = %d\n",prm->random_with_time );
	printf( "hist_var = %d\n",prm->hist_var );
	printf( "histogram_bins = %d\n\n\n",prm->histogram_bins );

	
	 
	if ( prm->M0<1e6 || prm->M0>1e15 || prm->mmin<1e6 || prm->mmin>1e15 )
		printf( "Warning, all masses should be in the range [1e6,1e15] M_sun/h \n\n" ) ;
}


// ---------------------------------------------------------------------------------

void init_tree( tree_type **tree, prm_type *prm ) 

{
	// This function initialize the tree array

	int n ;

	// maximum tree size !
	n = ceil( prm->nsteps*prm->M0/prm->mmin ) ; 

	// allocate memory for the tree array
	*tree = malloc( n * sizeof(tree_type) ) ;
	if (*tree==NULL) printf("malloc failure, tree size is too big") ;

	// keep track of the last index used
	prm->last_tree_indx = 0 ;
	
	// initial values for the tree trunk
	(*tree)[0].descendant_id = -1 ;
	(*tree)[0].M = prm->M0 ;
	(*tree)[0].S = prm->S0 ;
	(*tree)[0].step = 0 ;

}



// ---------------------------------------------------------------------------------

void construct_tree( tree_type *tree, prm_type *prm, int ind )
{
	// Construct a tree starting from the index: 'ind' using a recursive scheme
	// This function generate all progenitors of the halo mass: tree[ind].M

	int n_progs, i, flag=0, desc_ind, flag2, count ;
	double x, Mleft, Sleft, f, logS0, logS0_2, xmax, a, b ;
	double sig_k, mu_k, sig_a, mu_a, M=0, S, M1, S0, Mp=0 ;
	double I, I_max ;

	long IDUM;

	// prepare variables
	S0 = tree[ind].S ;
	logS0 = log10( S0 ) ;
	logS0_2 = logS0*logS0 ;
	desc_ind = prm->last_tree_indx ; // save the index of the descendant

	// Kernel definitions, sig_k & mu_k
	sig_k = prm->model_prm[1] + prm->model_prm[2]*logS0 + prm->model_prm[3]*logS0_2 ;
	mu_k  = prm->model_prm[4] + prm->model_prm[5]*logS0 + prm->model_prm[6]*logS0_2 ;

	// get main progenitor 'S'
	//x = random_number( prm, 1 ) ;
	x = gasdev(&(prm->random_generator)); // random gaussian deviate
	S = S0 + exp( x*sig_k + mu_k ) ;

	// main progenitor is less than the minimum mass?
	if ( S>prm->smin )
		return ;

	// update values for the main progenitor mass
	M = ms_transform( S, prm, 2 ) ;
	M1 = M ;
	update_tree_array( tree, S, M, ind, prm ) ;
	n_progs = 1 ;

	// get small progs
	if ( prm->main_flag==0 ) {

		// Needed values for Mleft, mu_a & sig_a
		f = prm->model_prm[10] + prm->model_prm[11]*logS0 ;
		a = prm->model_prm[12] + prm->model_prm[13]*logS0 ;
		b = prm->model_prm[14] + prm->model_prm[15]*logS0 + prm->model_prm[16]*logS0_2 ;

		// add small progenitors until there is no mass left
		while ( flag==0 ) {

			// update Mleft
			Mp = Mp + M ;
			Mleft = f*tree[ind].M - Mp ;			
			if (Mleft>M1) Mleft=M1 ;

			if (Mleft<=prm->mmin*1.0) flag = 1 ;
			else {

				// update Sleft, sig_k & mu_k
				Sleft = ms_transform( Mleft, prm, 1 ) ;
				sig_a = sig_k + (Sleft-S0)*a ;
				mu_a  = mu_k  + (Sleft-S0)*b ;

				// the straight-forward approach, very time-consuming
				/*
				// this part ensures S will be below smin (or M>Mmin)
				xmax = ( log( prm->smin-Sleft ) - mu_a )/sig_a ;
				flag2 = 0 ; count=0 ;
				while (flag2==0) {
					x = random_number( prm, 1 ) ;
					count = count+1 ;
					if (count>1000000) x=xmax-0.0001 ; // only rare cases
					if (x<xmax) flag2=1 ;
				} 
				*/
				
				// Inverse probability approach, more fast 
				//I = random_number( prm, 0 ) ;
				I = drand48(); //uniform [0,1)
				I_max = 0.5 + 0.5*erf( (log(prm->smin-Sleft)-mu_a)/sqrt(2.0)/sig_a ) ;
				I = I*I_max ; // the integral is between 0 to maximum value
				x = sqrt(2.0)*erfinv(2.0*I-1.0) ; 
				

				// new progenitor
				S = Sleft + exp( x*sig_a + mu_a ) ;
				M = ms_transform( S, prm, 2 ) ;
				update_tree_array( tree, S, M, ind, prm ) ;
				n_progs = n_progs+1 ;
			
			}
		}
	}

	// check max number of steps
	if ( tree[ind].step==prm->nsteps-1 )
		return ;

	// repeat recursive for all progenitors
	for (i=1; i<=n_progs ; i++ ) 
		construct_tree( tree, prm, desc_ind+i ) ;
 
}



// ---------------------------------------------------------------------------------

void update_tree_array( tree_type *tree, double S, double M, int ind, prm_type *prm )
{
	// This funcation adds an additional progenitor to the tree
	// Writing location is at the end of the list

	int i_write ;

	i_write =  prm->last_tree_indx + 1 ;
	tree[i_write].descendant_id = ind ;
	tree[i_write].M = M ;
	tree[i_write].S = S ;
	tree[i_write].step = tree[ind].step+1 ;
	prm->last_tree_indx = prm->last_tree_indx + 1 ;

}


// ---------------------------------------------------------------------------------

void gather_save_output( tree_type *tree, prm_type *prm, int indx )
{
	// This function gather histogram from a tree
	// and save the tree into a file if needed

	FILE *f_ptr, *f_ptr1 ;
	int i, ind, ind0 ;
	char fname[100], fname1[100], istr[20] ;
	double x, d ;
	
	
	// Save all progenitors of last time-step?
	if ( prm->save_mode>0 ) {
		strcpy( fname1, "all_progenitors.txt" ) ;
		f_ptr1 = fopen( fname1, "a" ) ; 
		fprintf( f_ptr1, "%d ", indx ) ;
	}

	// check all tree (or last point if only main-progenitor)
	ind0 = 0 ;
	if ( prm->main_flag==1 ) ind0=prm->nsteps ;
	
	// go over tree indexes
	for (i=ind0; i<=prm->last_tree_indx; i++) {
		// does the progenitor have the appropriate time-step?
		if ( tree[i].step == prm->nsteps ) {
			// add to the histogram in S
			if ( prm->hist_var==0 ) {
				x = tree[i].S ;
				d = prm->hist_x[1] - prm->hist_x[0] ;
				ind = floor( (x-prm->S0)/d ) ;
			}
			if ( prm->hist_var==1 ) {	
				// histogram in M
				d = prm->hist_x[1] / prm->hist_x[0] ;
				x = tree[i].M/prm->M0 ;
				ind = floor( log(x/prm->hist_x[0]) / log(d) ) ;
			}

			// limits on bin index
			if ( ind<0 ) ind=0 ;
			if ( ind>=prm->histogram_bins ) ind=prm->histogram_bins-1 ;

			// add a count to the bin
			prm->hist_y[ind] = prm->hist_y[ind] + 1.0/prm->ntrees ;
			
			// write progenitor data into a file
			if ( prm->save_mode>0 ) fprintf( f_ptr1, ", %f ", x ) ;
		}
	}

	if ( prm->save_mode>0 ) {
		fprintf( f_ptr1," \n" ) ;
		fclose( f_ptr1 ) ;
	}

	// save full tree ?
	if ( prm->save_mode<2 )
		return ;	

	// write tree into a file
	strcpy( fname, "output_fulltree" ) ;
	sprintf( istr, "%d", indx ) ;
	strcat( fname , istr ) ;
	strcat( fname, ".txt" ) ;
	f_ptr = fopen( fname, "w" ) ;
	fprintf(f_ptr,"#indx Mass S step-number desc_indx \n" ) ;
	for (i=0; i<=prm->last_tree_indx; i++) {
	  fprintf(f_ptr, "%5d  %6.4e %10.4f %5d %5d %10.4f\n", i, tree[i].M, tree[i].S, tree[i].step, tree[i].descendant_id, get_redshift(tree[i].step*0.1+1.6865)) ;
	}
	fclose( f_ptr ) ;
	

}


// ---------------------------------------------------------------------------------

void save_hist( prm_type *prm )
{
	// Thie function writes the histogram into a file

	FILE *f_ptr ;
	int i ;
	char fname[100] ;

	// file name
	strcpy( fname, "output_histogram.txt" ) ;
	f_ptr = fopen( fname, "w" ) ;

	// heading
	if ( prm->hist_var==0 )
		fprintf(f_ptr,"#indx, S, #-of-progs/N_trees \n" ) ;
	if ( prm->hist_var==1)
		fprintf(f_ptr,"#indx, M/M0, #-of-progs/N_trees \n" ) ;

	// write data
	for (i=0; i<prm->histogram_bins; i++) {
		fprintf(f_ptr, "%d, %6.4g, %6.4g \n", i, prm->hist_x[i], prm->hist_y[i] ) ;
	}
	fclose( f_ptr ) ;
	
	free( prm->hist_x ) ;
	free( prm->hist_y ) ;


}


// ---------------------------------------------------------------------------------

void init_parameters( prm_type *prm )

{
	// Initialize all parameters

	int i, j, inds, lutsize ;
	double M, S, logS, dlogS, logM, xmin, d ;
	FILE *f_ptr ;
	char fname[100] ;
	
	// initialize progenitors file if needed
	if ( prm->save_mode>0 ) {
		strcpy( fname, "all_progenitors.txt" ) ;
		f_ptr = fopen( fname, "w" ) ;
		if ( prm->hist_var==0 )
			fprintf(f_ptr,"Tree indx, S1, S2, ... \n" ) ;
		if ( prm->hist_var==1)
			fprintf(f_ptr,"Tree indx, M1/M0, M2/M0, ...  \n" ) ;
		fclose( f_ptr ) ;
	}

	// model parameters, K_1
	prm->model_prm[1] =		1.367 ;
	prm->model_prm[2] =		0.012 ;
	prm->model_prm[3] =		0.234 ;
	prm->model_prm[4] =		-3.682 ;
	prm->model_prm[5] =		0.760 ;
	prm->model_prm[6] =		-0.360 ;

	// model parameters, fraction of mass lost, f
	prm->model_prm[10] =		0.967 ;
	prm->model_prm[11] =		-0.0245 ;

	// model parameters, K_a extension
	prm->model_prm[12] =		0.104 ;
	prm->model_prm[13] =		0.118 ;
	prm->model_prm[14] =		2.700 ;
	prm->model_prm[15] =		-4.760 ;
	prm->model_prm[16] =		2.900 ;

	// random seed
    if ( prm->random_with_time==1 ) {
		#ifdef HAVE_CONFIG_H 
			// Linux
			srand((unsigned int)time((time_t *)NULL)) ;
		#else
			// windows
			srand( (unsigned)time( NULL ) );
		#endif
		prm->random_generator = rand() ;
	} else 
      prm->random_generator = -(prm->random_with_time) ; // JLT


	//////////////////////////////////////////////////////////////////////
	//   S <--> M, two tables are created for fast pick-up
	//////////////////////////////////////////////////////////////////////
	lutsize = 1000 ;
	// create look-up-table for M -->>> S transormation
	for ( i=0; i<lutsize ; i++ ) {  
		M = pow(10, 6.01+0.01*i ) ;
		S = m2s( M, prm->Omega_M, prm->sigma8, prm->Gamma_shape ) ;
		prm->lut1_logM[i] = log10(M) ;
		prm->lut1_logS[i] = log10(S) ;
	}

	
	// create look-up-table for S -->>> M transormation
	// this is going from min(S) up to max(S) in equal step of logS
	dlogS = ( prm->lut1_logS[0]-prm->lut1_logS[lutsize-1] )/ (lutsize-1) ;
	prm->lut2_logM[0] = prm->lut1_logM[lutsize-1] ;
	prm->lut2_logS[0] = prm->lut1_logS[lutsize-1] ;
	prm->lut2_logM[lutsize-1] = prm->lut1_logM[0] ;
	prm->lut2_logS[lutsize-1] = prm->lut1_logS[0] ;
	for ( i=1; i<lutsize-1 ; i++ ) {  
		logS = prm->lut1_logS[lutsize-1] + i*dlogS ; 
		// find M according to previous lut
		inds = 0 ;
		for ( j=0; j<lutsize-1 ; j++ ) {
			if ( prm->lut1_logS[j]>logS && prm->lut1_logS[j+1]<=logS ) inds = j ;
		}

		logM = ( prm->lut1_logM[inds+1] - prm->lut1_logM[inds] ) /
			   ( prm->lut1_logS[inds+1] - prm->lut1_logS[inds] ) *
			   ( logS - prm->lut1_logS[inds] ) + prm->lut1_logM[inds] ;
		
		prm->lut2_logM[i] = logM ;
		prm->lut2_logS[i] = logS ;
	}

	// S stuff 
	prm->S0 = ms_transform( prm->M0, prm, 1 ) ;
	prm->smin = ms_transform( prm->mmin, prm, 1 ) ;

	/////////////////////////////////////////////////////////////////////////////

	// histogram bins
	prm->hist_x = malloc( prm->histogram_bins * sizeof(double) ) ;
	prm->hist_y = malloc( prm->histogram_bins * sizeof(double) ) ;
	memset( prm->hist_x, 0, prm->histogram_bins*sizeof(double)   ) ;
	memset( prm->hist_y, 0, prm->histogram_bins*sizeof(double)   ) ;

	// linear bins in S
	if ( prm->hist_var==0 ) {
		for ( i=0; i<prm->histogram_bins ; i++ ) 
			prm->hist_x[i] = prm->S0 + (i+0.5)*(prm->smin-prm->S0)/prm->histogram_bins ;
	}

	// log bins in M/M0
	if ( prm->hist_var==1 ) {
		xmin = prm->mmin/prm->M0 ;
		d = pow( xmin, -1.0/prm->histogram_bins ) ;
		for ( i=0; i<prm->histogram_bins ; i++ ) 
			prm->hist_x[i] = xmin * pow(d,i+0.5) ;
	}

}


// ---------------------------------------------------------------------------------

double ms_transform( double MS, prm_type *prm, int tmode )
{
	// transfrom M->S if tmode = 1
	// transfrom S->M if tmode = 2
	// Simple linear interpolation is used

	double logMS, d, x1, x2, y1, y2, mmin, smin, y, logy ;
	int ind ;

	logMS = log10( MS ) ;

	if ( tmode==1 ) {
		mmin = prm->lut1_logM[0] ;
		d = prm->lut1_logM[1]-prm->lut1_logM[0] ;
		ind = floor( (logMS-mmin)/d ) ;
		x1 = prm->lut1_logM[ind] ;
		x2 = prm->lut1_logM[ind+1] ;
		y1 = prm->lut1_logS[ind] ;
		y2 = prm->lut1_logS[ind+1] ;

	} else if ( tmode==2 ) {

		smin = prm->lut2_logS[0] ;
		d = prm->lut2_logS[1]-prm->lut2_logS[0] ;
		ind = floor( (logMS-smin)/d ) ;
		x1 = prm->lut2_logS[ind] ;
		x2 = prm->lut2_logS[ind+1] ;
		y1 = prm->lut2_logM[ind] ;
		y2 = prm->lut2_logM[ind+1] ;

	}

	logy = (y2-y1)/(x2-x1)*(logMS-x1) + y1 ;

	y = pow( 10, logy ) ;

	return y ;

}

// ---------------------------------------------------------------------------------

double m2s( double M, double Omega_M, double sigma8, double Gamma_shape )

{
	// See vdBosch (2002), Appendix
	// it is also described in Neistein & Dekel (2007)

	double S=0, c1, u8, fu8, u, fu, sigma ;
	double M1 ;
	
	c1 = 3.804e-4 ;
	u8 = 32 * Gamma_shape ;
	fu8 = fitting_fun( u8 ) ;

	M1 = M/Omega_M ;
	u =  c1 * Gamma_shape * pow( M1, 0.33333 ) ;
	fu = fitting_fun( u ) ;
	sigma = sigma8 * fu / fu8 ;
	S = sigma * sigma ;

	return S ;
}



// ---------------------------------------------------------------------------------

double fitting_fun( double x )

{
	double g, a2, a3, a4, a5, a6 ;

	a2 = 8.8e-4 * pow( x, 0.2 ) ; // small correction for vdBosch (2002) formula
	a3 = 1.074  * pow( x, 0.3 ) ;
	a4 = -1.581 * pow( x, 0.4 ) ;
	a5 = 0.954  * pow( x, 0.5 ) ;
	a6 = -0.185 * pow( x, 0.6 ) ;
	g = 64.087 * pow( 1 + a2 + a3 + a4 + a5 + a6 , -10.0 )  ;

	return g ;
}

