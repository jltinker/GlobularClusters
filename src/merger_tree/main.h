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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <malloc.h>
#include <memory.h>
#include <string.h>

typedef
struct _prm_type {
	double M0 ;
	double S0 ;
	int ntrees ;					
	int nsteps ;					
	double mmin ;
	double smin ;
	double Omega_M ;					
	double Gamma_shape ;			
	double sigma8 ;				
	int main_flag ;					
	int save_mode ;
	int random_with_time ;
	int histogram_bins ;
	int hist_var ;
	int last_tree_indx ;
	long random_generator ;
	double model_prm[20] ;
	double *hist_x ;
	double *hist_y ;
	double lut1_logM[1000] ;
	double lut1_logS[1000] ;
	double lut2_logM[1000] ;
	double lut2_logS[1000] ;
} prm_type ;

typedef
struct _tree_type {
	double M ;			
	double S ;					
	int descendant_id ;	
	int step ;
} tree_type ;


void user_parameters( prm_type *prm ) ;
void init_parameters( prm_type *prm ) ;
double m2s( double M, double Omega_M, double sigma8, double Gamma_shape ) ;
double fitting_fun( double x ) ;
void init_tree( tree_type **tree, prm_type *prm ) ;
void construct_tree( tree_type *tree, prm_type *prm, int ind ) ;
double ms_transform( double MS, prm_type *prm, int tmode ) ;
void update_tree_array( tree_type *tree, double S, double M, int ind, prm_type *prm ) ;
void gather_save_output( tree_type *tree, prm_type *prm, int indx ) ;
void save_hist( prm_type *prm ) ;
double random_number( prm_type *prm, int flag ) ; // random number generator
double my_erf( double x ) ; // Erf function, from GSL (http://www.gnu.org/software/gsl/) or Numerical recipes (http://www.nr.com)
double my_erfinv( double f ) ; // See e.g. http://home.online.no/~pjacklam/notes/invnorm/
                               // also: W. J. Cody, Rational Chebyshev approximations for the error function, Math. Comp., pp. 631-638, 1969
