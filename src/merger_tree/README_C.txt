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

General 
=========

This code is producing merger trees according to the algorithm described in
Neistein & Dekel (2007). It was verified with the Millennium simulation data 
(Springel et al. 2005). Specifically, it was verified with the cosmological 
parameters (O_matter=0.25, O_Lambda=0.75, h=0.73, sigma_8=0.9), and for haloes 
in the mass range ~ 10^10 - 10^15 M_sun/h. 

Using the code for different halo masses, and for different cosmologies can be done, 
but you should note that this was not compared to results from N-body simulations yet. 
The code computes the variance of the smoothed density field, S(M), in any LCDM 
cosmology, according to the fitting function of van den Bosch (2002). 



C-Language aspects 
===================
The code is written in ANSI-C, so it should be easily compiled with a standard C
compiler (gcc) on any operating system.

The are 3 mathematical functions that are NOT given in this release because
of licensing issues. These are:


1. 
double random_number( prm_type *prm, int flag ) 
This function should return a random number.  
If flag==0, the random number should be uniformly distributed between 0 & 1.
If flag==1, the random number should be drawn from a normal distribution
with std=1 and average=0.
This function can use the parameter: 'prm.random_generator' which is of 
type 'long'.


2. 
double my_erf( double x ) ;
This function should give the resulting Error function for any real x.
For the standard definition of the Error function see for example:
http://mathworld.wolfram.com/Erf.html



3. 
double my_erfinv( double f ) ;
This function should give the resulting Inverse Error function for f,
where -1<=f<=1.


You can get these functions from few sources:
- Numerical Recipe (a license is needed), see http://www.nr.com.
- GSL library (free) : http://www.gnu.org/software/gsl/
- For the inverse error function you might want to look at:
   http://home.online.no/~pjacklam/notes/invnorm/
   W. J. Cody, Rational Chebyshev approximations for the error function, 
     Math. Comp., pp. 631-638, 1969



User parameters
================
The user may change the parameters of the code by modifying
the file "user_prm.txt". Please do not change the format of this 
file if you are not familiar with the code, or do not want to change it.
Each line of the file should contain a string ""xxxx" followed by some spaces
and a number. You cannot change the order of the parameters in this file.

The parameters description:

M0 		  - the mass of the tree 'trunk', in units of M_sun/h
ntrees 	          - number of merger-trees to run
nsteps            - number of time-step to run, each time-step is with \Delta\omega=0.1
mmin              - the minimum mass of progenitors in the tree, in units of M_sun/h
Omega_M           - the value of the cosmological matter density ("Omega matter")
Gamma_shape       - the value of the power-spectrum shape parameter (See app. A in Neistein & Dekel 2007)
sigma8            - the value of sigma_8
main_flag         - 0: full trees are constructed;  1: only main-progenitor histories are constructed
save_mode         - 0: save only the histogram of progenitors at the last time-step;  
                    1: save a table of all the progenitors at the last time-step (see "Output files" section); 
                    2: save a table as in (1), also save full trees  (see "Output files" section)
random_with_time  - 0: always use the same series of random numbers, 1: generate a new series in each run
hist_var          - variable to use in histograms and tables. 0: units of S(M), 1: units of M_sun/h
histogram_bins    - number of bins in the histogram above




Output files
=============
There are 3 types of output files the code can generate.


1. "output_histogram.txt"
This file is always updated.
It include the histogram of progenitors at the last 
time-step. Units are M_sun/h if hist_var=1, or S(M) if
hist_var=0. The number of bins in the histogram is fixed
by the parameter: histogram_bins


2. "all_progenitors.txt"
This file includes specific values for all progenitors at the last time-step.
It is created if save_mode=1 or save_mode=2.
Each line in the file includes the progenitors of only one tree.
Units are according to the parameter: hist_var


3. "output_fulltree1.txt"
This is a file that contains the full tree. 
There are "ntrees" such files saved in each run.
It is saved only if save_mode=2. It is significantly slowing the code.
An example file:

"
#indx, Mass, S, step-number, desc_indx 
0, 1e+014, 0.9704, 0, -1 
1, 9.367e+013,  0.999, 1, 0 
2, 2.797e+010,  13.17, 1, 0 
3, 1.802e+010,  14.55, 1, 0 
4, 1.13e+011,  9.406, 1, 0 
5, 3.916e+010,  12.18, 1, 0 
6, 6.6e+010,  10.75, 1, 0 
7, 2.309e+011,  7.816, 1, 0 
8, 6.882e+011,  5.776, 1, 0 
9, 2.004e+011,  8.114, 1, 0 
10, 2.574e+010,  13.42, 1, 0 
"




Performance
============

The performance stated below was verified with a computer working with a 
CPU of 1.4 GHz, and 512MB RAM.  The tests were done with nsteps=100, 
save_mode=0, and the cosmological parameters as quoted in the paper.

- main-progenitor histories:  ~1e-4 sec per 1 history
- Full tree :   ~ 4e-6 * (M0/mmin)  sec per tree 

As an example, with mmin=1.72e10, and M0=1e14, a full tree is produced at: 0.02 sec.
This tree has a total number of ~15,000 progenitors, and it describe a tree from
z=0 to z~8.






