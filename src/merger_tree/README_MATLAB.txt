/***************************************************************************
                                 Merger tree algorithm,  version 1.0
                             ------------------------------------------
    written                 : Monday, July 23, 2007
    last updated            : Wednesday, August 1, 2007
    by                      : Eyal Neistein, Hebrew University, Israel
    email                   : eyal_n@phys.huji.ac.il
    reference               : Neistein & Dekel, 2007, astro-ph 0708.1599
    web-page                : http://www.phys.huji.ac.il/~eyal_n/merger_tree/
    File list               : merger_tree_algo.m
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
 



MATLAB-Language aspects 
=======================
The code was written with MATLAB version 7.0 (R14). 



User parameters
================
The user may change the parameters of the code by modifying
the input structure: "user_prm". 

The parameters description:

M0 		  - the mass of the tree 'trunk', in units of M_sun/h
ntrees 	          - number of merger-trees to run
nsteps            - number of time-step to run, each time-step is with \Delta\omega=0.1
mmin              - the minimum mass of progenitors in the tree, in units of M_sun/h
Omega_M           - the value of the cosmological matter density ("Omega matter")
Gamma_shape       - the value of the power-spectrum shape parameter (See app. A in Neistein & Dekel 2007)
sigma8            - the value of sigma_8
main_flag         - 0: full trees are constructed;  1: only main-progenitor histories are constructed
random_with_time  - 0: always use the same series of random numbers, 1: generate a new series in each run


Remarks:

- 'M0' may be a vector with length=ntrees. In this case a set of trees will be produced, according to 'M0' values.
-  If 'nsteps' is a vector of length 1, output will be save at all time-steps. If 'nsteps' is a vector with more
than one number, the trees will be produced up to max(nsteps), and will be saved only at values of 'nsteps'.
To save output only at the last time-step use for example:  nsteps=[10 10] ;



Output Structure
=================

The code output is as follows:

out.s_array{i}[j,:]       = the values of S(M), for all progenitors at the time-step 'i', and for the tree 'j'
out.m_array{i}[j,:]       = the values of M (M_sun/h), for all progenitors at the time-step 'i', and for the tree 'j'
out.desc_ind{i}[j,:]      = the index of the descendant at the previous time-step for the progenitor described above
out.prm                   = parameter sturcture used in the code



Performance
============


The performance stated below was verified with a computer working with a 
CPU of 1.4 GHz, and 512MB RAM.  The tests were done with nsteps=100, 
save_mode=0, and the cosmological parameters as quoted in the paper.

- main-progenitor histories:  ~5e-4 sec per 1 history
- Full tree :   ~ 4e-4 * (M0/mmin)  sec per tree 

As an example, with mmin=1.72e10, and M0=1e14, a full tree is produced at: ~2 sec.
This tree has a total number of ~15,000 progenitors, and it describe a tree from
z=0 to z~8.




