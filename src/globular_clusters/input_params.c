#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hdr.h"

/* read in the parameters from the input file.
 * THe order of the parameters will be the same order as
 * listed in the hdr.c file
 */

/* 
//cosomology
extern float HUBBLE = 70;
extern float OMEGA_M = 0.3;
extern float BARYON_FRACTION = 0.15;

// options
extern int SWITCH_FGAS = 1;
extern int SWITCH_DYNAMICAL_FRICTION = 1;
extern int SINGLE_GC_POPULATION = 0;
extern int MERGER_GROWTH = 1;
extern int DISK_GROWTH = 0;
extern int MURATOV_LOOKUP = 0;
extern int MSTAR_SCATTER = 1;
extern int MGAS_SCATTER = 1;

float MERGER_EFFICIENCY = 5; // boost factor for making GCs in mergers
float MERGER_RATIO = 0.2; // galaxy mass ratio required for making GCs in merger
float MERGER_EFF_EVOLUTION_PIVOT = 3.0E11; // mass dependence of merger efficiency, pivot point
float MERGER_EFF_EVOLUTION_SLOPE = -0.5; // mass dependence of merger efficiency, slope
float DISK_EFFICIENCY = 5; // boost factor for making GCs in gas-rich disks
float DISK_THRESHOLD = 0.9; // threshold in fgas for making GCs in disks
float MSTAR_SCATTER_VALUE = 0.2; // scatter in stellar mass at fixed halo mass (log10)
float MGAS_SCATTER_VALUE = 0.2; // scatter in gas mass at fixed halo mass (log10)

// other globals 
extern long IDUM = 555;
extern long ISEED = 555;
*/


#define gl() fgets(aa,1000,fp)

int input_params(char **argv)
{
  FILE *fp;
  char aa[1000];

  fp = openfile(argv[1]);

  // get the cosmology parameters
  fscanf(fp,"%f",&HUBBLE);gl(); 
  fscanf(fp,"%f",&OMEGA_M);gl(); 
  fscanf(fp,"%f",&BARYON_FRACTION);gl(); 

  // get the various GC options
  fscanf(fp,"%d",&SWITCH_FGAS);gl(); 
  fscanf(fp,"%d",&SWITCH_DYNAMICAL_FRICTION);gl(); 
  fscanf(fp,"%d",&SINGLE_GC_POPULATION);gl(); 
  fscanf(fp,"%d",&MERGER_GROWTH);gl(); 
  fscanf(fp,"%d",&DISK_GROWTH);gl(); 
  fscanf(fp,"%d",&MURATOV_LOOKUP);gl(); 
  fscanf(fp,"%d",&MSTAR_SCATTER);gl(); 
  fscanf(fp,"%d",&MGAS_SCATTER);gl(); 

  //free parameters
  fscanf(fp,"%f",&MERGER_EFFICIENCY);gl(); 
  fscanf(fp,"%f",&MERGER_RATIO);gl(); 
  fscanf(fp,"%f",&MERGER_EFF_EVOLUTION_PIVOT);gl(); 
  fscanf(fp,"%f",&MERGER_EFF_EVOLUTION_SLOPE);gl(); 
  fscanf(fp,"%f",&DISK_EFFICIENCY);gl(); 
  fscanf(fp,"%f",&DISK_THRESHOLD);gl(); 
  fscanf(fp,"%f",&MSTAR_SCATTER_VALUE);gl(); 
  fscanf(fp,"%f",&MGAS_SCATTER_VALUE);gl(); 

  // random seeds
  fscanf(fp,"%ld",&IDUM);gl(); 
  fscanf(fp,"%ld",&ISEED);gl(); 
  

}
