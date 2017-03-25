#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hdr.h"

/* read in the parameters from the input file.
 * THe order of the parameters will be the same order as
 * listed in the hdr.c file
 */

#define gl() fgets(aa,1000,fp)

int input_params(char **argv)
{
  int n;
  FILE *fp;
  char aa[1000];
  float x1;

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
  fscanf(fp,"%f",&DISK_EFFICIENCY);gl(); 
  fscanf(fp,"%f",&DISK_THRESHOLD);gl(); 
  fscanf(fp,"%f",&MERGER_EFF_EVOLUTION_PIVOT);gl(); 
  fscanf(fp,"%f",&MERGER_EFF_EVOLUTION_SLOPE);gl(); 
  fscanf(fp,"%f",&MERGER_EFF_MASS_PIVOT);gl(); 
  fscanf(fp,"%f",&MERGER_EFF_MASS_SLOPE);gl(); 
  fscanf(fp,"%f",&MSTAR_SCATTER_VALUE);gl(); 
  fscanf(fp,"%f",&MGAS_SCATTER_VALUE);gl(); 

  // random seeds
  fscanf(fp,"%ld",&IDUM);gl(); 
  fscanf(fp,"%ld",&ISEED);gl(); 
  
  // check input:
  if(MERGER_EFF_MASS_PIVOT < 100)
    MERGER_EFF_MASS_PIVOT = pow(10.0,MERGER_EFF_MASS_PIVOT);

}
