#include <stdlib.h>
#include "hdr.h"

//cosomology
float HUBBLE = 70;
float OMEGA_M = 0.3;
float BARYON_FRACTION = 0.15;

// options
int SWITCH_FGAS = 1;
int SWITCH_DYNAMICAL_FRICTION = 1;
int SINGLE_GC_POPULATION = 0;
int MERGER_GROWTH = 1;
int DISK_GROWTH = 0;
int MURATOV_LOOKUP = 0; // 0 means used behroozi_lookup table
int MSTAR_SCATTER = 1;
int MGAS_SCATTER = 1;

// free parameters
float MERGER_EFFICIENCY = 5; // boost factor for making GCs in mergers
float MERGER_RATIO = 0.2; // galaxy mass ratio required for making GCs in merger
float MERGER_EFF_EVOLUTION_PIVOT = 7; // redshift dependence of merger efficiency, pivot point
float MERGER_EFF_EVOLUTION_SLOPE = 0.0; // redshift dependence of merger efficiency, slope
float MERGER_EFF_MASS_PIVOT = 3.0E11; // mass dependence of merger efficiency, pivot point
float MERGER_EFF_MASS_SLOPE = -0.5; // mass dependence of merger efficiency, slope
float DISK_EFFICIENCY = 5; // boost factor for making GCs in gas-rich disks
float DISK_THRESHOLD = 0.9; // threshold in fgas for making GCs in disks
float MSTAR_SCATTER_VALUE = 0.2; // scatter in stellar mass at fixed halo mass (log10)
float MGAS_SCATTER_VALUE = 0.2; // scatter in gas mass at fixed halo mass (log10)

// random seeds
long IDUM = 555; // this one for gasdev()
long ISEED = 555; // this one for drand48()

// other globals
float MGC;
float MIN_GC_MASS = 1.0E+5;

