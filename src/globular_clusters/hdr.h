#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Some function prototypes
 */
int filesize(FILE *fp);
FILE *openfile(char *ff);
int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);


/* Global variables
 */

//cosomology
extern float HUBBLE;
extern float OMEGA_M;
extern float BARYON_FRACTION;

// options
extern int SWITCH_FGAS;
extern int SWITCH_DYNAMICAL_FRICTION;
extern int SINGLE_GC_POPULATION;
extern int MERGER_GROWTH;
extern int DISK_GROWTH;
extern int MURATOV_LOOKUP;
extern int MSTAR_SCATTER;
extern int MGAS_SCATTER;

// free parameters
extern float MERGER_EFFICIENCY; // boost factor for making GCs in mergers
extern float MERGER_RATIO ; // galaxy mass ratio required for making GCs in merger
extern float MERGER_EFF_EVOLUTION_PIVOT; // redshift dependence of merger efficiency, pivot point
extern float MERGER_EFF_EVOLUTION_SLOPE; // redshift dependence of merger efficiency, slope
extern float MERGER_EFF_MASS_PIVOT; // mass dependence of merger efficiency, pivot point
extern float MERGER_EFF_MASS_SLOPE; // mass dependence of merger efficiency, slope
extern float DISK_EFFICIENCY; // boost factor for making GCs in gas-rich disks
extern float DISK_THRESHOLD; // threshold in fgas for making GCs in disks
extern float MERGER_RATIO_POWER;
extern float MSTAR_SCATTER_VALUE; // scatter in stellar mass at fixed halo mass (log10)
extern float MGAS_SCATTER_VALUE; // scatter in gas mass at fixed halo mass (log10)
extern float EFFICIENCY_SCATTER;

// other globals 
extern long IDUM;
extern long ISEED;

// other globals (not to be read in from batch file)
extern float MIN_GC_MASS;
extern float MGC;
