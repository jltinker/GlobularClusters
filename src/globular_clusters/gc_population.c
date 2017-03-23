#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cutil.h"
#include "nrutil.h"
#include "hdr.h"

#define NUM_GC 100000


struct TREE1 {
  int nhalo;
  int nz;
  float *redshift;
  float *lookback_time;
  int *id;
  int *nd;
  int *desc;
  int *main_merger;
  float **mass;
  float **gasmass;
  float **mstar;
  float *gc;
  float *gc_age;
  float *metallicity;
  float *eta_circ;
} t;

struct GLOBULAR_CLUSTER {
  int n;
  float minit[NUM_GC];
  float mass[NUM_GC];
  float FeH[NUM_GC];
  float age[NUM_GC];
  float mgal[NUM_GC];
  float redshift[NUM_GC];
  int ihalo[NUM_GC];
  int mode[NUM_GC];
  float mhalo[NUM_GC];
} gc;

/* External functions
 */
float qromo(float (*func)(float), float a, float b,
	    float (*choose)(float(*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);
float zbrent(float (*func)(float), float x1, float x2, float tol);
void test_bimodality(float *x, float n);
float gasdev(long *idum);


/* function prototypes
 */
void read_tree(char *fname);
float primordial_gc(float hmass, float redshift);
void evolve_tree(void);
float fgas_lookup(float redshfhit, float mass);
float mstar_lookup(float redshift, float mass);
float func_calc_lookback_time(float z);
float random_circularity_parameter(void);
float lookback_time(float redshift);
int does_merger_happen(float m1, float m2, float time, float redshift, float eta_circ);
float metallicity(float mstar, float time);
float func_mmax(float lnmass);
void monte_carlo_gc_population(float mgc, float mstar, float time, float redshift, int mode, int ihalo, int iz);
float dynamically_evolve(float mass, float time);


int main(int argc, char **argv)
{
  // seed the randoms
  srand48(555);
  IDUM = 555;

  //get the input parameters
  input_params(argv);

  // read in the parsed tree
  read_tree(argv[2]);
  
  //evolve the tree forward in time
  evolve_tree();
}

/* at each time, keep track of all the halos. 
 * - put GCs in the halos when applicable (z=4?)
 * - determine when accretion events happen
 * - if accretion event happens, will dyn friction bring halo to center?
 */
void evolve_tree()
{
  FILE *fp;
  int i, ip, flag=1, j, nhalo, nmerge=0,itarg;
  float m1, m2, dm, mgc;


  nhalo = t.nhalo;

  // start at the second redshift entry
  for(i=t.nz;i>=2;--i)
    {
      // determine if there have been any accretion events
      // skip the first halo
      for(j=2;j<=nhalo;++j)
	{
	  // is a halo not here at this time but as at the previous?
	  if(t.mass[j][i-1]<=0 && t.mass[j][i]>0)
	    {
	      //what halo did it go into?
	      itarg = t.desc[j];
	      //is the target defined yet?
	      if(t.mass[itarg][i]<0)continue;
	      // can we merge in a Hubble time?
	      //printf("%d %d\n",itarg,j);
	      if(does_merger_happen(t.mass[j][i],t.mass[itarg][i],t.lookback_time[i], t.redshift[i], t.eta_circ[i]))
		{
		  // check to see if this one merges with the main halo
		  if(itarg==1) t.main_merger[j] = 1;
		}
	    }
	}

    }
  fprintf(stderr,"here\n");

  // need to walk the tree again and figure out which halos make it into the main
  for(i=2;i<=nhalo;++i)
    {
      ip = i;
      // walk forward til we get to desc==1
      while(t.desc[ip]!=1)
	ip = t.desc[ip];
      if(t.main_merger[ip]) t.main_merger[i] = 1; 
      //printf("%d %d\n",i,t.main_merger[i]);
    }
  fflush(stdout);

  // start the loop over, now can maintain one population of GCs,
  // just don't add GCs to the main halo population

  gc.n = 0;

  // start at the second redshift entry
  for(i=t.nz;i>=2;--i)
    {
      // determine if there have been any accretion events
      // skip the first halo
      for(j=1;j<=nhalo;++j)
	{
	  // if main halo, skip to disk growth
	  if(j==1)goto SKIP_MERGERS;
	  if(!t.main_merger[j])continue;

	  // is a halo not here at this time but as at the previous?
	  if(t.mass[j][i-1]<=0 && t.mass[j][i]>0)
	    {
	      //what halo did it go into?
	      itarg = t.desc[j];
	      //is the target defined yet?
	      if(t.mass[itarg][i]<0)continue;
	      // can we merge in a Hubble time?
	      if(does_merger_happen(t.mass[j][i],t.mass[itarg][i],t.lookback_time[i], t.redshift[i], t.eta_circ[i]))
		{
		  t.gc[itarg]+=t.gc[j];
		  t.gc_age[itarg]+=t.gc_age[j];
		  t.metallicity[itarg]+=t.metallicity[j];
		  nmerge++;
		  //printf("MASS %f %e %e %e %e %d\n",
		  //	 t.redshift[i],t.mass[itarg][i],t.gc[itarg],t.gc[1],t.gc[j],itarg);
		  // do we need to produce GCs in this merger?
		  if(MERGER_GROWTH)
		    {
		      if(t.gasmass[itarg][i]/(BARYON_FRACTION*t.mass[itarg][i])<0.04)continue;
		      m1 = t.mstar[itarg][i] + t.gasmass[itarg][i];
		      m2 = t.mstar[j][i] + t.gasmass[j][i];

		      //printf("%3d %10d %10d %f %f %e %e\n",i,itarg,j,t.redshift[i],m2/m1,m1,m2);
		      // is mass ratio high enough? 
		      if(m2/m1>MERGER_RATIO)
			{
			  dm = 3.0E6*(t.gasmass[itarg][i] + t.gasmass[j][i])/BARYON_FRACTION/1.0E11*MERGER_EFFICIENCY*
			    pow((1+t.redshift[i])/(1+4.3),1.0)* // no z-evolution
			    pow(t.mass[itarg][i]/MERGER_EFF_EVOLUTION_PIVOT,MERGER_EFF_EVOLUTION_SLOPE); // yes mass evolution
			  if(dm<1.0E+5)dm=0;
			  if(dm==0)continue;
			  t.gc[itarg] += dm;
			  t.gc_age[itarg] += dm*t.lookback_time[i];
			  t.metallicity[itarg]+=dm*metallicity(t.mstar[itarg][i],t.lookback_time[i]);
			  monte_carlo_gc_population(dm,t.mstar[itarg][i]+t.mstar[j][i],t.lookback_time[i],t.redshift[i],1,itarg,i);
			  if(itarg==-1)
			    fprintf(stdout,"MERGER %f %e %d %e %e %e %e %e %e %f\n",
				    t.redshift[i],t.mass[itarg][i],j,m1,m2,
				    t.gasmass[itarg][i],t.gasmass[j][i],dm,t.gc[itarg],
				    t.gc[itarg]/t.mstar[itarg][i]);
			}
		      continue;
		    }
		SKIP_MERGERS:
		  if(DISK_GROWTH)
		    {
		      if(t.gasmass[j][i]/(BARYON_FRACTION*t.mass[j][i])>DISK_THRESHOLD)
			{
			  dm = 3.06E6*(t.gasmass[j][i])/BARYON_FRACTION/1.0E11*DISK_EFFICIENCY*
			    fabs(t.lookback_time[i-1]-t.lookback_time[i])/0.1;
			  if(dm<1.0E+5)dm=0;
			  if(dm==0)continue;
			  t.gc[itarg] += dm;
			  t.gc_age[itarg] += dm*t.lookback_time[i];
			  t.metallicity[itarg]+=dm*metallicity(t.mstar[itarg][i],t.lookback_time[i]);
			  monte_carlo_gc_population(dm,t.mstar[itarg][i],t.lookback_time[i],t.redshift[i],2,itarg,i);
			}
		    }
		}
	    }
	}

    }

  // let output the population
  fp = fopen("gc_population.out","w");
  for(i=1;i<=gc.n;++i)   
    {
      gc.mass[i]=dynamically_evolve(gc.minit[i],gc.age[i]);
      if(gc.mass[i]>3.0E+3)
	fprintf(fp,"%e %e %f %f %e %e %f %d %7d\n",gc.mass[i], 
		gc.minit[i],gc.FeH[i],gc.age[i],gc.mgal[i],gc.mhalo[i],gc.redshift[i], gc.mode[i], gc.ihalo[i]);
    }
  mgc =0 ;
  for(i=1;i<=gc.n;++i)   
    if(gc.mass[i]>0)mgc+=gc.mass[i];
  fclose(fp);

  if(t.gc[1]>0)
    fprintf(stdout,"%e %e %e %e %f %f\n",t.mass[1][1],mgc,t.mstar[1][1],t.gasmass[1][1],
	    t.gc_age[1]/t.gc[1],t.metallicity[1]/t.gc[1]);
  else
    fprintf(stdout,"%e %e %e %e %f %f\n",t.mass[1][1],t.gc[1],t.mstar[1][1],t.gasmass[1][1],0.0,0.0);

}

float dynamically_evolve(float mass, float time)
{
  float vev_inv, d0=1.0/3.0;
  vev_inv = 10*pow(mass/2e5,((1+3*d0)/2));
  return mass*(0.6)*pow(1-(1+3*d0)/2/vev_inv*time,2./(1.+3*d0));
}

float func_mmax(float lnmass)
{
  return MGC - exp(lnmass)*(lnmass - log(MIN_GC_MASS));
}

/* randomly generate a GC population.
 * equations (11)-(13) in MG10 for generation.
 * 
 */
void monte_carlo_gc_population(float mgc, float mstar, float time, float redshift, int mode, int ihalo, int iz)
{
  float mmax, r, m, mtot =0;
  int n;

  n = gc.n;

  // set global
  MGC = mgc;

  // get most massive cluster
  mmax = zbrent(func_mmax, log(1.0E3), log(1.0E9), 1.0E-4);
  mmax = exp(mmax);
  n++;
  gc.minit[n] = mmax;
  gc.age[n] = time;
  gc.FeH[n] = metallicity(mstar, time);
  gc.mgal[n] = mstar;
  gc.redshift[n] = redshift;
  gc.mode[n] = mode;
  gc.ihalo[n] = ihalo;
  gc.mhalo[n] = t.mass[ihalo][iz];
  mtot += mmax;


  // start randomly assiging the additional clusters
  while (mtot < mgc)
    {
      r = drand48();
      m = MIN_GC_MASS/(1-r*(1-MIN_GC_MASS/mmax));
      if(mtot+m>mgc)
	if(drand48()>0.5)break;
      n++;
      gc.minit[n] = m;
      gc.age[n] = time;
      gc.FeH[n] = metallicity(mstar, time);
      gc.mgal[n] = mstar;
      gc.redshift[n] = redshift;
      gc.mode[n] = mode;
      gc.ihalo[n] = ihalo;
      gc.mhalo[n] = t.mass[ihalo][iz];
      mtot += m;
    }      
  fflush(stdout);
  gc.n = n;
}

/* equations (17_ and (18) in Muratov & Gnedin
 */
float metallicity(float mstar, float time)
{
  return -1.8 + 0.4*log10(mstar/1.0E+6) - 0.03*time;
}


// m1 is the larger one. the target. 
//(but let's make sure: sometimes these trees do weird things)

int does_merger_happen(float m2, float m1, float time, float redshift, float eta_circ)
{
  float t_dyn, t_merge, mt;
  
  // if we're not taking into account dynamical friction, then just assume
  // that the halos merge
  if(!SWITCH_DYNAMICAL_FRICTION) return 1;

  if(m2>m1) { mt = m1; m1 = m2; m2 = mt; }
  //eta_circ = random_circularity_parameter(); // now hard-wired per halo
  //get t_dyn
  t_dyn = 1.4*pow(1+redshift,-1.5);
  //get t_merge      
  t_merge = t_dyn*0.216*pow(m1/m2,1.3)/log(1+m1/m2)*
    exp(1.9*eta_circ)*(0.5/(1-eta_circ*eta_circ/2));
  // printf("MERGE %e %e %f %f %f\n",m1,m2,t_merge,time,redshift);
  //t_merge *= 2;

  //if circularity low, assume immediate merging.
  if(eta_circ<0.1)t_merge = 0;
  
  //see if sat will merge before z=0
  if(t_merge>time)
    return 0;
  return 1;
}

float merger_growth()
{
  return 1;
}

float primordial_gc(float hmass, float redshift)
{
  if(SWITCH_FGAS)
    return 0.006*BARYON_FRACTION*hmass*fgas_lookup(redshift,hmass);
  return 0.006*BARYON_FRACTION*hmass;
}

void read_tree(char *fname)
{
  FILE *fp;
  int i, idum, nz, nhalo, j;
  float xdum, x1;

  fp = openfile(fname);

  // get the number of halos and number of redshift bins
  fscanf(fp,"%d %d",&t.nhalo,&t.nz);
  nhalo = t.nhalo;
  nz = t.nz;
  t.redshift = vector(1,nz);
  t.lookback_time = vector(1,nz);

  // get the redshift information
  for(i=1;i<=t.nz;++i)
    fscanf(fp,"%d %f",&idum,&t.redshift[i]);
  for(i=1;i<=t.nz;++i)
    t.lookback_time[i] = lookback_time(t.redshift[i]);

  // allocate memory for the full tree
  fprintf(stderr,"Allocation [%.2f] MBytes for tree\n",2*t.nhalo*t.nz*sizeof(float)/1024./1024.);
  t.mass = matrix(1,nhalo,1,nz);
  t.mstar = matrix(1,nhalo,1,nz);
  t.gasmass = matrix(1,nhalo,1,nz);
  t.gc = vector(1,nhalo);
  t.gc_age = vector(1,nhalo);
  t.metallicity = vector(1,nhalo);
  t.eta_circ = vector(1,nhalo);

  t.id = ivector(1,nhalo);
  t.desc = ivector(1,nhalo);
  t.nd = ivector(1,nhalo);
  t.main_merger = ivector(1,nhalo);

  for(i=1;i<=t.nhalo;++i)
    {
      t.eta_circ[i] = random_circularity_parameter();
      t.main_merger[i] = t.gc[i] = t.gc_age[i] = t.metallicity[i] = 0;
      // get header information about halo
      fscanf(fp,"%d %d %d %f %f %f %f %d",&idum,&t.id[i],&t.nd[i],&xdum,&xdum,&xdum,&xdum,&t.desc[i]);
      //fprintf(stderr, "Halo %d %d\n",t.id[i],t.nd[i]);
      for(j=1;j<=t.nz;++j)t.mass[i][j]=-1;
      for(j=1;j<=t.nd[i];++j)
	{
	  fscanf(fp,"%f %f %d\n",&xdum,&x1,&idum);
	  t.mass[i][idum+1] = x1*0.7;
	  t.gasmass[i][idum+1] = fgas_lookup(t.redshift[idum+1],x1*0.7)*x1*0.7*BARYON_FRACTION;
	  t.mstar[i][idum+1] = mstar_lookup(t.redshift[idum+1],x1*0.7)*BARYON_FRACTION;
	  //printf("%f %e %e %e\n",t.redshift[idum+1],t.mass[i][idum+1],
	  //	 t.gasmass[i][idum+1]/BARYON_FRACTION/t.mass[i][idum+1],t.mstar[i][idum+1]);
	  if(SINGLE_GC_POPULATION==2)
	    if(j==t.nd[i])t.gc[i] = primordial_gc(t.mass[i][idum+1], t.redshift[idum+1]);
	}
    }
}

float fgas_lookup(float redshift, float mass)
{
  FILE *fp;
  float fg1,fg2,fg,x1,x2,x3,mdev;
  int i,j,iz;
  static int flag = 1;
  static float *zz, *mh, **fgas;
  static int nz, nm;

  if(flag)
    {
      if(MURATOV_LOOKUP)
	fp = openfile("fgas_lookup.dat");
      else
	fp = openfile("behroozi_lookup.dat");

      fscanf(fp,"%d %d",&nz,&nm);
      zz = vector(1,nz);
      mh = vector(1,nm);
      fgas = matrix(1,nz,1,nm);
      for(i=1;i<=nz;++i)
	fscanf(fp,"%f",&zz[i]);
      for(i=1;i<=nz;++i)
	for(j=1;j<=nm;++j)
	  {
	    fscanf(fp,"%f %f %f",&x1,&x2,&x3);
	    mh[j]=log(x1);
	    fgas[i][j] = log(x2);
	  }
      flag = 0;
      //fprintf(stderr,"Done with fgas readin\n");
    }

  mass = log(mass);
  // what is the upper and lower redshift?
  for(i=2;i<=nz;++i)
    if(zz[i]>redshift)break;
  iz = i;
  if(iz>nz) { iz = nz; }

  // NGP in mass;
  for(i=1;i<=nm;++i)
    if(mh[i]>mass || mh[nm]<mass) 
      { 
	fg1 = fgas[iz][i]; 
	fg2 = fgas[iz-1][i];
	break; 
      }
  fg = (fg1-fg2)/(zz[iz]-zz[iz-1])*(redshift-zz[iz-1]) + fg2;  if(MSTAR_SCATTER)

  // add in scatter in the gas fractions if needed
  if(MGAS_SCATTER)
    mdev = gasdev(&IDUM)*MGAS_SCATTER_VALUE*log(10);
  else
    mdev = 0;
  return(exp(fg+mdev));

}

float mstar_lookup(float redshift, float mass)
{
  FILE *fp;
  float fg1,fg2,fg,x1,x2,x3,mdev;
  int i,j,iz;
  static int flag = 1;
  static float *zz, *mh, **fgas;
  static int nz, nm;

  if(flag)
    {
      if(MURATOV_LOOKUP)
	fp = openfile("fgas_lookup.dat");
      else
	fp = openfile("behroozi_lookup.dat");
	
      fscanf(fp,"%d %d",&nz,&nm);
      zz = vector(1,nz);
      mh = vector(1,nm);
      fgas = matrix(1,nz,1,nm);
      for(i=1;i<=nz;++i)
	fscanf(fp,"%f",&zz[i]);
      for(i=1;i<=nz;++i)
	for(j=1;j<=nm;++j)
	  {
	    fscanf(fp,"%f %f %f",&x1,&x2,&x3);
	    mh[j]=log(x1);
	    fgas[i][j] = log(x3); // not fgas!-- now fstar!
	    //printf("%d %d %f %e %e\n",i,j,zz[i],exp(mh[j]),exp(fgas[i][j]));
	  }
      flag = 0;
      //fprintf(stderr,"Done with fstar readin\n");
    }

  mass = log(mass);
  // what is the upper and lower redshift?
  for(i=2;i<=nz;++i)
    if(zz[i]>redshift)break;
  iz = i;
  if(iz>nz) { iz = nz; }

  // NGP in mass;
  for(i=1;i<=nm;++i)
    if(mh[i]>mass || mass>mh[nm]) 
      { 
	fg1 = fgas[iz][i]; 
	fg2 = fgas[iz-1][i];
	break; 
      }
  fg = (fg1-fg2)/(zz[iz]-zz[iz-1])*(redshift-zz[iz-1]) + fg2;

  // add in scatter in the SHMR if required
  if(MSTAR_SCATTER)
    mdev = gasdev(&IDUM)*MSTAR_SCATTER_VALUE*log(10);
  else
    mdev = 0;

  //printf("LOOKUP %f %e %e %e %e %d %d\n",redshift, exp(mass),exp(fg1),exp(fg2),exp(fg),iz,i);
  return(exp(fg+mass+mdev));

}

float func_calc_lookback_time(float z)
{
  return 1/((1+z)*sqrt(OMEGA_M*(1+z)*(1+z)*(1+z)+(1-OMEGA_M)));
}
//returns t_lookback in Gyr (assumes H0=70)
float lookback_time(float redshift)
{
  return 13.9*qromo(func_calc_lookback_time,0,redshift,midpnt);
}

/* using the distribution function of
 * zentner+05 (eq. 5, a=2.22)
 */
float random_circularity_parameter(void)
{
  float p,pmax,eta;

  pmax = 8.692*pow(0.5*(1-0.5),1.22);
  p = 0;
  while(drand48()>p)
    {
      eta = drand48();
      p = 8.692*pow(0.5*(1-0.5),1.22)/pmax;
    }
  return eta;
}
