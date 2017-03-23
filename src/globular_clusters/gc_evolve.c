#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cutil.h"
#include "nrutil.h"

//cosomology
float HUBBLE = 70;
float OMEGA_M = 0.3;
float BARYON_FRACTION = 0.15;
// options
int SWITCH_FGAS = 1;
int SWITCH_DYNAMICAL_FRICTION = 1;
int SINGLE_GC_POPULATION = 0;
int MERGER_GROWTH = 1;

struct TREE1 {
  int nhalo;
  int nz;
  float *redshift;
  float *lookback_time;
  int *id;
  int *nd;
  int *desc;
  float **mass;
  float **gasmass;
  float **mstar;
  float *gc;
  float *gc_age;
  float *metallicity;
} t;

/* External functions
 */
float qromo(float (*func)(float), float a, float b,
	    float (*choose)(float(*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);


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
int does_merger_happen(float m1, float m2, float time, float redshift);
float metallicity(float mstar, float time);

int main(int argc, char **argv)
{
  // read in the parsed tree
  read_tree(argv[1]);
  
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
  int i, flag=1, j, nhalo, nmerge=0,itarg;
  float m1, m2, dm;


  nhalo = t.nhalo;

  // start at the second redshift entry
  for(i=t.nz;i>=2;--i)
    {
      // first case: all progenitor GCs at z=4. no further growth
      if(SINGLE_GC_POPULATION==1)
	{
	  if(t.redshift[i]>6)continue;
	  if(flag) 
	    {
	      flag = 0;
	      //fprintf(stderr,"Starting with primordial GCs\n");
	      for(j=1;j<=nhalo;++j)
		if(t.mass[j][i]>0) t.gc[j] = primordial_gc(t.mass[j][i], t.redshift[i]);
	      //for(j=1;j<=nhalo;++j)
	      //fprintf(stderr,"%d %f %e %e\n",i,t.redshift[i],t.mass[j][i],t.gc[j]);
	      //fprintf(stderr,"Done with primordial GCs\n");
	      continue;
	    }
	}
      // determine if there have been any accretion events
      // skip the first halo
      for(j=2;j<=nhalo;++j)
	{
	  // is a halo not here at this time but as at the previous?
	  if(t.mass[j][i-1]<=0 && t.mass[j][i]>0)
	    {
	      //if(j==2)
	      //	printf("HALO %d %d %e %e %d %d\n",j,i,t.mass[j][i-1],t.mass[j][i],i,
	      //       does_merger_happen(t.mass[j][i],t.mass[itarg][i],t.lookback_time[i], t.redshift[i]));

	      //what halo did it go into?
	      itarg = t.desc[j];
	      //is the target defined yet?
	      if(t.mass[itarg][i]<0)continue;
	      // can we merge in a Hubble time?
	      if(does_merger_happen(t.mass[j][i],t.mass[itarg][i],t.lookback_time[i], t.redshift[i]))
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
		      m1 = t.mstar[itarg][i] + t.gasmass[itarg][i];
		      m2 = t.mstar[j][i] + t.gasmass[j][i];

		      //printf("%3d %10d %10d %f %f %e %e\n",i,itarg,j,t.redshift[i],m2/m1,m1,m2);
		      // is mass ratio high enough? 
		      if(m2/m1>0.2)
			{
			  dm = 3.0E6*(t.gasmass[itarg][i] + t.gasmass[j][i])/BARYON_FRACTION/1.0E11*5;
			  if(dm<1.0E+5)dm=0;
			  t.gc[itarg] += dm;
			  t.gc_age[itarg] += dm*t.lookback_time[i];
			  t.metallicity[itarg]+=dm*metallicity(t.mstar[itarg][i],t.lookback_time[i]);
			  if(itarg==-1)
			    fprintf(stdout,"MERGER %f %e %d %e %e %e %e %e %e %f\n",
				    t.redshift[i],t.mass[itarg][i],j,m1,m2,
				    t.gasmass[itarg][i],t.gasmass[j][i],dm,t.gc[itarg],
				    t.gc[itarg]/t.mstar[itarg][i]);
			}
		    }
		}
	    }
	}

    }
  if(t.gc[1]>0)
    fprintf(stdout,"%e %e %e %e %f %f\n",t.mass[1][1],t.gc[1],t.mstar[1][1],t.gasmass[1][1],
	    t.gc_age[1]/t.gc[1],t.metallicity[1]/t.gc[1]);
  else
    fprintf(stdout,"%e %e %e %e %f %f\n",t.mass[1][1],t.gc[1],t.mstar[1][1],t.gasmass[1][1],0.0,0.0);

}

/* equations (17_ and (18) in Muratov & Gnedin
 */
float metallicity(float mstar, float time)
{
  return -1.8 + 0.4*log10(mstar/1.0E+6) - 0.03*time;
}


// m1 is the larger one. the target. 
//(but let's make sure: sometimes these trees do weird things)

int does_merger_happen(float m2, float m1, float time, float redshift)
{
  float eta_circ, t_dyn, t_merge, mt;
  
  // if we're not taking into account dynamical friction, then just assume
  // that the halos merge
  if(!SWITCH_DYNAMICAL_FRICTION) return 1;

  if(m2>m1) { mt = m1; m1 = m2; m2 = mt; }
  eta_circ = random_circularity_parameter();
  //get t_dyn
  t_dyn = 1.4*pow(1+redshift,-1.5);
  //get t_merge      
  t_merge = t_dyn*0.216*pow(m1/m2,1.3)/log(1+m1/m2)*
    exp(1.9*eta_circ)*(0.5/(1-eta_circ*eta_circ/2));
  // TESTING
  //t_merge = 20;
  //printf("MERGE %e %e %f %f %f\n",m1,m2,t_merge,time,redshift);

  //if circularity low, assume immediate merging.
  if(eta_circ<0.2)t_merge = 0;
  
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

  t.id = ivector(1,nhalo);
  t.desc = ivector(1,nhalo);
  t.nd = ivector(1,nhalo);

  for(i=1;i<=t.nhalo;++i)
    {
      t.gc[i] = t.gc_age[i] = t.metallicity[i] = 0;
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
	  //printf("%f %e %e %e\n",t.redshift[idum+1],t.mass[i][idum+1],t.gasmass[i][idum+1],t.mstar[i][idum+1]);
	  if(SINGLE_GC_POPULATION==2)
	    if(j==t.nd[i])t.gc[i] = primordial_gc(t.mass[i][idum+1], t.redshift[idum+1]);
	}
    }
}

float fgas_lookup(float redshift, float mass)
{
  FILE *fp;
  float fg1,fg2,fg,x1,x2,x3;
  int i,j,iz;
  static int flag = 1;
  static float *zz, *mh, **fgas;
  static int nz, nm;

  if(flag)
    {
      fp = openfile("fgas_lookup.dat");
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

  // NGP in mass;
  for(i=1;i<=nm;++i)
    if(mh[i]>mass || mh[nm]<mass) 
      { 
	fg1 = fgas[iz][i]; 
	fg2 = fgas[iz-1][i];
	break; 
      }
  fg = (fg1-fg2)/(zz[iz]-zz[iz-1])*(redshift-zz[iz-1]) + fg2;
  //if(exp(fg)>1)
  //printf("ERROR: %f %e %f\n",redshift,exp(mass),exp(fg));
  //return 1;
  return(exp(fg));

}

float mstar_lookup(float redshift, float mass)
{
  FILE *fp;
  float fg1,fg2,fg,x1,x2,x3;
  int i,j,iz;
  static int flag = 1;
  static float *zz, *mh, **fgas;
  static int nz, nm;

  if(flag)
    {
      fp = openfile("fgas_lookup.dat");
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

  // NGP in mass;
  for(i=1;i<=nm;++i)
    if(mh[i]>mass || mass>mh[nm]) 
      { 
	fg1 = fgas[iz][i]; 
	fg2 = fgas[iz-1][i];
	break; 
      }
  fg = (fg1-fg2)/(zz[iz]-zz[iz-1])*(redshift-zz[iz-1]) + fg2;
  //printf("LOOKUP %f %e %e %e %e %d %d\n",redshift, exp(mass),exp(fg1),exp(fg2),exp(fg),iz,i);
  return(exp(fg+mass));

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
