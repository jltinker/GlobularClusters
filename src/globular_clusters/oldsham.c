#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "header.h"

double g7_nhalo,
  g7_mass,
  g7_mass2,
  g7_msub,
  g7_alpha=0.175,
  g7_mstar,
  g7_redshift,
  g7_dz,
  g7_nmax,
  g7_mbar,
  g7_mhalo;

double func_nhalo(double m);
double subhalo_abundance(double m);
double func_subhalo1(double mhost);
double func_nstar(double m);
double stellar_mass_function(double m);
double func_find_mass(double m);
double get_random_ratio(double nmax);
double func_merger_ratio(double mr);
double n_probability(double n);
double n_mean_mass(double nmax);
double random_circularity_parameter(void);
double tabulated_halo2stellar_mass(double mhalo);
void analytic_stellar_mass_growth(double start_mass, double start_redshift);
double merger_probability(double mass, double redshift);
double func_merger_probability(double m);

double halo2stellar_mass(double m)
{
  static int flag=1, n=100, prev_cosmo=-1, call_count=1;
  static double *mh, *ms, *mx, mlo, mhi, dlogm, mmax;
  int i;
  double a;

  //return 1;

  if(flag || RESET_COSMOLOGY != prev_cosmo)
    {
      prev_cosmo = RESET_COSMOLOGY;
      flag = 0;
      mh = dvector(1,n);
      ms = dvector(1,n);
      mx = dvector(1,n);

      mlo = 1.0e10;
      mhi = 1.0e15;
      mmax = mhi*10;
      dlogm = log(mhi/mlo)/(n-1);
      
      for(i=1;i<=n;++i)
	{
	  mh[i] = exp(dlogm*(i-1))*mlo;
	  g7_nhalo = qromo(func_nhalo,log(mh[i]),log(mmax),midpnt);

	  // get L* stellar mass at z=0;
	  /*
	  RESET_COSMOLOGY++;
	  REDSHIFT = 0;
	  SIGMA_8 = 0.8;

	  fmuh(qromo(func_nhalo,log(1.0e14),log(mmax),midpnt));
	  REDSHIFT = 2;
	  SIGMA_8 = 0.8*growthfactor(2);
	  RESET_COSMOLOGY++;
	  fmuh(qromo(func_nhalo,log(1.0e13),log(mmax),midpnt));


	  g7_nhalo = 3.16e-3;
	  printf("%e %e\n",g7_nhalo,zbrent(func_nstar,log(1.0e-13),log(1.0e13),1.0E-4));
	  exit(0);
	  */
	  //fmuh(g7_nhalo);
	  /*
	  if(func_nstar(log(1.0E-3)>0))
	    {
	      ms[i] = log(1.0E-3);
	      mh[i] = log(mh[i]);
	      printf("CC %e %e\n",g7_nhalo,func_nstar(log(1.0E-3)));
	      continue;
	    }
	  */
	  ms[i] = zbrent(func_nstar,log(1.0e-13),log(1.0e13),1.0E-4);
	  printf("SHAM%d %f %e %e %e %e\n",call_count,REDSHIFT,mh[i],exp(ms[i]),g7_nhalo,exp(ms[i])*0.7/mh[i]/(0.044/0.26));
	  //printf("MF %e %e %e %e\n",mh[i],dndM_interp(mh[i]),mh[i],stellar_mass_function(log(mh[i]))/mh[i]);
	  //printf("CUMU %e %e %e\n",mh[i],g7_nhalo,qromo(stellar_mass_function,log(mh[i]),log(1.0e13),midpnt));
	  //g7_msub = mh[i];
	  //printf("%e\n",qromo(func_subhalo1,log(mh[i]),log(HOD.M_max),midpnt));
	  mh[i] = log(mh[i]);
	}
      //exit(0);
      call_count++;
      spline(mh,ms,n,1.0E+30,1.0E+30,mx);
    }
  m = log(m);
  if(m<mh[1])return 1;
  splint(mh,ms,mx,n,m,&a);
  return exp(a);

}


double tabulated_halo2stellar_mass(double mhalo)
{
  double zhi, zlo, mhi, mlo, dlogm, dz, mh1, ms1, mmax, z, m, m1, m2, w2;
  int i,j;

  static int flag=1, nz=100, nm=100;
  static double **ms, **mh, **zz, *redshift;
  FILE *fp, *outf;
  
  if(flag)
    {
      flag = 0;
      if(!(fp=fopen("tabulated_SHAM.dat","r")))
	{
	  outf = fopen("tabulated_SHAM.dat","w");
	  zhi = 2.5;
	  zlo = 0;
	  dz = (zhi-zlo)/(nz-1);
	  mlo = 1e9;
	  mhi = 1e15;
	  mmax = mhi*10;
	  dlogm = log(mhi/mlo)/(nz-1);
	  for(i=1;i<=nz;++i)
	    {
	      REDSHIFT = (i-1)*dz+zlo;
	      RESET_COSMOLOGY++;
	      SIGMA_8 = 0.8*growthfactor(REDSHIFT);
	      
	      for(j=1;j<=nm;++j)
		{
		  mh1 = exp(dlogm*(j-1))*mlo;
		  g7_nhalo = qromo(func_nhalo,log(mh1),log(mmax),midpnt);
		  ms1 = zbrent(func_nstar,log(1.0e-13),log(1.0e13),1.0E-4);
		  fprintf(outf,"%e %e %e\n",REDSHIFT,mh1,exp(ms1));
		}
	    }
	  fclose(outf);
	  fp = openfile("tabulated_SHAM.dat");
	}
      ms = dmatrix(1,nz,1,nm);
      mh = dmatrix(1,nz,1,nm);
      zz = dmatrix(1,nz,1,nm);
      redshift = dvector(1,nz);

      for(i=1;i<=nz;++i)
	for(j=1;j<=nm;++j)
	  {
	    fscanf(fp,"%lf %lf %lf",&redshift[i],&mh[i][j],&ms[i][j]);
	    ms[i][j] = log(ms[i][j]);
	    mh[i][j] = log(mh[i][j]);
	    redshift[i] = log(1+redshift[i]);
	  }

      for(i=1;i<=nz;++i)
	spline(mh[i],ms[i],nm,1.0E+30,1.0E+30,zz[i]);
    }

  z = log(1+REDSHIFT);
  for(i=2;i<=nz;++i)
    if(redshift[i]>z)break;

  splint(mh[i],ms[i],zz[i],nm,log(mhalo),&m2);
  splint(mh[i-1],ms[i-1],zz[i-1],nm,log(mhalo),&m1);

  w2 = (redshift[i]-z)/(redshift[i]-redshift[i-1]);
  // printf("%e %e %e %f %e %f %f %f %e\n",mhalo,m2,m1,redshift[i],w2,redshift[i],redshift[i-1],redshift[i]-redshift[i-1],exp(((1-w2)*m2 + w2*m1)));
  return exp(((1-w2)*m2 + w2*m1));

}


double func_nhalo(double m)
{
  m = exp(m);
  return (dndM_interp(m)+subhalo_abundance(m))*m;
}

/* integrate over the parent halo mass function
 * to get the density of subhalos of mass m
 */
double subhalo_abundance(double m)
{  
  g7_msub = m;
  return qromo(func_subhalo1,log(m),log(HOD.M_max),midpnt);
}

double func_subhalo1(double mhost)
{
  mhost = exp(mhost);
  return pow(g7_msub/mhost,-0.7)*exp(-9.9*pow(g7_msub/mhost,2.5))*0.3*dndM_interp(mhost)*mhost/g7_msub;

  return pow(g7_msub/mhost,-0.8)*exp(-g7_msub/mhost*1.25)*0.2*dndM_interp(mhost)*mhost/g7_msub;
}
 
/* note, need to integrate 
 * the stellar mass function over log10(m),
 * but m is passed as lnM
 */
double func_nstar(double m)
{
  //return g7_nhalo - qromo(stellar_mass_function,m,log(1.0e14),midpnt);
  return g7_nhalo - qromo(stellar_mass_function,m/log(10),log10(1.0e14),midpnt);
}

/* 
 * this is dN/dlog_10(M). (and there's no 'h' in the stellar mass unit).
 *
 * taking the Marchesini etal SMFs and interpolating as a function of log(1+z).
 * the only parameter that changes to any degree is phi*.
 */
double stellar_mass_function(double m)
{
  static int flag=1, n=4, prev_cosmo=-1;
  static double *redshift, *m1, *zm1, *a1, *za1, *p1, *zp1;
  static double phi, ms, alpha;

  //m = m/log(10);
  
  if(flag)
    {
      flag = 0;
      redshift = dvector(1,n);
      m1 = dvector(1,n);
      zm1 = dvector(1,n);
      a1 = dvector(1,n);
      za1 = dvector(1,n);
      p1 = dvector(1,n);
      zp1 = dvector(1,n);

      redshift[1] = log(1+0.1);
      redshift[2] = log(1+1.6);
      redshift[3] = log(1+2.5);
      redshift[4] = log(1+3.5);

      m1[1] = 10.96;
      m1[2] = 10.91;
      m1[3] = 10.96;
      m1[4] = 11.38;

      a1[1] = -1.18;
      a1[2] = -0.99;
      a1[3] = -1.01;
      a1[4] = -1.39;
      
      p1[1] = 0.003087*2.91; //2.91 = 0.7^-3
      p1[2] = 0.001017*2.91; //2.91 = 0.7^-3
      p1[3] = 0.000395*2.91; //2.91 = 0.7^-3
      p1[4] = 0.000053*2.91; //2.91 = 0.7^-3

      spline(redshift,m1,n,1.0E+30,1.0E+30,zm1);
      spline(redshift,a1,n,1.0E+30,1.0E+30,za1);
      spline(redshift,p1,n,1.0E+30,1.0E+30,zp1);


    }
  if(RESET_COSMOLOGY!=prev_cosmo)
    {
      prev_cosmo = RESET_COSMOLOGY;
      splint(redshift,m1,zm1,n,log(1+REDSHIFT),&ms);
      splint(redshift,p1,zp1,n,log(1+REDSHIFT),&phi);
      splint(redshift,a1,za1,n,log(1+REDSHIFT),&alpha);
      fprintf(stdout,"Resetting SMF: %e %e %e %f\n",ms,phi,alpha,REDSHIFT);

    }

  // charlie's version, z=0
  /*
  alpha = -1.25;
  phi = 0.003*2.91;
  ms = 10.95;
  */
  return log(10.0)*phi*pow(10.0,(m-ms)*(1+alpha))*exp(-(pow(10.0,m-ms)));
}


/* Using the data that Andrew sent to me, figure out the fraction of
 * 3:1 stellar masss mergers as a function of halo mass
 */
void stellar_mass_merger_rate()
{
  int i,j,k,n=100;
  double m, dlogm, mlo, mhi, mstar1, m2, m1, merger_rate;

  mlo = 1e11;
  mhi = 1e14;
  dlogm = log(mhi/mlo)/(n-1);

  for(i=1;i<=n;++i)
    {
      m1 = exp((i-1)*dlogm)*mlo;
      mstar1 = halo2stellar_mass(m1);
      g7_mass = mstar1/3;

      printf("D %e %e\n",m1,mstar1);
      //find halo mass of mstar/3
      m2 = zbrent(func_find_mass,log(m1/100.0),log(m1),1.0E-4);
      m2 = exp(m2);

      //old fit with wrong 2:1 rate
      //merger_rate = 0.046*(pow(m2/m1,-1.1)-1);
      //set y = 10**-1.17/(-0.8)*(1-x**-0.8)
      merger_rate = pow(10.0,-1.17)/0.8*(pow(m2/m1,-0.8)-1);

      //print out merger rate
      printf("MM %e %e %e %e\n",m1,mstar1,m2,merger_rate);

    }
  exit(0);
}

double ncen_by_mergers(double m)
{
  int i,j,k;
  double dlogm, mlo, mhi, mstar1, m2, m1, merger_rate;

  static int flag=1, n=100;
  static double *x, *y, *z;

  if(flag)
    {
      flag = 0;
  
      x = dvector(1,n);
      y = dvector(1,n);
      z = dvector(1,n);

      mlo = 1e11;
      mhi = HOD.M_max;
      dlogm = log(mhi/mlo)/(n-1);
      
      for(i=1;i<=n;++i)
	{
	  m1 = exp((i-1)*dlogm)*mlo;
	  mstar1 = halo2stellar_mass(m1);
	  g7_mass = mstar1/3;
	  
	  m2 = zbrent(func_find_mass,log(m1/100.0),log(m1),1.0E-4);
	  m2 = exp(m2);
	  
	  merger_rate = pow(10.0,-1.17)/0.8*(pow(m2/m1,-0.8)-1);

	  x[i] = log(m1);
	  y[i] = log(merger_rate*3.23);
	  //y[i] = log(merger_rate);
	  if(y[i]>0)y[i]=0;
	}
      spline(x,y,n,1.0E+30,1.0E+30,z);
    }
  splint(x,y,z,n,log(m),&m1);
  return exp(m1);

}

double func_find_mass(double m)
{
  return g7_mass - halo2stellar_mass(exp(m));
}


/* project halo forward from z to 0, and see if growth rate
 * equals the known growth rate for that mass
 */
double func_halo_growth_rate(double alpha)
{
  double m;
  m = g7_mass*exp(alpha*g7_redshift);
  return alpha - 1.15*pow(m/1e14,0.1);
}


/* 
 */
double merger_probability(double mass, double redshift)
{
  double s1,mhi;
  REDSHIFT = redshift;
  SIGMA_8 = 0.8*growthfactor(redshift);
  RESET_COSMOLOGY++;
  g7_mass2 = mass;
  //printf("\n");
  mhi = mass*(g7_nmax+1);
  if(mhi>HOD.M_max)mhi = HOD.M_max;
  s1 = qromo(func_merger_probability,log(mass),log(mhi),midpnt);
  //printf("\n>%e %e %e %e\n",s1,dndM_interp(mass),mass,s1/dndM_interp(mass)/mass);
  return s1/dndM_interp(mass)/mass;
}

double func_merger_probability(double m)
{
  double alpha, dz, ptot;

  g7_mass = exp(m);
  alpha = zbrent(func_halo_growth_rate,0.0,2.0,1.0E-4);
  dz = -log(1+(g7_mbar))/alpha;
  ptot = 9 + (pow(g7_nmax+1,1-g7_alpha) - pow(10.0,1-g7_alpha))/(1-g7_alpha);

  //printf(">> %e %e %e %e %e %e\n",
  //	 g7_mass,dndM_interp(g7_mass),n_probability(g7_mass/g7_mass2)/ptot,g7_mass/g7_mass2,g7_dz,dz);

  return dndM_interp(g7_mass)*n_probability(g7_mass/g7_mass2)/ptot*g7_dz/dz*g7_mass;
}

/* generate a bunch of halos at z=1 or 2 and project them forward.
 */
void galaxy_evolution_wrapper2()
{
  int i,n=0,nhalo=10000;
  double redshift, mass, mmin=12, p, pmax;
  

  redshift = REDSHIFT = 1.0;
  SIGMA_8 = 0.8*growthfactor(REDSHIFT);
  RESET_COSMOLOGY++;

  pmax = dndM_interp(pow(10.0,mmin))*pow(10.0,mmin);

  srand48(1000);

  while(n<nhalo)
    {
      //srand48((int)second()*1000);
      mass = pow(10.0,mmin + drand48()*3);      
      p = dndM_interp(mass)*mass;
      //printf("%e %e\n",mass,p/pmax);
      if(drand48()>p/pmax)continue;
      n++;
      //fmuh(mass);

      //NB!! need to recover this function!!!
      //analytic_stellar_mass_growth(mass,redshift);


    }


  exit(0);

}

/*
 * Analytic merger trees from the van der Wel paper.
 * P(n) = constant. n = m/M.
 */
void analytic_merger_tree(double M0)
{
  int i, j, k;
  double n, nmax=999;
  double M_now, redshift, dz, alpha, msub, M_final, mbar, mstar, mstar_now, m, merger_ratio,sham_mstar;

  redshift = 6.0;

  alpha = 1.15; //1e14 z=0
  alpha = 1.15*pow(M0/1e14,0.1);

  //alpha = 0.6; // ?
  M_now = M0*exp(-alpha*redshift);
  //mbar = (1-g7_alpha)/g7_alpha*(1-pow(nmax+1,-g7_alpha))/(pow(nmax+1,1-g7_alpha)-1);
  mbar = n_mean_mass(nmax+1);
  dz = -log(1+(mbar))/alpha;
  fmuh(mbar);

  srand48(atoi(ARGV[3]));

  while(redshift>0)
    {
      //n = drand48()*nmax+1;
      n = get_random_ratio(nmax);
      msub = M_now/n;
      redshift += dz;
      M_now += msub;
    }
  M_final = M_now;
  fmuh(M_final);

  redshift = 6.0;
  M_now = M0*exp(-alpha*redshift);
  srand48(atoi(ARGV[3]));
  REDSHIFT = redshift;
  RESET_COSMOLOGY++;
  SIGMA_8 = 0.8*growthfactor(REDSHIFT);

  // get original stellar mass
  mstar_now = tabulated_halo2stellar_mass(M_now);

  while(redshift>0)
    {
      //n = drand48()*nmax+1;
      n = get_random_ratio(nmax);
      msub = M_now/n;

      // find the halo mass that's minimum for major merger (assuming sham masses for everything)
      m = M_now;
      merger_ratio = 1;      
      sham_mstar = tabulated_halo2stellar_mass(M_now);
      g7_mstar = sham_mstar;
      g7_mhalo = M_now;
      merger_ratio = zbrent(func_merger_ratio,1.0,1.0E5,1.0E-4);
      m = M_now/merger_ratio;
      /*
      while (merger_ratio<3)
	{
	  m = m/1.01;
	  mstar = halo2stellar_mass(m);
	  merger_ratio = mstar_now/sham_mstar;
	  //printf("%e %e %e %e\n",m,mstar,merger_ratio,M_now/m);
	}
      */
      //exit(0);
      mstar = tabulated_halo2stellar_mass(msub);
      printf("MACC %f %e %e %e %e %e %e %e %e %e %e %e %e %e\n",redshift,M_now,msub,M_now/msub,M_final/msub,mstar_now,mstar,mstar_now/mstar,mstar_now/mstar/(M_now/msub),m,M_now/m,sham_mstar,mstar_now/sham_mstar,sham_mstar/mstar);

      redshift += dz;
      M_now += msub;
      mstar_now += mstar;
      REDSHIFT = redshift;
      RESET_COSMOLOGY++;
      SIGMA_8 = 0.8*growthfactor(REDSHIFT);
    }
  exit(0);

}
double func_merger_ratio(double mr)
{
  return 3 - g7_mstar/tabulated_halo2stellar_mass(g7_mhalo/mr);
}

double get_random_ratio(double nmax)
{
  double pmax = 1.0, n, p;

  p = 0;
  while (drand48()>p)
    {
      n = drand48()*nmax+1;
      p = n_probability(n);
      //p = pow(n,-g7_alpha);
    }
  return n;
}

double n_probability(double n)
{

  if(n<10)return 1;
  return pow(n/10,-g7_alpha);
}

double n_mean_mass(double nmax)
{
  int i;
  double n,s0,s1,dn;

  dn = (nmax-1)/(nmax*10);
  s0 = s1 = 0;
  for(i=1;i<=(int)(nmax*10);++i)
    {
      n = (i-0.5)*dn+1;
      s0 += dn*n_probability(n);
      s1 += dn*n_probability(n)*(1/n);
    }
  return s1/s0;

}

void tree_post_processing()
{
  FILE *fp;
  int nz,i,N_sats=0;
  char aa[1000];
  double M_now, redshift, dz, alpha, msub, M_final, mbar, mstar, mstar_now, m, 
    merger_ratio,sham_mstar,x1,mstar1,sigma_now,sigma2,sigma_FJ, radius_now, eta, eps;
  double eta_circ, t_dyn, t_merge, f_energy=1;

  double FJ_evolution=0.5;

  fp = openfile(ARGV[3]);
  nz = filesize(fp);

  for(i=1;i<=nz;++i)
    {
      fscanf(fp,"%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",aa,&redshift,&M_now,&msub,
	     &x1,&x1,&mstar1,&mstar,&x1,&x1,&m,&x1,&sham_mstar,&x1,&x1);
      if(i==1)
	{
	  mstar_now = mstar1;
	  sigma_now = pow(mstar_now/1e12,0.333)*pow(1+redshift,FJ_evolution)*300;
	  sigma_FJ = sigma_now;
	  radius_now = 1; // 1kpc from bezanson+
	  printf("%f %e %e %e %e %e %e %e %e %e %d\n",redshift,M_now,msub,mstar_now,
		 mstar,mstar1,sham_mstar,sigma_now,radius_now,sigma_FJ,N_sats);
	}
      if(i==nz)
	{
	  printf("%f %e %e %e %e %e %e %e %e %e %d\n",redshift,M_now,msub,mstar_now,
		 mstar,mstar1,sham_mstar,sigma_now,radius_now,sigma_FJ,N_sats);
	  continue;
	}



      //*****************************************
      //this is is PROPER MERGER TIMESCALES MODEL
      // - at each time of accretion, calc t_merge, see if <t_hub(z)

      //random sample circularity parameter
      eta_circ = random_circularity_parameter();
      //hardwire eta_circ
      //eta_circ = 1;
      //get t_dyn
      t_dyn = 1.4*pow(1+redshift,-1.5);
      //get t_merge      
      t_merge = t_dyn*0.216*pow(M_now/msub,1.3)/log(1+M_now/msub)*
	exp(1.9*eta_circ)*(0.5/(1-eta_circ*eta_circ/2));
      //if circularity low, assume immediate merging.
      if(eta_circ<0.2)t_merge = 0;
      // temp for testing
      //t_merge = 0;

      //see if sat will merge before z=0
      if(t_merge>lookback_time(redshift))
	{
	  if(mstar>6.81e10 && redshift>0.25)N_sats++; //keep track of remaining L* sats
	  continue;
	}
      //if it would have merged in between z=0 and z=0.25, incude as sat at z=0.25
      if(mstar>6.81e10 && redshift>0.25)
	if(lookback_time(redshift)-t_merge<lookback_time(0.25))
	  N_sats++;

      //if(drand48()>0.333)continue;

      //NB assume that mass is reduced by factor of 2
      //mstar*=0.07;
      //f_energy = 100;

      sigma2 = pow(mstar/1.0e12,0.333)*pow(1+redshift,FJ_evolution)*300;
      eta = mstar/mstar_now;
      eps = sigma2/sigma_now;
      //if high circularity, assume all energy lost to dynamical friction
      if(t_merge>0.001)
	{
	  radius_now *= (1+eta)*(1+eta)/(1+eta*eps+eta/f_energy);
	  sigma_now *= sqrt((1+eta*eps+eta/f_energy)/(1+eta));
	}
      else //use the maximal values
	{
	  radius_now *= (1+eta)*(1+eta)/(1+eta*eps);
	  sigma_now *= sqrt((1+eta*eps)/(1+eta));
	}
      sigma_FJ = pow(mstar_now/1.0e12,0.333)*pow(1+redshift,FJ_evolution)*300;

      printf("%f %e %e %e %e %e %e %e %e %e %d\n",redshift,M_now,msub,mstar_now,
	     mstar,mstar1,sham_mstar,sigma_now,radius_now,sigma_FJ,N_sats);
      mstar_now += mstar;
      continue;
      //********************************



      //********************************
      //this is for the ALL RADIAL MODEL
      // - so here everything mergers with maximum energy.
      //if(drand48()>0.2)continue;
      sigma2 = pow(mstar/1.0e12,0.333)*pow(1+redshift,FJ_evolution)*300;
      eta = mstar/mstar_now;
      eps = sigma2/sigma_now;
      radius_now *= (1+eta)*(1+eta)/(1+eta*eps);
      sigma_now *= sqrt((1+eta*eps)/(1+eta));
      sigma_FJ = pow(mstar_now/1.0e12,0.333)*pow(1+redshift,0.5)*300;
      printf("%f %e %e %e %e %e %e %e %e %e %d\n",redshift,M_now,msub,mstar_now,
	     mstar,mstar1,sham_mstar,sigma_now,radius_now,sigma_FJ,N_sats);
      mstar_now += mstar;
      continue;
      //********************************





      //********************************
      //this is is ALL CIRCULAR MODEL
      // - so here only things with Mhalo/Msub<10 merge with the central
      // - want to keep track of satellites left over as well.
      // - want to keep track of energy put into halo
      if(i==nz)
	{
	  printf("%f %e %e %e %e %e %e %e %e %e %d\n",redshift,M_now,msub,mstar_now,
		 mstar,mstar1,sham_mstar,sigma_now,radius_now,sigma_FJ,N_sats);
	  continue;
	}
      if(M_now/msub>100000)
	{
	  if(mstar>9.74e9)N_sats++; //keep track of remaining L* sats
	  continue;
	}
      sigma2 = pow(mstar/1.0e12,0.333)*pow(1+redshift,0.5)*300;
      eta = mstar/mstar_now;
      eps = sigma2/sigma_now;
      radius_now *= (1+eta)*(1+eta)/(1+eta*eps+eta);
      sigma_now *= sqrt((1+eta*eps+eta)/(1+eta));
      mstar_now += mstar;
      sigma_FJ = pow(mstar_now/1.0e12,0.333)*pow(1+redshift,0.5)*300;
      M_now += msub;
      
      //added_heat += 

      printf("%f %e %e %e %e %e %e %e %e %e %d\n",redshift,M_now,msub,mstar_now,
	     mstar,mstar1,sham_mstar,sigma_now,radius_now,sigma_FJ,N_sats);
      continue;
      //********************************

      //if Mhost/Msub > 10, no merger.
      if(M_now/msub<=10 || i==1 || i==nz)
	printf("%f %e %e %e %e %e %e\n",redshift,M_now,msub,mstar_now,mstar,mstar1,sham_mstar);
      
    }
  exit(0);

}

double func_calc_lookback_time(double z)
{
  return 1/((1+z)*sqrt(OMEGA_M*(1+z)*(1+z)*(1+z)+(1-OMEGA_M)));
}
//returns t_lookback in Gyr (assumes H0=70)
double lookback_time(double redshift)
{
  return 13.9*qromo(func_calc_lookback_time,0,redshift,midpnt);
}

/* using the distribution function of
 * zentner+05 (eq. 5, a=2.22)
 */
double random_circularity_parameter(void)
{
  double p,pmax,eta;

  pmax = 8.692*pow(0.5*(1-0.5),1.22);
  p = 0;
  while(drand48()>p)
    {
      eta = drand48();
      p = 8.692*pow(0.5*(1-0.5),1.22)/pmax;
    }
  return eta;
}
