/*  Gaussian Mixture Modeling
    -------------------------

    This routine calculates probabilities of a one-dimensional data
    set being modeled as a sum of gaussian distributions with
    different means and variances.  The probabilities of a
    given data point belonging to each of the gaussian modes are
    calculated using the Gaussian Mixture Modeling (GMM) method, as
    described in Numerical Recipes, 3rd edition, section 16.1.  It is
    similar to the KMM method adapted by Ashman et al. (1994),
    http://adsabs.harvard.edu/abs/1994AJ....108.2348A

    GMM finds the means and variances of the gaussians and their
    mixing proportions by maximizing the likelihood of the data given
    the parametric model.  The desired number of modes must be
    specified on input.

    How to decide if the data favor a multimodal distribution to a
    unimodal distribution?  Ashman et al. suggest that the log of the
    ratio of the likelihoods, 2 ln(L2/L1), approximately obeys a chi2
    distribution with a number of degrees of freedom equal to "twice
    the difference between the number of parameters of the two models,
    not including the mixing proportions".  Alternatively, it can be
    done using bootstrap of a unimodal distribution with the same
    mean, variance, and number of points as the data, then repeating
    the GMM analysis, and calculating how likely the observed ratio
    L2/L1 appears in the bootstrap sample.

    Alternatively, the Odds Ratio test can be used, where the
    likelihoods of each hypothesis (unimodal vs. multimodal) are
    integrated over all values of the parameters, 2 and 3*nk-1,
    respectively, with the assumed prior probabilities for the
    parameters.  The multidimensional integration is done using Monte
    Carlo methods.  I assume constant priors of each parameter over a
    finite range of parameter values that are most likely given the
    input sample means and variances.  Note that the result depends on
    how the intervals are defined.

    Additional interesting ideas are discussed in Nemec and Nemec (1991),
    http://adsabs.harvard.edu/abs/1991PASP..103...95N

    If you use this code for a publication, please acknowledge the original paper 
    A. L. Muratov & O. Y. Gnedin, 2010, ApJ, 718, 1266
    "Modeling the Metallicity Distribution of Globular Clusters"

    Written by Oleg Gnedin on May 4, 2009
    Last modified on March 19, 2012
*/

#include "inc.h"

#define lnsqrt2pi 0.9189385  // = ln(sqrt(2*pi))

Int n;
VecDoub z(99999);


// statistics function: mean, standard deviation, and kurtosis of an array
Doub fmean( Doub *x, Int nd, Doub &sig, Doub &kurt ) {
  Doub mean=0., s=0., ku=0.;
  for(Int i=0; i<nd; i++) mean += x[i]/nd;
  for(Int i=0; i<nd; i++) s += SQR(x[i]-mean)/(nd-1.);
  for(Int i=0; i<nd; i++) ku += pow(x[i]-mean,4)/nd;
  sig = sqrt(s);
  kurt = ku/(s*s) - 3.;
  return mean;
}



int main( int argc, char *argv[] ) {
  Int i, k, nk, iter, Var, Ndof;
  Doub d, mean, sig, sig2, kurt, chi2, chi2m, chi2v, p, s, ns, DD, DDv,
    _mean[2], _dbmean[2], _sig[2], _dbsig[2], frac2, dbfrac, dbDD, pchi2, pDD, pku,
    sigmax, sigfrac=0.0625;
  FILE *in, *out;

  if(argc < 4) { 
    fprintf(stderr, "Syntax: gmm file var peak1 peak2 ...\n"); return 1; }

  printf("Gaussian Mixture Model of a univariate sample\n");
  printf("looking for %d peaks with ", argc-3);

  // common variance for all peaks if Var=1
  Var = atoi(argv[2]);
  if(Var) printf("the same variance\n"); else printf("different variances\n");

  // read data file
  if( (in = fopen( argv[1], "r" )) == NULL )
      { printf("Can't open input file <%s>\n", argv[1]); return 1; }
  i=0;
  while(!feof(in)) {
    fscanf(in, "%le", &z[i]); i++;
  }
  fclose(in);
  n = i-1;  if(n > 99999) { printf("increase array size\n"); return 1; }

  // test dependence of the significance on the number of input points
  if(0) {
    for(k=1; k<2; k++) 
      for(i=0; i<n; i++) z[i+k*n] = z[i];
    n *= 2;
  }

  // calculate kurtosis of the data: 
  // a unimodal Gaussian has kurt=0, two Gaussians tend to have kurt<0
  mean = fmean(&z[0], n, sig, kurt);
  printf("number of data points = %d  kurtosis = %5.3f\n", n, kurt);
 
  // define data structures
  nk = argc-3;
  MatDoub zz(n,1), zmean(nk,1), zmean1(1,1);

  // fill data matrix
  for(i=0; i<n; i++) zz[i][0] = z[i];

  // set initial guesses for mode peaks
  for(k=0; k<nk; k++) zmean[k][0] = atof(argv[k+3]);
  zmean1[0][0] = zmean[0][0];

  // initialize Gaussian Mixture Model structures
  Gaumixmod gmm(zz, zmean);   // unequal means and variances
  Gaumixmod gmmv(zz, zmean);  // equal variances
  Gaumixmod gmmm(zz, zmean);  // equal means
  Gaumixmod gmm1(zz, zmean1); // unimodal distribution

  // calculate unimodal distribution
  printf("...running unimodal Gaussian\n");
  d=9.; iter=0;

  while(fabs(d) > 1.e-6 && iter < 9999) {
    d = gmm1.estep();
    gmm1.mstep();
    iter++;
    if(fabs(d) <= 1.e-6 || iter>=9999) {
      printf("iter=%d err=%7.1e:  ", iter, d);
      printf("peak=%6.3f (n=%4.1f sig=%5.3f)  ", 
	     gmm1.means[0][0], n*gmm1.frac[0], sqrt(gmm1.sig[0][0][0]));
      printf("logL1=%g\n", gmm1.loglike-n*lnsqrt2pi);
    }
  }


  // run Gaussian mixture model until convergence
  printf("...running Gaussian mixture with different variances\n");
  d=9.; iter=0;

  while(fabs(d) > 1.e-6 && iter < 9999) {
    d = gmm.estep();
    gmm.mstep();
    iter++;

    // constraint on deviation of the variances, see Lo (2008)
    sigmax = gmm.sig[0][0][0];
    for(k=1; k<nk; k++) 
      if(gmm.sig[k][0][0] > sigmax) sigmax = gmm.sig[k][0][0];
    for(k=0; k<nk; k++) 
      if(gmm.sig[k][0][0] < sigfrac*sigmax) gmm.sig[k][0][0] = sigfrac*sigmax;

    if(fabs(d) <= 1.e-6 || iter>=9999) {
      printf("iter=%d err=%7.1e:  ", iter, d);
      for(k=0; k<nk; k++)
	printf("peak%d=%6.3f (n=%4.1f sig=%5.3f)  ", 
	       k+1, gmm.means[k][0], n*gmm.frac[k], sqrt(gmm.sig[k][0][0]));
      printf("logL=%g\n", gmm.loglike-n*lnsqrt2pi);
      if(nk==2) {frac2=gmm.frac[1]; for(k=0; k<2; k++) {_mean[k]=gmm.means[k][0]; _sig[k]=sqrt(gmm.sig[k][0][0]);}}
    }
  }

  // Analytical chi-square statistic
  chi2 = 2.*(gmm.loglike-gmm1.loglike);
  Ndof = 4*(nk-1);  // see McLachlan (1987)
  p = 1. - Chisqdist(Ndof).cdf(chi2);
  printf("Chi-square statistic (null=unimodal): chi2=%4.2f Ndof=%d p=%8.2e\n", chi2, Ndof, p);

  // Separation of the peaks
  DD = fabs(gmm.means[1][0]-gmm.means[0][0])/sqrt((gmm.sig[1][0][0]+gmm.sig[0][0][0])/2.);
  printf("Peak separation DD = %4.2f\n", DD);


  printf("TEST %f %f %f\n",DD,gmm.means[0][0], gmm.means[1][0]);
  exit(0);

  // run Gaussian mixture with same variances
  printf("...running Gaussian mixture with same variances\n");
  d=9.; iter=0;

  while(fabs(d) > 1.e-6 && iter < 9999) {
    d = gmmv.estep();
    gmmv.mstep();
    iter++;

    for(s=0., k=0; k<nk; k++) s += gmmv.sig[k][0][0]/nk;
    for(k=0; k<nk; k++) gmmv.sig[k][0][0] = s;

    if(fabs(d) <= 1.e-6 || iter>=9999) {
      printf("iter=%d err=%7.1e:  ", iter, d);
      for(k=0; k<nk; k++)
	printf("peak%d=%6.3f (n=%4.1f sig=%5.3f)  ", 
	       k+1, gmmv.means[k][0], n*gmmv.frac[k], sqrt(gmmv.sig[k][0][0]));
      printf("logL=%g\n", gmmv.loglike-n*lnsqrt2pi);
    }
  }
  chi2v = 2.*(gmmv.loglike-gmm1.loglike);
  Ndof = 2*(nk-1);  // see McLachlan (1987)
  if(chi2v>0) p = 1. - Chisqdist(Ndof).cdf(chi2v); else p=1;
  printf("Chi-square statistic (null=unimodal): chi2=%4.2f Ndof=%d p=%8.2e\n", chi2v, Ndof, p);
  DDv = fabs(gmmv.means[1][0]-gmmv.means[0][0])/sqrt((gmmv.sig[1][0][0]+gmmv.sig[0][0][0])/2.);
  printf("Peak separation DD = %4.2f\n", DDv);
  

  // run Gaussian mixture with same means
  printf("...running Gaussian mixture with same means and different variances\n");
  d=9.; iter=0;

  while(fabs(d) > 1.e-6 && iter < 9999) {
    d = gmmm.estep();
    gmmm.mstep();
    if(!iter && nk>=2) { gmmm.sig[0][0][0] *= 2.; gmmm.sig[1][0][0] /= 2.; }
    iter++;
    
    for(s=0., k=0; k<nk; k++) s += gmmm.means[k][0]/nk;
    for(k=0; k<nk; k++) gmmm.means[k][0] = s;
    // constraint on deviation of the variances, see Lo (2008)
    sigmax = gmmm.sig[0][0][0];
    for(k=1; k<nk; k++) 
      if(gmmm.sig[k][0][0] > sigmax) sigmax = gmmm.sig[k][0][0];
    for(k=0; k<nk; k++) 
      if(gmmm.sig[k][0][0] < sigfrac*sigmax) gmmm.sig[k][0][0] = sigfrac*sigmax;

    if(fabs(d) <= 1.e-6 || iter>=9999) {
      printf("iter=%d err=%7.1e:  ", iter, d);
      for(k=0; k<nk; k++)
	printf("peak%d=%6.3f (n=%4.1f sig=%5.3f)  ", 
	       k+1, gmmm.means[k][0], n*gmmm.frac[k], sqrt(gmmm.sig[k][0][0]));
      printf("logL=%g\n", gmmm.loglike-n*lnsqrt2pi);
    }
  }
  chi2m = 2.*(gmm.loglike-gmmm.loglike);
  Ndof = 1;
  p = 1. - Chisqdist(Ndof).cdf(chi2m);
  printf("Chi-square statistic (null=equal means): chi2=%4.2f Ndof=%d p=%8.2e\n", chi2m, Ndof, p);

  chi2m = 2.*(gmm.loglike-gmmv.loglike);
  Ndof = 1;
  p = 1. - Chisqdist(Ndof).cdf(chi2m);
  printf("Chi-square statistic (null=equal variances): chi2=%4.2f Ndof=%d p=%8.2e\n", chi2m, Ndof, p);


  // Output parameters and probabilities
  if(Var) gmm = gmmv;

  if( (out = fopen( "peakprob.out", "w" )) == NULL )
      { printf("Can't open output file\n"); return 1; }
  // output peak values and fraction of data points in each
  for(k=0; k<nk; k++) 
    fprintf(out, "%6.3f %6.4f %6.4f  ", 
	    gmm.means[k][0], sqrt(gmm.sig[k][0][0]), gmm.frac[k]);
  fprintf(out, "\n");
  // output probabilities of each data point assigned to a given peak
  for(i=0; i<n; i++) {
    for(k=0; k<nk; k++) fprintf(out, "%5.3f ", gmm.resp[i][k]);
    fprintf(out, " %g\n", z[i]);
  }
  fclose(out);
  

  // Various tests of a unimodal vs. a multimodal distribution

  if(nk>1) {

    // Non-parametric Bootstrap to estimate errors of the parameters
    Int ib, decompok=1, nb=100;
    Doub bmean1_,dbmean1,bsig1_,dbsig1,bmean_,dbmean,bsig_,dbsig,bfrac_,bDD_,ku;
    Doub bmean1[nb], bsig1[nb], bDD[nb];
    MatDoub bmean(nk,nb), bsig(nk,nb), bfrac(nk,nb);
    Ran ran(1001);

    printf("...running bootstrap to estimate errors of best-fit parameters\n");

    for(ib=0; ib<nb; ib++) {
      for(i=0; i<n; i++) zz[i][0] = z[(int)(n*ran.doub())];

      Gaumixmod gmm(zz, zmean);
      Gaumixmod gmm1(zz, zmean1);

      d=9.; iter=0;
      while(fabs(d) > 1.e-6 && iter < 9999) {
	d = gmm1.estep(); gmm1.mstep(); iter++;
      }   
      d=9.; iter=0;
      while(fabs(d) > 1.e-6 && iter < 9999 && decompok) {
	d = gmm.estep(); gmm.mstep(); iter++;
	if(Var) {
	  for(s=0., k=0; k<nk; k++) s += gmm.sig[k][0][0]/nk;
	  for(k=0; k<nk; k++) gmm.sig[k][0][0] = s;
	}
	// constraint on deviation of the variances, see Lo (2008)
	sigmax = gmm.sig[0][0][0];
	for(k=1; k<nk; k++) 
	  if(gmm.sig[k][0][0] > sigmax) sigmax = gmm.sig[k][0][0];
	for(k=0; k<nk; k++) 
	  if(gmm.sig[k][0][0] < sigfrac*sigmax) gmm.sig[k][0][0] = sigfrac*sigmax;
	// check for bad matrices
	for(k=0; k<nk; k++) if(gmm.sig[k][0][0] < 1.e-10) decompok=0;
      }
      if(!decompok) { decompok=1; ib--; continue; }
      bmean1[ib] = gmm1.means[0][0];
      bsig1[ib] = sqrt(gmm1.sig[0][0][0]);
      for(k=0; k<nk; k++) {
	bmean[k][ib] = gmm.means[k][0];
	bsig[k][ib] = sqrt(gmm.sig[k][0][0]);
	bfrac[k][ib] = gmm.frac[k];
      }
      bDD[ib] = fabs(gmm.means[1][0]-gmm.means[0][0])/sqrt((gmm.sig[1][0][0]+gmm.sig[0][0][0])/2.);
    }
    bmean1_ = fmean(bmean1, nb, dbmean1, ku);
    bsig1_ = fmean(bsig1, nb, dbsig1, ku);
    printf("Bootstrap unimodal: mean = %5.3f +- %5.3f  sig = %5.3f +- %5.3f\n", 
	   bmean1_, dbmean1, bsig1_, dbsig1);
    for(k=0; k<nk; k++) {
      bmean_ = fmean(&bmean[k][0], nb, dbmean, ku);
      bsig_ = fmean(&bsig[k][0], nb, dbsig, ku);
      bfrac_ = fmean(&bfrac[k][0], nb, dbfrac, ku);
      printf("Bootstrap peak%d: mean = %5.3f +- %5.3f  sig = %5.3f +- %5.3f  n = %3.1f +- %3.1f\n", 
	     k+1, bmean_, dbmean, bsig_, dbsig, n*bfrac_, n*dbfrac);
      if(nk==2) {_dbmean[k]=dbmean; _dbsig[k]=dbsig;}
    }
    bDD_ = fmean(bDD, nb, dbDD, ku);
    printf("Bootstrap DD = %4.2f +- %4.2f\n", bDD_, dbDD);

    // Parametric Bootstrap of a unimodal Gaussian sample 
    // with the same mean and variance as the data:
    // a test of whether a multimodal distribution is an improvement over a unimodal
    Int ichi2=0, iDD=0, iku=0, nbb;
    Doub bchi2;
    Normaldev gau(mean,sig,1001);

    printf("...running parametric bootstrap to rule out unimodal distribution\n");
    nbb = 1000;
    if(n > 1000) nbb = 100; // to prevent long runtime
    if(n < 20) nbb = 40;    // to prevent numerical failures
    decompok=1;
    for(ib=0; ib<nbb; ib++) {
      for(i=0; i<n; i++) zz[i][0] = gau.dev();
      Gaumixmod gmm(zz, zmean);
      Gaumixmod gmm1(zz, zmean1);
      d=9.; iter=0;
      while(fabs(d) > 1.e-6 && iter < 9999) {
	d = gmm1.estep(); gmm1.mstep(); iter++;
      }

      d=9.; iter=0;
      while(fabs(d) > 1.e-6 && iter < 9999 && decompok) {
	d = gmm.estep(); gmm.mstep(); iter++;
	if(Var) {
	  for(s=0., k=0; k<nk; k++) s += gmm.sig[k][0][0]/nk;
	  for(k=0; k<nk; k++) gmm.sig[k][0][0] = s;
	}

	// constraint on deviation of the variances, see Lo (2008)
	sigmax = gmm.sig[0][0][0];
	for(k=1; k<nk; k++) 
	  if(gmm.sig[k][0][0] > sigmax) sigmax = gmm.sig[k][0][0];
	for(k=0; k<nk; k++) 
	  if(gmm.sig[k][0][0] < sigfrac*sigmax) gmm.sig[k][0][0] = sigfrac*sigmax;
	// check for bad matrices
	for(k=0; k<nk; k++) if(gmm.sig[k][0][0] < 1.e-10) decompok=0;
      }
      if(decompok==0) { decompok=1; ib--; continue; }
      // likelihood ratio 
      bchi2 = 2.*(gmm.loglike-gmm1.loglike);
      if(chi2 <= bchi2) ichi2++;
      // peak separation
      bDD_ = fabs(gmm.means[1][0]-gmm.means[0][0])/sqrt((gmm.sig[1][0][0]+gmm.sig[0][0][0])/2.);
      if(DD <= bDD_) iDD++;
      // kurtosis
      bmean_ = fmean(&zz[0][0], n, dbmean, ku);
      if(kurt > ku) iku++;
    }
    if(ichi2)
      printf("Parametric bootstrap: p(chi2) = %g\n", pchi2=(double)ichi2/(double)nbb);
    else 
      printf("Parametric bootstrap: p(chi2) < %g\n", pchi2=1./(double)nbb);
    if(iDD)
      printf("Parametric bootstrap: p(DD) = %g\n", pDD=(double)iDD/(double)nbb);
    else 
      printf("Parametric bootstrap: p(DD) < %g\n", pDD=1./(double)nbb);
    if(iku)
      printf("Parametric bootstrap: p(kurt) = %g\n", pku=(double)iku/(double)nbb);
    else 
      printf("Parametric bootstrap: p(kurt) < %g\n", pku=1./(double)nbb);
  }

  if(nk==2) 
    printf("summary %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %4d %5.3f %5.3f %4.2f %4.2f %5.3f %5.3f %5.3f %s\n", 
	   _mean[0],_dbmean[0],_mean[1],_dbmean[1],_sig[0],_dbsig[0],_sig[1],_dbsig[1],
	   n, frac2, dbfrac, DD, dbDD, pchi2, pDD, pku, argv[1]);

  return 0;
}
