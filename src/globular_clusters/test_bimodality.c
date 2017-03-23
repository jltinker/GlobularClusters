#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cutil.h"
#include "nrutil.h"

void sort(unsigned long n, float arr[]);

void test_bimodality(float *x, float n)
{
  float mu1, mu2, s1, s2, p, mlo, mhi, slo, shi, dm, ds, dp;
  int i, j, k, nx=10;

  // sort the distribution for cumulative statistics
  sort(x,n);

  // use the ranges to set the ranges of the means and dispersions
  mlo = x[1];
  mhi = x[n];
  shi = 0.5*(mhi-mlo);
  slo = shi/100;

  dm = (mhi-mlo)/(nx-1);
  ds = (shi-slo)/(nx-1);
  dp = 1.0/(nx-1);

  f

}
