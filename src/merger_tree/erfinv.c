#include <stdlib.h>
#include <math.h>

double func1_inv(double x);
double zbrent(double (*func)(double), double x1,double x2, double tol);

double x_g1;

/* input x is between -1 and 1
 */
double erfinv(double x)
{
  x_g1 = x;
  return zbrent(func1_inv,-1000.0,1000.0,1.0E-5);
}

double func1_inv(double x)
{
  return erf(x) - x_g1;
}
