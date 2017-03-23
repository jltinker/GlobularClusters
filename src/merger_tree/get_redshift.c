#include <stdlib.h>
#include <math.h>

double func1_z(double x);
double zbrent(double (*func)(double), double x1,double x2, double tol);

double x_g2;

/* input x is between -1 and 1
 */
double get_redshift(double w)
{
  x_g2 = w;
  return zbrent(func1_z,0.0,10.0,1.0E-5);
}

double func1_z(double z)
{
  return 1.260*(1 + z + 0.09/(1+z) + 0.24*exp(-1.16*z)) - x_g2;
}
