#include "linear_algebra.h"

void difference_vector(double * restrict a, double * restrict b, double * restrict c, int size)
{
  int i;

  double *x = __builtin_assume_aligned(a, 16);
  double *y = __builtin_assume_aligned(b, 16);
  double *z = __builtin_assume_aligned(c, 16);

  for (i = 0; i < size; i++)
    {
      x[i] = fabs(y[i] - z[i]);
    }
}

void l2_distance_loop(double * restrict a, double * restrict b, double * restrict c, int size)
{
  int i;

  double *x = __builtin_assume_aligned(a, 16);
  double *y = __builtin_assume_aligned(b, 16);
  double *z = __builtin_assume_aligned(c, 16);

  for (i = 0; i < size; i++)
    {
      x[0] += (y[i] - z[i]) * (y[i] - z[i]);
    }
}


double l2_distance(double * restrict a, double * restrict b, int size)
{
  double * restrict c;
  c = (double*)malloc(sizeof(double));
  l2_distance_loop(c,a,b,size);
  return sqrt(c[0]);
  free(c);
}

void norm_inner_product_loop(double * restrict a, double * restrict b, int size)
{
  int i;

  double *x = __builtin_assume_aligned(a, 16);
  double *y = __builtin_assume_aligned(b, 16);

  for (i = 0; i < size; i++)
    {
      x[0] += (y[i] * y[i]);
    }
}

double norm(double * restrict a, int size)
{
  double * restrict c;
  c = (double*)malloc(sizeof(double));

  norm_inner_product_loop(c,a,size);
  free(c);
  return sqrt(c[0]);
}
