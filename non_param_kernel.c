#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define SIZE 1000000000
#define MAX 1000
#define M_1_SQRT2PI 0.398942280401432677939946059
#define M_1_2 0.5
#define M_3_4 0.75
#define M_70_81 0.86419753086

#define MAD_constant 1.4826

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
  return sqrt(c[0]);

}

double boxcar(double a)
{
  if (fabs(a) <= 1) {
      return M_1_2;
    }
  else {
    return 0;
  }
}

double guassian(double a)
{
  return M_1_SQRT2PI * exp(-1 * M_1_2 * a);
}

double tricube(double a)
{
  if ( fabs(a) <= 1) {
    return M_70_81 * (1 - pow(fabs(a),3));
    }
  else {
    return 0;
     }
}

double epanechnikov(double a)
{
  if ( fabs(a) <= 1 ) {
    return M_3_4 * (1 - pow(a,2));
  }
  else {
    return 0;
  }
}

double local_guassian_kernel_approximate(double p_x ,double p_y, double * restrict x_points, double * restrict y_points, int n_points, double h)
{
  int i;
  double u;
  double k_factor;
  double k_tmp;

  u = 0;
  k_factor = 0;
  h = 0.01;
  for(i = 0; i< n_points; i++) {
    k_tmp = guassian(fabs(p_x - x_points[i]) / h);
    u += y_points[i] * k_tmp;
    k_factor += k_tmp;
  }
  return u / k_factor;
}

void guassian_kernel_regression(double * restrict result, double * restrict x_points, double * restrict y_points, int n_points, double h)
{
  int i;

  for(i = 0; i< n_points; i++) {
    result[i] = local_guassian_kernel_approximate(x_points[i],y_points[i],x_points,y_points,n_points, h);
  }
}

int main(int argc, char* argv[]) {
  int n; // How large of a curve                                                                                                                                                  

  if( argc > 1 ) {
    n = atoi(argv[1]);
  }
  else {
    n = SIZE;
  }


  double *restrict x_points; // Test data                                                                                                                                         
  double *restrict y_points; // Test data                                                                                                                                         
  double *restrict c; // The result, m_hats, wiil go ehre                                                                                                                         
  double h; // Bandwidth                                                                                                                                                          
  h = .001;

  x_points =  (double*)malloc(n*sizeof(double));
  y_points =  (double*)malloc(n*sizeof(double));
  c =  (double*)malloc(n*sizeof(double));

  // Create a linear test set                                                                                                                                                     
  for( int i = 0; i < n; ++i ) x_points[i] = i;
  for( int i = 0; i < n; ++i ) y_points[i] = i;

  guassian_kernel_regression(c, x_points, y_points, n, h);

  printf("Done Testing\n");
  return 0;

}
