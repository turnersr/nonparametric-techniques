#include "kernels.h"

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

void guassian_fft_kde(double *freq_array, size_t length, double h, int gridsize, double RANGE, double *smoothed_output){

  double FAC1 = 2*pow(((M_PI * h)/RANGE), 2);

  int loopsize = (gridsize/2);// + 1;

  double BC;
  double JFAC;
  double FAC;

  double N1 = (1.0 / gridsize) * M_PI;
  
  for(int j = 0; j < loopsize; ++j) {
    
    JFAC = pow(j,2) * FAC1;
    BC = 1.0 - (1.0 / 3.0) * pow(j*N1,2);
    FAC = exp(-1.0 * JFAC) / BC;
    
    smoothed_output[j] = FAC * freq_array[j];
    smoothed_output[loopsize+j] = FAC * freq_array[loopsize+j];
    
  }
}
