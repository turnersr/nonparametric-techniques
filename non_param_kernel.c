#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "fft.h"

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


void linear_binning(double * gcnts, double * X, int n_points, double a, double b, int M) {

  double* lxi = (double *)calloc(n_points, sizeof(double));
  int* li = (int *)calloc(n_points, sizeof(int));
  double* rem = (double *)calloc(n_points, sizeof(double));

  double delta = (b - a) / (M-1);

  double normalization_factor = delta * n_points;
    
  int li_i;

  for(int i = 0; i < n_points; ++i) {
    lxi[i] = (X[i] - a) / delta;
    li[i] = floor(lxi[i]);
    rem[i] = lxi[i] - li[i];
    li_i = li[i];

    if(li_i > 1 && li_i < M) {
      gcnts[li_i] = (gcnts[li_i] + 1 - rem[i]);
      gcnts[li_i + 1] = (gcnts[li_i + 1] + rem[i]);
    }

  }

  for(int i = 0; i < M; ++i) {
    gcnts[i] /=  normalization_factor;
  }

  free(lxi);
  free(li);
  free(rem);
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


void min_max(const double *array, size_t length, double* min, double* max) {
    // returns the minimum and max value of array
    size_t i;
    double minimum = array[0];
    double maximum = array[0];
    
    for (i = 1; i < length; ++i) {
        if (minimum > array[i]) {
            minimum = array[i];
        }
	
	if (maximum < array[i]) {
            maximum = array[i];
        }

    }
    *min = minimum;
    *max = maximum;
}





int main(int argc, char* argv[]) {
  int n; // How large of a curve

  double min;
  double max; 

  if( argc > 1 ) {
    n = atoi(argv[1]);
  }
  else {
    n = SIZE;
  }

  FILE *fp = fopen("smoothed_data.bin","wb");

  FILE *test_data = fopen("test_data.bin","rb");

  FILE *binned_data = fopen("binned_data.bin","wb");
  

  if(fp == NULL) {
        printf("error creating file\n");
        return -1;
    }

  double *restrict x_points; // Test data                                                                                                                                         
  double *restrict y_points; // Test data                                                                                                                                         
  double *restrict c; // The result, m_hats, wiil go ehre                                                                                                                         
  double h; // Bandwidth                                                                                                                                                          
  h = .1;

  x_points =  (double*)malloc(n*sizeof(double));
  y_points =  (double*)malloc(n*sizeof(double));
  c =  (double*)malloc(n*sizeof(double));

  double* test_points = (double*)malloc(10000*sizeof(double));
  
  fread(test_points, sizeof(double), 10000, test_data);

  fclose(test_data);

  // Create a linear test set                                                                                                                                                     
  for( int i = 0; i < n; ++i ) x_points[i] = i;
  for( int i = 0; i < n; ++i ) y_points[i] = i;

  printf("Start Testing\n");

  min_max(test_points, 10000, &min, &max);

  int cut = 3; 

  int gridsize = 16384;

  double a = min - cut * h;
  double b = max + cut * h;

  double* gcnts = (double *)calloc(gridsize, sizeof(double));

  double* rfft_vector = (double *)calloc(gridsize, sizeof(double));

  printf("%f, %f, %f, %f, %d, %f, %d\n", a, b, min, max, cut, h, gridsize);
  linear_binning(gcnts, test_points, 10000, a, b, gridsize);

  //fwrite(gcnts, sizeof(double), 512, binned_data);
  //fclose(binned_data);

  printf("POST Binnning\n");

  for(int i = 0; i < 512; ++i) { 
      printf("%f\n", gcnts[i]); 
   } 
  
  RFFT(gcnts, gridsize, rfft_vector);

  printf("RFFT Done\n");

  printf("------------------------\n");

 for(int i = 0; i < 512; ++i) {
   printf("%f\n", rfft_vector[i]);
   }

  printf("------------------------\n");
  
  printf("Done Linear Binning\n");

  double RANGE = b-a;

  double* smoothed_data = (double *)calloc(gridsize, sizeof(double));
  guassian_fft_kde(rfft_vector, gridsize, h,  gridsize, RANGE, smoothed_data);

   /* for(int i = 0; i < 512; ++i) {   */
   /*   printf("%f\n", smoothed_data[i]);   */
   /* }  */

  int test_size = 512;
  double* test_in =  (double*)malloc(test_size*sizeof(double));
  double* test_fft_out =  (double*)malloc(test_size*sizeof(double));
  double* test_ifft_out =  (double*)malloc(gridsize*sizeof(double));

  /* for(int i = 0; i < test_size; ++i) { */
  /*   test_in[i] = 2.0; */
  /* } */

  /* RFFT(test_in, test_size, test_fft_out); */

  /* RIFFT(test_fft_out, test_size, test_ifft_out); */
  /* printf("------------------------\n"); */
  
  /* for(int i = 0; i < test_size; ++i) { */
  /*   printf("%d %f, %f\n", i, test_in[i], test_ifft_out[i]); */
  /* } */


  RIFFT(smoothed_data, gridsize, test_ifft_out); 
  
  //  for(int i = 0; i < test_size; ++i) {
  // printf("%f\n", test_ifft_out[i]);
  ///}


  printf("Done Testing\n");

  //  guassian_kernel_regression(c, x_points, y_points, n, h);


  fwrite(test_ifft_out, sizeof(double), gridsize, binned_data);
  
  fclose(fp);
  fclose(test_data);
  fclose(binned_data);

  free(gcnts);
  free(rfft_vector);
  free(x_points);
  free(y_points);
  free(c);
  return 0;

}
