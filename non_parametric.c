#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "libfft.h"
#include "kernels.h"
#include "binning.h"
#include "utilities.h"

int main(int argc, char* argv[]) {

  int return_code;
  
  int n; // How large of an input curve
  double h; // Bandwidth

  double min;
  double max;

  double *restrict univariate_input; // Input data

  double* discretized_vector;
  double* rfft_vector;
  double* irfft_vector;
  double* smoothed_data;


  
  char* input_file;
  char* output_file;

  if(argc == 5) {
    n = atoi(argv[1]);
    h = atof(argv[2]);
    input_file = argv[3];
    output_file = argv[4];
  }
  else {
    fprintf(stderr, "error reading parameters. Call the program as ./non_parametric NPOINTS BANDWIDTH INPUT OUTPUT\n");
    return -1;
  }


  printf("Expecting %d packed doubles in %s and writing density estimation with bandwidth %f to %s\n",
	 n, input_file, h, output_file);


  int gridsize = determine_gridsize(n);
  
  univariate_input = (double*) malloc(n *sizeof(double));
  discretized_vector = (double *)calloc(gridsize, sizeof(double));
  rfft_vector = (double *)calloc(gridsize, sizeof(double));
  irfft_vector = (double *)calloc(gridsize, sizeof(double));
  smoothed_data = (double *)calloc(gridsize, sizeof(double));


  
  if (univariate_input == NULL) {
    fprintf(stderr, "error creating univariate input array");
    return -1;
  }
  
  if (discretized_vector == NULL) {
    fprintf(stderr, "error creating univariate input array");
    return -1;
  }
  
  if (rfft_vector == NULL) {
    fprintf(stderr, "error creating rfft_vector");
    return -1;
  }
  
  if (irfft_vector == NULL) {
    fprintf(stderr, "error creating irfft_vector");
    return -1;
  }
  
  if (smoothed_data == NULL) {
    fprintf(stderr, "error creating smoothed_data");
    return -1;
  }

    
  
  return_code = read_double_vector(input_file, n, univariate_input);

  if(return_code == -1) {
    return -1;
  }

  min_max(univariate_input, n, &min, &max);

  int cut = 3; 
  


  double a = min - cut * h;
  double b = max + cut * h;
  double RANGE = b-a;

  //  printf("%f, %f, %f, %f, %d, %f, %d\n", a, b, min, max, cut, h, gridsize);
  
  linear_binning(discretized_vector, univariate_input, n, a, b, gridsize);
  RFFT(discretized_vector, gridsize, rfft_vector);
  guassian_fft_kde(rfft_vector, gridsize, h,  gridsize, RANGE, smoothed_data);
  IRFFT(smoothed_data, gridsize, irfft_vector);
  
  return_code = write_double_vector(output_file, gridsize, irfft_vector);

  if(return_code == -1) {
    return -1;
  }


  printf("wrote %d points of the estimate density function to %s\n", gridsize, output_file);
  printf("KDE Domain Parameters\n");
  printf("import numpy \n");
  printf("grid, delta = numpy.linspace(%f, %f, %d, retstep=True)\n", a, b, gridsize);
  
  free(rfft_vector);
  free(irfft_vector);
  free(discretized_vector);
  return 0;

}
