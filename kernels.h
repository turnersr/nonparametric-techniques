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

double boxcar(double a);
double guassian(double a);
double tricube(double a);
double epanechnikov(double a);

double local_guassian_kernel_approximate(double p_x ,double p_y, double * restrict x_points, double * restrict y_points, int n_points, double h);

void guassian_kernel_regression(double * restrict result, double * restrict x_points, double * restrict y_points, int n_points, double h);

void guassian_fft_kde(double *freq_array, size_t length, double h, int gridsize, double RANGE, double *smoothed_output);
