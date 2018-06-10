#include <stdlib.h>
#include <math.h>
#include <stdio.h>


void difference_vector(double * restrict a, double * restrict b, double * restrict c, int size);
void l2_distance_loop(double * restrict a, double * restrict b, double * restrict c, int size);
double l2_distance(double * restrict a, double * restrict b, int size);
void norm_inner_product_loop(double * restrict a, double * restrict b, int size);
double norm(double * restrict a, int size);
