#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void min_max(const double *array, size_t length, double* min, double* max);
int read_double_vector(char* filepath, int vector_size, double* data_vector);
int write_double_vector(char* filepath, int vector_size, double* data_vector);
int next_largest_power2(int x); 
