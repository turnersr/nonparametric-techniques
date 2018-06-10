#include "utilities.h"

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


int read_double_vector(char* filepath, int vector_size, double* data_vector) {
  
  FILE *input_file = fopen(filepath,"r");

  size_t result;

  if(input_file == NULL) {
    fprintf(stderr, "error opening file\n");
    return -1;
    }

  result = fread(data_vector, sizeof(double), vector_size, input_file);

  if(result != vector_size) {
    fprintf(stderr, "error reading input file\n");
    return -1;
  }

  fclose(input_file);

  return 0; 
}


int write_double_vector(char* filepath, int vector_size, double* data_vector) {
  
  FILE *output_file = fopen(filepath,"w");

  size_t result;

  if(output_file == NULL) {
    fprintf(stderr, "error opening file\n");
    return -1;
    }

  result = fwrite(data_vector, sizeof(double), vector_size, output_file);
  
  if(result != vector_size) {
    fprintf(stderr, "error writing to file\n");
    return -1;
  }

  fclose(output_file);

  return 0; 

}


int next_largest_power2(int x) {

  int p2 = pow(2,ceil(log2(x)));
  return p2;

}
