#include "binning.h"
#include "utilities.h"

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


int determine_gridsize(int n){

  int default_value = 512;

  int gridsize;

  if(n > default_value) {

    gridsize = n;

  }
  else {
    gridsize = default_value;
  }

  gridsize = next_largest_power2(gridsize);

  return gridsize;
}
