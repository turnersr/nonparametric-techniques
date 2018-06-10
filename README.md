# nonparametric 

This library implements several methods from nonparametric statistics. Currently can compute kernel density estimator based on the FFT method. It can also compute nonparametric regressions for several kernels using the simple computationally expensive approach. Future plans will include spatial data structures for nearest neighbor methods as well as incorporating numerical integration to compute probabilities. 


# Getting started

The current command line tool will compute the KDE using the Gaussian kernel using the linear binning described in Fan, J. and Marron, J. S. (1994), “Fast Implementations of Nonparametric Curve Estimators,” and a FFT method for the Gaussian kernel described in Silverman, B. W. "Algorithm AS 176: Kernel density estimation using the fast Fourier transform." FFTs are computed using the [Eigen3 linear algebra library](https://eigen.tuxfamily.org). 



Here is an example usage. 

```bash
$ make && ./non_parametric 10000 .01 tests/test_data.bin output.bin
g++ -I /usr/local/include/eigen3 -std=c++11 -shared -fPIC -o libfft.so libfft.cpp
gcc -O3 -ftree-vectorize -ffast-math -msse2 -funroll-loops -march=native -mfpmath=sse -fstrict-aliasing -std=c99 -L. -lfft -lm  -o non_parametric non_parametric.c kernels.c binning.c utilities.c
Expecting 10000 packed doubles in tests/test_data.bin and writing density estimation with bandwidth 0.010000 to output.bin
wrote 16384 points of the estimate density function to output.bin
KDE Domain Parameters
import numpy
grid, delta = numpy.linspace(-2.795456, 2.978763, 16384, retstep=True)
```

Now check our estimate by integration.

```python
>>> import numpy as np
>>> from numpy import linspace, trapz, fromfile
>>> grid, delta = linspace(-2.795456, 2.978763, 16384, retstep=True)
>>> kde = fromfile("output.bin", dtype=np.double)
>>> trapz(kde, grid)
0.9999999423722401
```
