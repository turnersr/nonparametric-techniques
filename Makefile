C_FILES = non_parametric.c kernels.c binning.c utilities.c
CPP_FILES = libfft.cpp
OUT_FFT_LIB = libfft.so
OUT_EXE = non_parametric
EIGEN_PATH = /usr/local/include/eigen3

CPP = g++ -I $(EIGEN_PATH) -std=c++11 -shared -fPIC
CC = gcc -O3 -ftree-vectorize -ffast-math -msse2 -funroll-loops -march=native -mfpmath=sse -fstrict-aliasing -std=c99 -L. -lfft -lm 



build: $(FILES)
	$(CPP) -o $(OUT_FFT_LIB) $(CPP_FILES)
	$(CC) -o $(OUT_EXE) $(C_FILES) 

clean:
	rm -f *.so *.o *.asm non_parametric

rebuild:
	clean build

asm_dump:
	objdump -d non_param_kernel > non_param_kernel.asm
