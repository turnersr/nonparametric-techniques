CCv5 = gcc -O3 -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5 -funroll-loops -march=native -mfpmath=sse -fstrict-aliasing -std=c99 
CCv0 = gcc -O3 -ftree-vectorize -msse2 -funroll-loops -march=native -mfpmath=sse -fstrict-aliasing -std=c99 
CCv1 = gcc -O3 -ftree-vectorize -msse2 -march=native -mfpmath=sse -fstrict-aliasing -std=c99 

CCv2 = gcc -O3 -ffast-math -std=c99
CC = gcc -O3 -ftree-vectorize -ffast-math -msse2 -funroll-loops -march=native -mfpmath=sse -fstrict-aliasing -std=c99 -L. -lfft

#gcc main.c -L. -lfoo -o main


FILES = non_param_kernel.c
OUT_EXE = non_param_kernel

build: $(FILES)
	g++ -I /usr/local/Cellar/eigen/3.3.4/include/eigen3/ -std=c++11 fft.cpp -shared -fPIC -o libfft.so
	$(CC) -o $(OUT_EXE) $(FILES) -lm

clean:
	rm -f *.o *.asm *.bin non_param_kernel

rebuild:
	clean build

asm_dump:
	objdump -d non_param_kernel > non_param_kernel.asm
