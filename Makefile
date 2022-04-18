BLAS_SRC = blas
STENCIL_SRC = stencil
LIB = lib
CC = gcc
CFLAGS = -O3 -fopenmp -mfma -fopenmp -fopenmp-simd -mavx2
all:${BLAS_SRC}/*.c ${LIB}/*.c
	$(CC) $(CFLAGS) ${LIB}/*.c $(BLAS_SRC)/*.c -o bin/cblas
	$(CC) $(CFLAGS) ${LIB}/*.c $(STENCIL_SRC)/*.c -o bin/stencil

clean: 
	rm *.o 

