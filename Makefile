BLAS_SRC = blas
STENCIL_SRC = stencil
MEMORY_SRC = memoryBench
CPU_SRC = cpuBench
LIB = lib
CC = gcc
CFLAGS = -O3 -mfma -fopenmp -fopenmp-simd -mavx2 -lm
all:${BLAS_SRC}/*.c ${LIB}/*.c ${STENCIL_SRC}/*.c ${MEMORY_SRC}/*.c ${CPU_SRC}/*.c
	$(CC) $(CFLAGS) ${LIB}/*.c $(BLAS_SRC)/*.c -o bin/cblas
	$(CC) $(CFLAGS) ${LIB}/*.c $(STENCIL_SRC)/*.c -o bin/stencil
	$(CC) $(CFLAGS) ${LIB}/*.c $(MEMORY_SRC)/*.c -o bin/memoryBench
	$(CC) $(CFLAGS) ${LIB}/*.c $(CPU_SRC)/*.c -o bin/cpuBench

clean: 
	rm *.o 

