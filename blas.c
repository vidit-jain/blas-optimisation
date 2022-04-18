#include "benchmark_blas.h"
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
int funcCount = 10;
char* functionNames[] = {"sscal", "dscal", "sdot", "ddot", "saxpy", "daxpy",
						"sgemv", "dgemv", "sgemm", "dgemm"};
void (*fp[])(int, int) = {&bench_cblas_sscal, &bench_cblas_dscal, &bench_cblas_sdot, 
						&bench_cblas_ddot, &bench_cblas_saxpy, &bench_cblas_daxpy, 
						&bench_cblas_sgemv, &bench_cblas_dgemv, &bench_cblas_sgemm,
						&bench_cblas_dgemm}; 

int main(int argc, char* argv[]) {
	if (argc != 6) {
		printf("Incorrect number of arguments passed\n");
		return 0;
	}
	char* function = argv[1];
	int style = atoi(argv[2]);
	long startingN = atoi(argv[3]);
	int factor = atoi(argv[4]);
	int iterations = atoi(argv[5]);
	srand((unsigned)time(NULL));

	if (style == 1)
		printf("Benchmarking %s\n", function);

	void (*f)(int, int) = NULL;
	for (int i = 0; i < funcCount; i++) {
		if (!strcmp(functionNames[i], function)) f = fp[i]; 
	}

	long N = startingN;
	for (int i = 0; i < iterations; i++, N *= factor) {
		if (f)
			f(N, style);
	}	
}
