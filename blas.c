#include "benchmark_blas.h"
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
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

	void (*f)(int, int);
	if (!strcmp("sscal", function)) f = &bench_cblas_sscal;
	else if (!strcmp("dscal", function)) f = &bench_cblas_dscal;
	else if (!strcmp("sdot", function)) f = &bench_cblas_sdot;
	else if (!strcmp("saxpy", function)) f = &bench_cblas_saxpy;
	else if (!strcmp("sgemv", function)) f = &bench_cblas_sgemv;
	else if (!strcmp("sgemm", function)) f = &bench_cblas_sgemm;
	else f = NULL;

	long N = startingN;
	for (int i = 0; i < iterations; i++, N *= factor) {
		if (f)
			f(N, style);
	}	
}
