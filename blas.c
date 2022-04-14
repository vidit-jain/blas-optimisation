#include "benchmark_blas.h"
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
int main(int argc, char* argv[]) {
	if (argc != 5) {
		printf("Incorrect number of arguments passed\n");
		return 0;
	}
	char* function = argv[1];
	long startingN = atoi(argv[2]);
	int factor = atoi(argv[3]);
	int iterations = atoi(argv[4]);
	srand((unsigned)time(NULL));

	printf("Benchmarking %s\n", function);
	if (!strcmp("sscal", function)) {
		long N = startingN;
		for (int i = 0; i < iterations; i++, N *= factor) {
			bench_cblas_sscal(N);
		}	
	}
}
