#include <stdlib.h>
#include <stdio.h>
#include "benchmark_stencil.h"
int main(int argc, char* argv[]) {
	if (argc != 4) {
		printf("Incorrect number of arguments passed\n");
		return 0;
	}
	int style = atoi(argv[1]);
	int resolution = atoi(argv[2]);
	int k = atoi(argv[3]);
	bench_stencil(k, resolution, style);	
}
