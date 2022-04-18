#include <stdlib.h>
#include <stdio.h>
#include "benchmark_stencil.h"
#include "../lib/helper.h"
void bench_stencil(int k, int resolution, int printStyle) {
	float *X, *Y;
	long operations = 0;
	if (!resolution) {
		X = (float*) malloc(1920 * 1080 * sizeof(float));
		Y = (float*) malloc(1920 * 1080 * sizeof(float));
	 	operations = 1920ll * 1080;	
	}
	else {
		X = (float*) malloc(3840 * 2160 * sizeof(float));
		Y = (float*) malloc(3840 * 2160 * sizeof(float));
	 	operations = 3840ll * 2160;	
	}

	int count = k * k + operations;
	operations *= 1ll * k * k;
	float *S = randomMatrixf(k, k);
	CALCTIME(stencil, X, resolution, k, S, Y);
	if (printStyle == 1)
		printf("Resolution: %d x %d", (resolution) ? 3840 : 1920, (resolution) ? 2160 : 1080);
	printBenchmark(k, operations, count, calcTime, 4, printStyle);
}
