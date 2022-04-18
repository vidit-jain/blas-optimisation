#include <stdlib.h>
#include <stdio.h>
#include "benchmark.h"
void tick(struct timeval *t) {
    gettimeofday(t, NULL);
}

double tock(struct timeval *t) {
    struct timeval now;
    gettimeofday(&now, NULL);
    return(double) (now.tv_sec - t->tv_sec) + 
    ((double)(now.tv_usec - t->tv_usec)/1000000.);
}


void printBenchmark(long N, long operations, long count, float calcTime, int size, int printStyle) {
	double gflops = 1ll * operations / calcTime * 1e-9;
	double memoryBandwidth = 1ll * size * count * 1e-9 / calcTime; 
	double timeMS = calcTime * 1000;
	if (printStyle == 1) {
		printf("\nValue of N:\t\t\t%ld\n", N);
		printf("Time (ms):\t\t\t%-3.3lf\n", timeMS);
		printf("Memory bandwidth (GB/s):\t%-3.3lf\n", memoryBandwidth);
		printf("Computing Throughput (GFLOPS):\t%-3.3lf\n\n", gflops);
	}
	else {
		printf("%ld,%f,%f,", N, gflops, memoryBandwidth);
	}
}
