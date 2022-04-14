#include <stdlib.h>
#include "cblas.h"
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include "helper.h"
#define CALCTIME(functionName, args...) \
	struct timeval calc;\
	tick(&calc);\
	functionName(args);\
	double calcTime = tock(&calc);\

	
void tick(struct timeval *t) {
    gettimeofday(t, NULL);
}

double tock(struct timeval *t) {
    struct timeval now;
    gettimeofday(&now, NULL);
    return(double) (now.tv_sec - t->tv_sec) + 
    ((double)(now.tv_usec - t->tv_usec)/1000000.);
}


void printBenchmark(int N, int operations, int count, float calcTime, int size) {
	float gflops = operations / calcTime * 1e-9;
	float memoryBandwidth = size * count * 1e-9 / calcTime; 
	float timeMS = calcTime * 1000;
	printf("\nValue of N:\t\t\t%d\n", N);
	printf("Time (ms):\t\t\t%-3.3f\n", timeMS);
	printf("Memory bandwidth (GB/s):\t%-3.3f\n", memoryBandwidth);
	printf("Computing Throughput (GFLOPS):\t%-3.3f\n\n", gflops);
}

void bench_cblas_sscal(const int N) {
	const float alpha = randomAlphaf(); 
	float *X = randomVectorf(N) ;
	CALCTIME(cblas_sscal, N, alpha, X, 1);
	int operations = N;
	int count = N;
	int size = 4;
	printBenchmark(N, operations, count, calcTime, size);
}
