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


void printBenchmark(int N, int operations, int count, float calcTime, int size, int printStyle) {
	float gflops = 1ll * operations / calcTime * 1e-9;
	float memoryBandwidth = 1ll * size * count * 1e-9 / calcTime; 
	float timeMS = calcTime * 1000;
	if (printStyle == 1) {
		printf("\nValue of N:\t\t\t%d\n", N);
		printf("Time (ms):\t\t\t%-3.3f\n", timeMS);
		printf("Memory bandwidth (GB/s):\t%-3.3f\n", memoryBandwidth);
		printf("Computing Throughput (GFLOPS):\t%-3.3f\n\n", gflops);
	}
	else {
		printf("%d,%f,%f,", N, gflops, memoryBandwidth);
	}
}

void printArrays(int N, int operations, int count, float calcTime, int size) {
	float gflops = operations / calcTime * 1e-9;
	float memoryBandwidth = size * count * 1e-9 / calcTime; 
	printf("%d,%f,%f,", N, gflops, memoryBandwidth);
}

void bench_cblas_sscal(const int N, int printStyle) {
	const float alpha = randomAlphaf(); 
	float *X = randomVectorf(N) ;
	CALCTIME(cblas_sscal, N, alpha, X, 1);
	printBenchmark(N, N, N, calcTime, 4, printStyle);
}

void bench_cblas_sdot(const int N, int printStyle) {
	float *X = randomVectorf(N);
	float *Y = randomVectorf(N);
	CALCTIME(cblas_sdot, N, X, 1, Y, 1);
	int operations = 2 * N;
	int count = 2 * N;
	printBenchmark(N, operations, count, calcTime, 4, printStyle);
}

void bench_cblas_saxpy(const int N, int printStyle) {
	const float alpha = randomAlphaf(); 
	float *X = randomVectorf(N) ;
	float *Y = randomVectorf(N) ;
	CALCTIME(cblas_saxpy, N, alpha, X, 1, Y, 1);
	int operations = 2 * N;
	int count = 2 * N;
	printBenchmark(N, operations, count, calcTime, 4, printStyle);
}
