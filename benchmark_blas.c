#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include "helper.h"
#include "benchmark_blas.h"
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

void bench_cblas_sscal(const int N, int printStyle) {
	const float alpha = randomAlphaf(); 
	float *X = randomVectorf(N) ;
	CALCTIME(cblas_sscal, N, alpha, X, 1);
	printBenchmark(N, N, N, calcTime, 4, printStyle);
}

void bench_cblas_dscal(const int N, int printStyle) {
	const double alpha = randomAlphad(); 
	double *X = randomVectord(N) ;
	CALCTIME(cblas_dscal, N, alpha, X, 1);
	printBenchmark(N, N, N, calcTime, 8, printStyle);
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

void bench_cblas_sgemv(const int N, int printStyle) {
	const float alpha = randomAlphaf();
	const float beta = randomAlphaf();
	float *X = randomVectorf(N);
	float *Y = randomVectorf(N);
	float *A = randomMatrixf(N, N);
	CALCTIME(cblas_sgemv, CblasRowMajor, CblasNoTrans, N, N, alpha, A, 1, X, 1, beta, Y, 1);
	int operations = N * (3 + 2 * N);
	int count = N * N + N + N;
	printBenchmark(N, operations, count, calcTime, 4, printStyle);
}

void bench_cblas_sgemm(const int N, int printStyle) {
	const float alpha = randomAlphaf();	
	const float beta = randomAlphaf();	
	float *A = randomMatrixf(N, N);
	float *B = randomMatrixf(N, N);
	float *C = randomMatrixf(N, N);
	CALCTIME(cblas_sgemm, CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, 
			N, alpha, A, N, B, N, beta, C, N);
	long operations = 1ll * 2 * N * N * N + 3 * N * N;
	long count = 3ll * N * N;
	printBenchmark(N, operations, count, calcTime, 4, printStyle);

}
