#pragma GCC target("avx2")
#pragma GCC target("sse")
#pragma GCC optimize("O3")
#include <x86intrin.h>
#include <immintrin.h>
#include "bench.h"
#include "../lib/helper.h"
#include "../lib/benchmark.h"
#include <stdio.h>
double f2(int loops, int vecSize) {
	float alpha = randomAlphaf();
	float* X = randomVectorf(vecSize);
	float* Y = randomVectorf(vecSize);
	__m256 _alpha = _mm256_set1_ps(alpha);
	struct timeval t;
	tick(&t);

	for (int i = 0; i < loops; i++) {
		for (int j = 0; j < vecSize; j += 8) {
			__m256 x = _mm256_loadu_ps(&X[j]);
			__m256 y = _mm256_loadu_ps(&Y[j]);
			__m256 z = _mm256_fmadd_ps(x, _alpha, y);
			_mm256_storeu_ps(&Y[j], z);
		}

	}
	return tock(&t);
}
void benchmark_multi(long threads, long loops) {
	float totalTime = 0;
	const int vecSize = 64;
#pragma omp parallel for 
	for (int i = 0; i < threads; i++) {
		totalTime = f2(loops, vecSize);	
	}

	long operations = 2ll * threads * loops * vecSize;
	long count = 1000; 
	printBenchmark(vecSize, operations, count, totalTime, 4, 1);
}

void benchmark_single(long loops) {
	float totalTime = 0;
	const int vecSize = 128;
	totalTime = f2(loops, vecSize);	
	long operations = 2ll * loops * vecSize;
	long count = 1000; 
	printBenchmark(vecSize, operations, count, totalTime, 4, 1);
}

