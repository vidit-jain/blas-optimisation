#include "bench.h"
#include <stdlib.h>
#include "../lib/benchmark.h"
#include <stdio.h>

void bench_memory(long vecSize) {
	long* a = (long*) malloc(vecSize * sizeof(long));
	long* b = (long*) malloc(vecSize * sizeof(long));
	struct timeval t;
	for (int i = 0; i < vecSize; i++) {
		if (i % 2)  {
			a[i] = -1;
			b[i] = 2;
		}
		else {
			a[i] = 2;
			b[i] = -1;
		}
	}
	long c = 0;

	tick(&t);
	for (int i = 0; i < vecSize; i++) {
		c = a[i] + b[i];
		
	}
	printf("%ld", c);
	double calcTime = tock(&t);
	long count = 2ll * vecSize;
	printBenchmark(vecSize, 1, count, calcTime, 8, 1);
}
