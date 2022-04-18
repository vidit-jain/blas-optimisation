#ifndef BENCH_H
#define BENCH_H
#include <sys/time.h>
#define CALCTIME(functionName, args...) \
	struct timeval calc;\
	tick(&calc);\
	functionName(args);\
	double calcTime = tock(&calc);\

void tick(struct timeval *t);
double tock(struct timeval *t);
void printBenchmark(long N, long operations, long count, float calcTime, int size, int printStyle);
#endif
