#include <stdio.h>
#include "bench.h"
#include "../lib/benchmark.h"
int main() {
	const int iterations = 6;
	benchmark_multi(iterations, 1000000000);
	benchmark_single(1000000000);
}
