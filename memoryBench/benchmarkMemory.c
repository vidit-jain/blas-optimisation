#include <stdio.h>
#include "bench.h"
int main() {
	long vecSize = (1 << 28);
	bench_memory(vecSize);
}
