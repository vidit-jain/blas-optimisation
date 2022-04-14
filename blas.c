#include <stdlib.h>
#include <stdio.h>
#include "benchmark_blas.h"
int main() {
	bench_cblas_sscal(10000000);
}
