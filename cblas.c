#include "cblas.h"
#include <math.h>

void cblas_sscal(const int N, const float alpha, float *X, const int incX) {
//#pragma omp parallel for 
	for (int i = 0; i < N; i++) {
		X[i * incX] *= alpha;
	}
}

float cblas_sdot(const int N, const float  *X, const int incX,
                  const float  *Y, const int incY) {
	float value = 0.0;	
	for (int i = 0; i < N; i++) {
		value += X[i * incX] + Y[i * incY];
	}
	return value;
}

void cblas_saxpy(const int N, const float alpha, const float *X,
                 const int incX, float *Y, const int incY) {
	for (int i = 0; i < N; i++) {
		Y[i * incY] = X[i * incX] * alpha + Y[i *incY];
	}
}

