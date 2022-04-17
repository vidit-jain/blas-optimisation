#include "cblas.h"
#include <math.h>
// BLAS Level 1
void cblas_sscal(const int N, const float alpha, float *X, const int incX) {
	if (incX == 1) {
#pragma omp simd 
		for (int i = 0; i < N; i++) {
			X[i] *= alpha;
		}
	}
	else {
		for (int i = 0; i < N; i++) {
			X[i * incX] *= alpha;
		}
	}
}

float cblas_sdot(const int N, const float  *X, const int incX,
                  const float  *Y, const int incY) {
	float value = 0.0;	
//#pragma omp parallel for 
	for (int i = 0; i < N; i++) {
		value += X[i * incX] * Y[i * incY];
	}
	return value;
}

void cblas_saxpy(const int N, const float alpha, const float *X,
                 const int incX, float *Y, const int incY) {
	for (int i = 0; i < N; i++) {
		Y[i * incY] = X[i * incX] * alpha + Y[i *incY];
	}
}

// BLAS Level 2

void cblas_sgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY) {
	if (order == CblasRowMajor && TransA == CblasNoTrans) {
		for (int i = 0; i < M; i++) {
			Y[i * incY] *= beta;
			float value = 0.0;
			for (int j = 0; j < N; j++) {
				value += X[j * incX] * A[(lda * i + j)];
			}
			Y[i * incY] += alpha * value;
		}
	}
	else if (order == CblasRowMajor && TransA == CblasNoTrans){
		for (int i = 0; i < M; i++) Y[i *incY] *= beta;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < M; j++) {
				Y[i * incY] += alpha * X[i * incY] * A[(lda * i + j)];
			}
		}
	}
	else if (order == CblasRowMajor && (TransA == CblasTrans || TransA == CblasConjTrans)) {
		for (int i = 0; i < N; i++) {
			Y[i * incY] *= beta;
		}
	}
}
