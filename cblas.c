#pragma GCC target("avx2")
#pragma GCC target("sse")
#pragma GCC optimize("O3")
#include "cblas.h"
#include <math.h>
#include <x86intrin.h>
#include <immintrin.h>
#include <stdio.h>
// BLAS Level 1
void cblas_sscal(const int N, const float alpha, float *X, const int incX) {
	if (incX == 1) {
		__m256 _alpha = _mm256_set1_ps(alpha);
		for (int i = 0; i < N; i+=8) {
			__m256 x = _mm256_loadu_ps(&X[i]);
			__m256 z = _mm256_mul_ps(x, _alpha);
			_mm256_storeu_ps(&X[i], z);
		}
	}
	else {
		for (int i = 0; i < N; i++) {
			X[i * incX] *= alpha;
		}
	}
}

void cblas_dscal(const int N, const double alpha, double *X, const int incX) {
	if (incX != 1) {
		__m256d _alpha = _mm256_set1_pd(alpha);
		for (int i = 0; i < N; i+=4) {
			__m256d x = _mm256_loadu_pd(&X[i]);
			__m256d z = _mm256_mul_pd(x, _alpha);
			_mm256_storeu_pd(&X[i], z);
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
	if (incX != 1 && incY == 1) {
		__m256 ans = _mm256_setzero_ps();
		for (int i = 0; i < N; i += 8) {
			__m256 x = _mm256_loadu_ps(&X[i]);
			__m256 y = _mm256_loadu_ps(&Y[i]);
			__m256 z = _mm256_mul_ps(x, y);
			__m256 a = _mm256_add_ps(ans, z);
			ans = a;
		}
		float* ansArr = (float*) &ans;
		for (int i = 0; i < 8; i++) value += ansArr[i];
	}
	else {
#pragma omp parallel for reduction(+:value)
		for (int i = 0; i < N; i++) {
			value += X[i] * Y[i];
		}
	}
	return value;
}

void cblas_saxpy(const int N, const float alpha, const float *X,
                 const int incX, float *Y, const int incY) {
	if (incY == 1 && incX == 1) {
		__m256 _alpha = _mm256_set1_ps(alpha);
		for (int i = 0; i < N; i += 8) {
			__m256 x = _mm256_loadu_ps(&X[i]);
			__m256 y = _mm256_loadu_ps(&Y[i]);
			__m256 z = _mm256_fmadd_ps(x, _alpha, y);
			_mm256_storeu_ps(&Y[i], z);
		}
	} else {
		for (int i = 0; i < N; i++) {
			Y[i * incY] = X[i * incX] * alpha + Y[i *incY];
		}
	}
}

// BLAS Level 2

void cblas_sgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY) {
	int n = N;
	int m = M;
	if (TransA != CblasNoTrans) {
		n = M;
		m = N;
	}
	cblas_sscal(m, beta, Y, incY);
	if ((order == CblasRowMajor && TransA == CblasNoTrans) || 
			(order == CblasRowMajor && TransA != CblasNoTrans)) {
#pragma omp parallel for
		for (int i = 0; i < m; i++) {
			float value = 0.0;
			for (int j = 0; j < n; j++) {
				value += X[j * incX] * A[(lda * i + j)];
			}
			Y[i * incY] += alpha * value;
		}
	}
	else {
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			float value = alpha * X[i * incX];
			for (int j = 0; j < m; j++) {
				Y[j * incY] += value * A[(lda * i + j)];
			}
		}
	}
}
// Implementation that uses sdot, sscal
/*
void cblas_sgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY) {
	if (order == CblasRowMajor && TransA == CblasNoTrans) {
		cblas_sscal(M, alpha, Y, incY);
#pragma omp parallel for
		for (int i = 0; i < M; i++) {
			float value = cblas_sdot(N, X, incX, A, 1);
			Y[i * incY] += alpha * value;
		}
	}
}*/

