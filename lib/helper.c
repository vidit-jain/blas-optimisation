#include <stdlib.h>
#include "helper.h"
float randomAlphaf() { return (float) rand() / RAND_MAX; }
double randomAlphad() { return (double) rand() / RAND_MAX; }


float* randomVectorf(const int N) {
	float *X = (float*) malloc(N * sizeof(float));
	for (int i = 0; i < N; i++) {
		X[i] = (float) rand() / RAND_MAX;
	}
	return X;
}

double* randomVectord(const int N) {
	double *X = (double*) malloc(N * sizeof(double));
	for (int i = 0; i < N; i++) {
		X[i] = (double) rand() / RAND_MAX;
	}
	return X;
}

float* randomMatrixf(const int M, const int N) {
	float *A = (float*) malloc(M * N * sizeof(float));
	for (int i = 0; i < M * N; i++) {
		A[i] = (float) rand() / RAND_MAX;
	}
	return A;
}

double* randomMatrixd(const int M, const int N) {
	double *A = (double*) malloc(M * N * sizeof(double));
	for (int i = 0; i < M * N; i++) {
		A[i] = (double) rand() / RAND_MAX;
	}
	return A;
}
