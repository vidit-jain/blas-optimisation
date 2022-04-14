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
