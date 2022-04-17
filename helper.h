#include <stdlib.h>
#ifndef HELPER_H
#define HELPER_H
	
float randomAlphaf(); 
double randomAlphad();


float* randomVectorf(const int N);
double* randomVectord(const int N);
float* randomMatrixf(const int M, const int N);
double* randomMatrixd(const int M, const int N);


#endif
