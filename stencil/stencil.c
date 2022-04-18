#include "stencil.h"
void stencil(float *X, const enum ImageType typ, int k, float *S, float *Y) {
    int N = 1920, M = 1080;

    if (typ == UHD) {
        N = 3840;
        M = 2160;
    }

#pragma omp parallel for
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            Y[i * N + j] = 0.0;
            for (int x = 0; x < k; x++) {
                if (i + x >= M) break;
                for (int y = 0; y < k; y++) {
                    if (j + y >= N) break;
                    Y[i * N + j] += X[(i + x) * N + (j + y)] * S[x * k + y];
                }
            }
        }
    }
}
