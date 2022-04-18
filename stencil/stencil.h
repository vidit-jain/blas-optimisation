#ifndef STENCIL_H
#define STENCIL_H
enum ImageType {HD=0, UHD=1};
void stencil(float *X, const enum ImageType typ, int k, float *S, float *Y);
#endif
