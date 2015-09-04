#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <complex.h>
#include <fftw3.h>

void conv(complex const *x,complex const *y,complex *c,int nx,int ny);
void Qsquare(complex const *FQ,complex *FQ2,int nx,int nz,int k);
void Qsym(complex const *FQ,complex *FQ2,int nx,int nz,int k);
void Qshift(complex const *FQ, complex *FQ2, int nx, int nz, int k, int shift);

void Qsquaref(float complex const *FQ,float complex *FQ2,int nx,int nz,int k);
void Qsymf(float complex const *FQ,float complex *FQ2,int nx,int nz,int k);
void Qshiftf(float complex const *FQ, float complex *FQ2, int nx, int nz, int k, int shift);
