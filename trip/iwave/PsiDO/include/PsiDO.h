#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <complex.h>
#include "fftw3.h"


void PsiDO(complex const *U ,complex const *FQ,complex *QU ,int nx, int nz, int k,
	   double dx,double dz,double m_ord);

void PsiDOf(float complex const *U ,float complex const *FQ,float complex *QU ,int nx, int nz, int k,
	   double dx,double dz,double m_ord);


void PsiDO_real(complex const  *U ,complex const *FQ,complex *QU ,int nx, int nz, int k, double dx,double dz,double m_ord);

void PsiDOf_real(float complex const  *U ,float complex const *FQ,float complex *QU ,int nx, int nz, int k, double dx,double dz,double m_ord);

void PsiDO_l(complex const  *U ,complex const *FQ,complex *QU ,int nx, int nz,double dx,double dz,double m_ord,int l);

void PsiDOf_l(float complex const  *U ,float complex const *FQ,float complex *QU ,int nx, int nz,double dx,double dz,double m_ord,int l);
