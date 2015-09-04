
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include "fftw3.h"
//#include "PsiDOf.h"


//void objective_real(complex<float> const *,float const *,int ,int ,int,double , int  , double , double , int , double , double , float , float *);

void objective_real_p(fftwf_complex const *Np_b, float const *coeff_c, int nxc, int nzc, int k, int p,float x0, int nx, float dx, float z0, int nz, float dz, float m_ord, float *Obj, float *Grad,int flag_grad,float lam);

void objective_real(fftwf_complex const *U,fftwf_complex const *QUtest, float const *coeff_c, int nxc, int nzc, int k, float x0, int nx, float dx, float z0, int nz, float dz, float m_ord, float *Obj, float *Grad,int flag_grad,float lam);

void make_inv_real(fftwf_complex const *U, fftwf_complex * inv, float const *coeff_c, int nxc, int nzc, int k, float x0, int nx, float dx, float z0, int nz, float dz, float m_ord);

void make_inv_real_p(fftwf_complex const *Np_b, fftwf_complex * inv, float const *coeff_c, int nxc, int nzc, int k,int p, float x0, int nx, float dx, float z0, int nz, float dz, float m_ord);
