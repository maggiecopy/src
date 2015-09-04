#include "objective_real.h"
#include <time.h>
int main(int argc, char **argv) {
  float  x0=-5;
  int  nx=101;
  float  dx=1;

  float  z0=-3;
  int  nz=101;
  float  dz=1;

  int  nxc=11;
  int nzc=11;
  float m_ord=-2;
  int k=3;
  float lam=0.01;

  float * Obj=(float *)malloc(sizeof(float));
  float * Obj2=(float *)malloc(sizeof(float));
  float *Grad=(float *)malloc(nxc*nzc*(2*k-1)*sizeof(float));
  float complex  *U = (float complex *)calloc(nx*nz,sizeof(float complex));
  float complex  *QUtest = (float complex *)calloc(nx*nz,sizeof(float complex));
  float x,z;
  srand ( time(NULL) );
  int m,n,l;
  for (m=0;m<nx;m++)
    {
      for (n=0;n<nz;n++)
	{
	  x=x0+m*dx;
	  z=z0+n*dz;
	 
	  U[m*nz+n]=(-4+4*x*x+4*z*z)*exp(-x*x-z*z)+ 0*I;
	  QUtest[m*nz+n]=exp(-x*x-z*z)+0*I;
	
	}
    }

  float * coeff_c = (float *)malloc(nxc*nzc*(2*k-1)*sizeof(float));

 

  for(m=0;m<nxc*nzc*(2*k-1);m++)
    {
      //coeff_c[m]=(float)rand()/RAND_MAX;//printf("%f \n",coeff_c[m]);
      coeff_c[m]=(float)sin(m);
    }
  
  objective_real(U,QUtest, coeff_c, nxc, nzc,  k, x0, nx,  dx, z0, nz, dz,  m_ord, Obj,Grad,1,lam);
  float * diff= (float *)malloc(nxc*nzc*(2*k-1)*sizeof(float));
  float h = 1e-3;
  printf("Objective = %f \n",Obj[0]);

  for (m=0;m<nxc*nzc*(2*k-1);m++)
    {

   coeff_c[m]=coeff_c[m]+h;
  objective_real(U,QUtest, coeff_c, nxc, nzc,  k, x0, nx,  dx, z0, nz, dz,  m_ord, Obj2,Grad,0,lam);
  diff[m]=(Obj2[0]-Obj[0])/h;
  coeff_c[m] = coeff_c[m]-h;
    }

  float error = 0.0;
  float norm_grad = 0.0;
  float min_grad = Grad[0]*Grad[0];
  int index=0;
for (m=0;m<nxc*nzc*(2*k-1);m++)
    {
      //printf("Grad %f \n",Grad[m]);
      error = error + (diff[m]-Grad[m])*(diff[m]-Grad[m]);
      norm_grad=norm_grad+Grad[m]*Grad[m];
      if(Grad[m]*Grad[m] < min_grad)
	{
	  min_grad = Grad[m]*Grad[m];
	  index=m;
	}
      
    }
 printf("Relative error between finite difference and gradient= %f \n",pow(error/norm_grad,0.5));
 printf("Norm Grad = %f \n",pow(norm_grad,0.5));
 printf("Min Grad  = %f \n", pow(min_grad,0.5));
 printf("At index= %d \n",index);
 printf("Finite difference at min = %f \n",diff[index]);

  //printf("finite diff  = %6.6f \n", (Obj2[0]-Obj[0])/1e-3);
  //printf("Grad         = %6.6f \n",Grad[242]);
  /*
  printf("\n");
  for ( m=1;m<= nxc*nzc*(2*k-1);m++)
    {
      printf("%f",Grad[m-1]);printf("   ");
      if(( m % nxc) ==0) printf("\n");
      if((m % (nxc*nzc)) == 0){ printf("\n");printf("\n");printf("\n");
      }
    } printf("\n");
  */
 free(Obj);free(Obj2);free(Grad);free(diff);free(U);free(QUtest);
}
