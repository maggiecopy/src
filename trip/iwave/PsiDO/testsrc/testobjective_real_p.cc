#include "objective_real.h"
#include <time.h>
#include<math.h>

int main(int argc, char **argv) {
  float  x0=-5; float xmax=5;
  int  nx=101;
  float  dx=(xmax-x0)/(nx-1);

  float  z0=-3;float zmax=4;
  int  nz=101;
  float  dz=(zmax-z0)/(nz-1);

  int  nxc=5;
  int nzc=5;
  float m_ord=-2;
  int k=3;
  int p=2;
  float lam=0.0;

  int loop_max=100;

  float * Obj=(float *)malloc(sizeof(float));
  float * Obj2=(float *)malloc(sizeof(float));
  float *Grad=(float *)calloc(nxc*nzc*(2*k-1)*p,sizeof(float));
  float *Grad2=(float *)calloc(nxc*nzc*(2*k-1)*p,sizeof(float));
  float *coeff_pert=(float *)calloc(nxc*nzc*(2*k-1)*p,sizeof(float));
  float *coeff_c2=(float *)calloc(nxc*nzc*(2*k-1)*p,sizeof(float));
  float complex  *Np_b = (float complex *)calloc(nx*nz*(p+1)*p,sizeof(float complex));
  float complex  *U11 = (float complex *)calloc(nx*nz,sizeof(float complex));
  float complex  *U12 = (float complex *)calloc(nx*nz,sizeof(float complex));
  float complex  *U13 = (float complex *)calloc(nx*nz,sizeof(float complex));
  float complex  *U21= (float complex *)calloc(nx*nz,sizeof(float complex));
  float complex  *U22 = (float complex *)calloc(nx*nz,sizeof(float complex));
  float complex  *U23 = (float complex *)calloc(nx*nz,sizeof(float complex));


  // printf("memory allocated for vectors \n");
  float x,z;
  srand ( time(NULL) );
  int m,n;
  
  for (m=0;m<nx;m++)
    {
      for (n=0;n<nz;n++)
	{
	  x=x0+m*dx;
	  z=z0+n*dz;
	 
	  U11[m*nz+n]=(-4+4*x*x+4*z*z)*exp(-x*x-z*z)+ 0*I;
	  U12[m*nz+n]=z*exp(-x*x-z*z)+0*I;
	  U13[m*nz+n]=x*exp(-x*x-z*z)+0*I;
	  U21[m*nz+n]=(x+z)*exp(-x*x-z*z)+0*I;
	  U22[m*nz+n]=x*x*exp(-x*x-z*z)+0*I;
	  U23[m*nz+n]=z*z*exp(-x*x-z*z)+0*I;
	  //printf("%d   ",m*nz+n);
	}
    }
  
  //printf("starting memcopy in test \n");
  
  memcpy(&Np_b[0],U11,nx*nz*sizeof(float complex));

  memcpy(&Np_b[nx*nz],U12,nx*nz*sizeof(float complex));
  
  if(p>1){
  memcpy(&Np_b[2*nx*nz],U13,nx*nz*sizeof(float complex));
  memcpy(&Np_b[3*nx*nz],U21,nx*nz*sizeof(float complex));
  memcpy(&Np_b[4*nx*nz],U22,nx*nz*sizeof(float complex));
  memcpy(&Np_b[5*nx*nz],U23,nx*nz*sizeof(float complex));
  }
 
 
  float * coeff_c = (float *)malloc(nxc*nzc*(2*k-1)*p*sizeof(float));

 

 
  for(m=0;m<nxc*nzc*(2*k-1)*p;m++)
    {
      coeff_c[m]=(float)rand()/RAND_MAX;//printf("%f \n",coeff_c[m]);
      coeff_pert[m]=(float)rand()/RAND_MAX;
      // coeff_c[m]=(float)sin(m);
    }

   /////Timing test

 clock_t start,end;
 double elapsed;
 start = clock();
 objective_real_p(Np_b,coeff_c, nxc, nzc,  k,p, x0, nx,  dx, z0, nz, dz,  m_ord, Obj,Grad,1,lam);
 end=clock();
 elapsed = ((double)(end-start))/CLOCKS_PER_SEC;
 printf(" time taken: %lf \n",elapsed);
 printf("Objective = %f \n",Obj[0]);

float * diff= (float *)malloc(nxc*nzc*(2*k-1)*p*sizeof(float));
 
 
  //printf("safely here \n");
  float h =1e-2;
 float norm_pert=0.0;
  for (m=0;m<nxc*nzc*(2*k-1)*p;m++)
    {
      //printf("%f \n",coeff_c[m]);printf("%f \n",coeff_pert[m]);
   coeff_c[m]=coeff_c[m]+h;
  
if(nxc*nzc<loop_max)
   {
  objective_real_p(Np_b, coeff_c, nxc, nzc,  k,p, x0, nx,  dx, z0, nz, dz,  m_ord, Obj2,Grad,0,lam);
  diff[m]=(Obj2[0]-Obj[0])/h;
  printf("%d",m); printf("  Diff %f",diff[m]);printf("  Grad %f",Grad[m]);printf("  Error %f\n",Grad[m]-diff[m]);
   }
  coeff_c[m] = coeff_c[m]-h;
   norm_pert=norm_pert+coeff_pert[m]*coeff_pert[m];
    }
  norm_pert=pow(norm_pert,0.5);

  float error = 0.0;
  float norm_grad = 0.0;
 
  float min_grad = Grad[0]*Grad[0];
  int index=0;
  float dot=0;
for (m=0;m<nxc*nzc*(2*k-1)*p;m++)
    {
     
      error = error + (diff[m]-Grad[m])*(diff[m]-Grad[m]);
      norm_grad=norm_grad+Grad[m]*Grad[m];
      
      coeff_pert[m]=coeff_pert[m]/norm_pert;
      coeff_c2[m]=coeff_c[m]+h*coeff_pert[m];
      dot=dot+Grad[m]*coeff_pert[m];
      if(Grad[m]*Grad[m] < min_grad)
	{
	  min_grad = Grad[m]*Grad[m];
	  index=m;
	}
      
    }
objective_real_p(Np_b, coeff_c2, nxc, nzc,  k,p, x0, nx,  dx, z0, nz, dz,  m_ord, Obj2,Grad2,1,lam);


//  FILE * f_grad;
//  f_grad=fopen("./grad.H@","w");
//  fwrite(Grad, sizeof(float), nxc*nzc*(2*k-1)*p, f_grad);
//  fclose(f_grad);

//  FILE * f_diff;
//  f_diff=fopen("./diff.H@","w");
//  fwrite(diff, sizeof(float), nxc*nzc*(2*k-1)*p, f_diff);
//  fclose(f_diff);
//   

 float grad_diff=0;
for (m=0;m<nxc*nzc*(2*k-1)*p;m++)
    {
      grad_diff=grad_diff+(Grad[m]-Grad2[m])*(Grad[m]-Grad2[m]);
    }
 grad_diff=pow(grad_diff,0.5);

if(nxc*nzc<loop_max)
   {
 printf("Relative error between finite difference and gradient= %f \n",pow(error/norm_grad,0.5));
 printf("Norm Grad = %f \n",pow(norm_grad,0.5));
 printf("Min Grad  = %f \n", pow(min_grad,0.5));
 printf("At index= %d \n",index);
 printf("Finite difference at min = %f \n",diff[index]);
   }
 printf("f(x+h) =%f \n", Obj2[0]);
 printf("f(x)   =%f \n",Obj[0]);
 printf("f(x)+(h,Grad) = %f \n",Obj[0]+h*dot);
 printf("Relative error = %f \n",fabs(((Obj2[0]-Obj[0])/h-dot)/Obj[0]));//this error is too optimistic!
 printf("Relative error2 = %f \n",fabs(((Obj2[0]-Obj[0])/h-dot)/(grad_diff/h)));//the denominator is my cheap attempt to approximate the norm of the hessian.
 printf("Relative error3 = %f \n",fabs(((Obj2[0]-Obj[0])/h-dot)/((Obj2[0]-Obj[0])/h)));

 free(Np_b);free(U11);free(U12);free(Grad);free(diff);free(coeff_pert);free(coeff_c2);free(U13);free(U21);free(U22);free(U23);
}
