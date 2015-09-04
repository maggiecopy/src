#include "PsiDO.h"
#include <time.h>

//tests the psido code using a special input vector U=(-4+4x^2+4z^2)*exp(-x^2 - z^2)
// with the symbol q = z*cos(x+theta)^ 2 , this combination of U and Q makes for a nice
//analytic formula for QU=2*z*(2*cos(x)^2*x^2-4*z*sin(x)*cos(x)*x+2*z^2-2*z%2*cos(x)^2-1)*exp(-x^2 - z^2)
//good other test are differentiation and integration operators.
int main(int argc, char **argv) {
  printf("starting...\n");

 int nx =1001;
 int nz=1001;
 int k=2;

 double xmin=-5.0;
 double xmax = 5.0;
 double zmin=-5.0;
 double zmax = 5.0;
 double x,z;

 double dx=(xmax-xmin)/((double)(nx-1));
 double dz = (zmax-zmin)/((double)(nz-1));
 double m_ord=0; /*for z*cos(x+theta)^2 */
 printf("declaring big vectors ...\n");
 
 float complex *U;
 float complex *QUtest;
 float complex *FQ;
 float complex *QU; 

  U = (float complex *) malloc(nx*nz*sizeof(float complex));printf("1 \n");
  FQ = (float complex *)  malloc(nx*nz*k*sizeof(float complex) );printf("2 \n");
  QUtest = (float complex *)  malloc(nx*nz*sizeof(float complex) );printf("3 \n");
  QU = (float complex *)  calloc(nx*nz,sizeof(float complex) );printf("4 \n");//initialize QU to zero using calloc
  
int m,n,l;

 printf("variables set...\n");

  /* U = (-4+4x^2+4z^2)*exp(-x^2 - z^2)  */ 
 
  for (m=0;m<nx;m++)
  {
    for (n=0;n<nz;n++)
    {
      x=xmin+m*dx;
      z=zmin+n*dz;
      //printf("%d",m*nz+n);
      U[m*nz+n]= (-4+4*x*x+4*z*z)*exp(-x*x-z*z);
      //U[m*nz+n]=exp(-x*x-z*z);
    }
  }
  

  



 printf("setting U ok...\n");
 // exit(0);
 /* Fourier Expansion of Q(x,z,theta)=z*cos(x+theta)^2, to construct a symbol */
 for (l=0;l <k;l++)
   { printf("l= %d \n",l);
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {
	 x=xmin+m*dx;
	 z=zmin+n*dz;
	 /*
	  if(l==0)
	 FQ[(k-1)/2+l+k*(m*nz+n)]=0;
	 else
	 FQ[(k-1)/2+l+k*(m*nz+n)]=(0+I*0.5);
	 */
	 
	 /* if (l ==-1 || l==1)
	 FQ[(k-1)/2+l+k*(m*nz+n)]=0;
         */
	 if (l ==1)
	   FQ[l+k*(m*nz+n)]=z*cpow(exp(1),I*2*l*x)/4.0;//cexp causing trouble use cpow instead!
	 if (l==0)
	   FQ[l+k*(m*nz+n)]=z*0.5;
	 

	  // printf("%3.3f ",creal(FQ[(k-1)/2+l+k*(m*nz+n)]));printf("+I %3.3f",cimag(FQ[(k-1)/2+l+k*(m*nz+n)]));printf("   ");
	
       }//printf("\n");
   }
   }
 /*
 for ( l=0;l<k;l++)
   {
   printf (" l = %d \n " ,l);
  for( m=0;m<nx;m++)
    {
      for(n=0;n<nz;n++)
	{
	  printf("%3.3f ",creal(FQ[l+k*(m*nz+n)]));printf("+I %3.3f",cimag(FQ[l+k*(m*nz+n)]));printf("   ");
	}
	printf("\n");
    }
   }
 */
 printf("setting FQ ok ... \n");

 /* QUtest=2*z*(2*cos(x)^2*x^2-4*z*sin(x)*cos(x)*x+2*z^2-2*z%2*cos(x)^2-1)*exp(-x^2 - z^2) the derivative in x */
 
 
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {
	 x=xmin+m*dx;
	 z=zmin+n*dz;
	 //QUtest[m*nz+n]=-2*x*U[m*nz+n];
	 // QUtest[m*nz+n]=-2*x*exp(-pow(x,2)-pow(z,2));
	 QUtest[m*nz+n] = z*2*(2*cos(x)*cos(x)*x*x-4*z*sin(x)*cos(x)*x+2*z*z-2*z*z*cos(x)*cos(x)-1)*exp(-x*x-z*z);
       }
   }

 printf("setting QUtest ok ...\n");
 
 clock_t start,end;
 double elapsed;
 start = clock();
 PsiDOf_real(U,FQ,QU,nx,nz,k,dx,dz,m_ord);
 end=clock();
 elapsed = ((double)(end-start))/CLOCKS_PER_SEC;
 printf(" time taken: %lf \n",elapsed);
 printf("PsiDO ok... \n");
 
 float error=0;


 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {
	 error = error + pow(cabs(QU[m*nz+n]-QUtest[m*nz+n]),2);
       }
   }
 float n_QUtest=0;
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {
	 n_QUtest = n_QUtest + pow(cabs(QUtest[m*nz+n]),2);
       }
   }

 printf("the relative error is %1.4e \n",pow(error/n_QUtest,0.5));

 /*  for( m=0;m<nx;m++)
    {
      for(n=0;n<nz;n++)
	{
	  printf("%f ",creal(QUtest[m*nz+n]));printf("+I %f",cimag(QUtest[m*nz+n]));printf("   ");
	}
	printf("\n");
	}
 */



 //cleaning
 free(U);free(FQ);free(QUtest);free(QU);//free(FQ1);free(FQ2);
  return; 
}

