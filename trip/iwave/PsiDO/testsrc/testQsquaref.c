#include "conv.h"
#include <time.h>


//#define VERBOSE 1

int main(int argc, char **argv) {
  printf("starting...\n");

 int nx =1001;
 int nz=1001;
 int k=5;

 double xmin=-5.0;
 double xmax = 5.0;
 double zmin=-5.0;
 double zmax = 5.0;
 double x,z;

 double dx=(xmax-xmin)/((double)(nx-1));
 double dz = (zmax-zmin)/((double)(nz-1));
 double m_ord=0; /*for cos(x+theta)^2 */

 float complex *FQ,*FQ1,*FQ2;
 int k1 = (k-1)/2+1;

 FQ = (float complex *)  malloc(nx*nz*k*sizeof(float complex));
 FQ1 = (float complex *) malloc(nx*nz*k1*sizeof(float complex));
 FQ2 = (float complex *) malloc(nx*nz*k*sizeof(float complex));

 int m,n,l;
 // Fourier Expansion of Q(x,z,theta)=cos(x+theta)^2, to construct a symbol 
 for (l= -(k-1)/2;l <=(k-1)/2;l++)
   {// printf("l= %d \n",l);
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {
	 x=xmin+m*dx;
	 z=zmin+n*dz;
	 /* if(l==0)
	 FQ[(k-1)/2+l+k*(m*nz+n)]=0;
	 else
	 FQ[(k-1)/2+l+k*(m*nz+n)]=(0+I*1/2);*/
	 
	 if (l ==-1 || l==1)
	 FQ[(k-1)/2+l+k*(m*nz+n)]=0;
	 if (l ==2 || l==-2)
	   FQ[(k-1)/2+l+k*(m*nz+n)]=cpow(exp(1),I*l*x)/4.0;//cexp causing trouble use cpow instead!
	 if (l==0)
	   FQ[(k-1)/2+l+k*(m*nz+n)]=0.5;
	 //#IFDEF VERBOSE
	 //printf("%3.3f ",creal(FQ[(k-1)/2+l+k*(m*nz+n)]));printf("+I %3.3f",cimag(FQ[(k-1)/2+l+k*(m*nz+n)]));printf("   ");
	 //#ENDIF	
       } //printf("\n");
   }
   }


 //Fourier expansion of cos(x+theta)

for (l= -(k1-1)/2;l <=(k1-1)/2;l++)
  {// printf("l= %d \n",l);
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {
	 x=xmin+m*dx;
	 z=zmin+n*dz;
	 /* if(l==0)
	 FQ[(k-1)/2+l+k*(m*nz+n)]=0;
	 else
	 FQ[(k-1)/2+l+k*(m*nz+n)]=(0+I*1/2);*/
	 
	 if (l ==0)
	 FQ1[(k1-1)/2+l+k1*(m*nz+n)]=0;
	 else
	   FQ1[(k1-1)/2+l+k1*(m*nz+n)]=cpow(exp(1),I*l*x)/2.0;//cexp causing trouble use cpow instead!
	 

	 //printf("%3.3f ",creal(FQ1[(k1-1)/2+l+k1*(m*nz+n)]));printf("+I %3.3f",cimag(FQ1[(k1-1)/2+l+k1*(m*nz+n)]));printf("   ");
	
       } //printf("\n");
   }
   }

//Now square FQ1
 clock_t start,end;
 double elapsed;
 printf("calling Qsquare \n");
 start = clock();
 Qsquaref(FQ1,FQ2,nx,nz,k1);
 end=clock();
 elapsed = ((double)(end-start))/CLOCKS_PER_SEC;
 printf(" time taken: %lf \n",elapsed);

double error = 0;

for (l= -(k-1)/2;l <=(k-1)/2;l++)
  { //printf("l=%d \n",l);
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {//printf("%3.3f ",creal(FQ2[(k-1)/2+l+k*(m*nz+n)]));printf("+I %3.3f",cimag(FQ2[(k-1)/2+l+k*(m*nz+n)]));printf("   ");

	 error = error + pow(cabs(FQ2[(k-1)/2+l+k*(m*nz+n)]-FQ[(k-1)/2+l+k*(m*nz+n)]),2);
       }//printf("\n");
   }
   }
 printf("the error in squaring FQ1 to obtain FQ is %1.4e \n",pow(error,0.5));

 //Begin test of Qsym
 float complex *FQ_half,*FQ_sym;
 FQ_half = (float complex*)malloc(nx*nz*k1*sizeof(float complex));
 FQ_sym = (float complex*)malloc(nx*nz*k*sizeof(float complex));
 //build half of the coefficients of cos(x+theta)^2
 
 for (l= 0;l <k1;l++)
   {// printf("l= %d \n",l);
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {
	 x=xmin+m*dx;
	 z=zmin+n*dz;
	
	 
	 if (l ==-1 || l==1)
	 FQ_half[l+k1*(m*nz+n)]=0;
	 if (l ==2 || l==-2)
	   FQ_half[l+k1*(m*nz+n)]=cpow(exp(1),I*l*x)/4.0;//cexp causing trouble use cpow instead!
	 if (l==0)
	   FQ_half[l+k1*(m*nz+n)]=0.5;
	 //#IFDEF VERBOSE
	 //printf("%3.3f ",creal(FQ[(k-1)/2+l+k*(m*nz+n)]));printf("+I %3.3f",cimag(FQ[(k-1)/2+l+k*(m*nz+n)]));printf("   ");
	 //#ENDIF	
       } //printf("\n");
   }
   }
 
 /*
 float complex *FQ_half2 = (float complex *) malloc(nx*nz*k1*sizeof(float complex));
 for (l= 0;l <k1;l++)
   {// printf("l= %d \n",l);
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {
	 x=xmin+m*dx;
	 z=zmin+n*dz;
	 
	 if (l ==-1 || l==1)
	 FQ_half[l+k1*(m*nz+n)]=0;
	 if (l ==2 || l==-2)
	   FQ_half2[l+k1*(m*nz+n)]=cpowf(exp(1),I*l*x)/4.0;//cexp causing trouble use cpow instead!
	 if (l==0)
	   FQ_half2[l+k1*(m*nz+n)]=0.5;
	 //#IFDEF VERBOSE
	 //printf("%3.3f ",creal(FQ[(k-1)/2+l+k*(m*nz+n)]));printf("+I %3.3f",cimag(FQ[(k-1)/2+l+k*(m*nz+n)]));printf("   ");
	 //#ENDIF	
       } //printf("\n");
   }
   }
 printf("FQhalf2 set\n");
 FQ_half = (double complex *) FQ_half2;

 */

 //Now symmetrize
 Qsymf(FQ_half,FQ_sym,nx,nz,k1);

 error = 0;

for (l= -(k-1)/2;l <=(k-1)/2;l++)
  {// printf("l=%d \n",l);
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {//printf("%3.3f ",creal(FQ_sym[(k-1)/2+l+k*(m*nz+n)]));printf("+I %3.3f",cimag(FQ_sym[(k-1)/2+l+k*(m*nz+n)]));printf("   ");

	 error = error + pow(cabs(FQ_sym[(k-1)/2+l+k*(m*nz+n)]-FQ[(k-1)/2+l+k*(m*nz+n)]),2);
       }//printf("\n");
   }
   }
 printf("the error in symmetrizing FQ_half to obtain FQ is %1.4e \n",pow(error,0.5));

 //Begin test of Qshift
 

 float complex *FQ3 = (float complex*)malloc(nx*nz*(2*k-1)*sizeof(float complex)); FQ3[0]=-13;
 Qshiftf(FQ,FQ3,nx,nz,k,-1);
 /*
for (l= -(2*k-2)/2;l <=(2*k-2)/2;l++)
  {printf("l=%d \n",l);
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {printf("%3.3f ",creal(FQ3[(2*k-2)/2+l+(2*k-1)*(m*nz+n)]));printf("+I %3.3f",cimag(FQ3[(2*k-2)/2+l+(2*k-1)*(m*nz+n)]));printf("   ");


       }printf("\n");
   }
   }
 */

 /*
 double *re_FQ = (double *) malloc(nx*nz*k*sizeof(double));
 double *im_FQ = (double *) malloc(nx*nz*k*sizeof(double));

for (l= 0;l <k;l++)
  {
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {
	 re_FQ[l+k*(m*nz+n)]= creal(FQ[l+k*(m*nz+n)]);
	 im_FQ[l+k*(m*nz+n)]= cimag(FQ[l+k*(m*nz+n)]);
       }
   }
  }


for (l= 0;l <k;l++)
  {printf("l=%d \n",l);
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {printf("%3.3f ",creal(re_FQ[l+k*(m*nz+n)]));printf("+I %3.3f",cimag(re_FQ[l+k*(m*nz+n)]));printf("   ");


       }printf("\n");
   }
   }
for (l= 0;l <k;l++)
  {printf("l=%d \n",l);
 for (m=0;m<nx;m++)
   {
     for (n=0;n<nz;n++)
       {printf("%3.3f ",creal(FQ[l+k*(m*nz+n)]));printf("+I %3.3f",cimag(FQ[l+k*(m*nz+n)]));printf("   ");


       }printf("\n");
   }
   }
 */
 //cleaning
 free(FQ);free(FQ1);free(FQ2);free(FQ_half);free(FQ_sym);free(FQ3);//free(re_FQ);free(im_FQ);

}
