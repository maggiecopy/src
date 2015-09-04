#include "conv.h"
#include <string.h>
#include <time.h>
//convolution of two vectors turns out to be useless, it is too slow to call 
//it in a loop to convolve along one dimension in a 3d vector.
void conv(fftw_complex const *x,fftw_complex const *y,fftw_complex *c,int nx,int ny)
     {

     fftw_complex *X,*Y;
     fftw_complex *FX,*FY,*FC;
     fftw_plan planx,plany,planc;

     X =(fftw_complex *) fftw_malloc((nx+ny-1)*sizeof(fftw_complex));
     Y= (fftw_complex *) fftw_malloc((nx+ny-1)*sizeof(fftw_complex));
     FX =(fftw_complex *) fftw_malloc((nx+ny-1)*sizeof(fftw_complex));
     FY =(fftw_complex *) fftw_malloc((nx+ny-1)*sizeof(fftw_complex));
     FC =(fftw_complex *) fftw_malloc((nx+ny-1)*sizeof(fftw_complex));

     planx = fftw_plan_dft_1d(nx+ny-1,X,FX,FFTW_FORWARD,FFTW_ESTIMATE);
     plany = fftw_plan_dft_1d(nx+ny-1,Y,FY,FFTW_FORWARD,FFTW_ESTIMATE);
     planc = fftw_plan_dft_1d(nx+ny-1,FC,c,FFTW_BACKWARD,FFTW_ESTIMATE);

     int m;
     for(m=0;m<nx;m++){
       X[m]=x[m];
     }

     fftw_execute(planx);

     for(m=0;m<ny;m++){
       Y[m]=y[m];
     }

     fftw_execute(plany);

     for(m=0;m<nx+ny-1;m++){
       FC[m] = FX[m]*FY[m]/(nx+ny-1);
     }

     fftw_execute(planc);

     //cleaning
     fftw_destroy_plan(planx);
     fftw_destroy_plan(plany);
     fftw_destroy_plan(planc);
     fftw_free(X);
     fftw_free(Y);
     fftw_free(FX);
     fftw_free(FY);
     fftw_free(FC);
     }

/*
void Qsquare(complex const *FQ,complex *FQ2,int nx,int nz,int k){

  complex *cl,*cl2;

  cl = (complex *) malloc(k*sizeof(complex));
  cl2 = (complex *) malloc((2*k-1)*sizeof(complex));

  int m,n,l;
 for (m=0;m<nx;m++)
	{
	  for(n=0;n<nz;n++)
	    {
	      for(l=0;l<k;l++)
		{
		  cl[l]=FQ[l+k*(m*nz+n)];
		}
	      conv(cl,cl,cl2,k,k);
	      for(l=0;l<2*k-1;l++)
		{
		  FQ2[l+(2*k-1)*(m*nz+n)]=cl2[l];
		}
	    }
	}

 free(cl);free(cl2);
}
*/

//squares the symbol by autoconvolving the fourier coefficients
//size of input FQ: nx*nz*k, size of outpout FQ2:nx*nz*(2k-1)
void Qsquare(complex const *FQ,complex *FQ2,int nx,int nz,int k){

  complex *FQ1;
  FQ1 = (complex *) calloc((2*k-1)*nx*nz,sizeof(complex));

  int m,n,l;
  // printf(" k = %d  \n",k);
  for (l=0;l<k;l++)
    {//printf("l= %d\n",l);
      for (m=0;m<nx;m++)
	{
	  for(n=0;n<nz;n++)
	    {//printf("%d ",l+(k)*(m*nz+n));
	    FQ1[l+(2*k-1)*(m*nz+n)] = FQ[l+k*(m*nz+n)]; //printf("set \n");
	    }
	}
    }


  // printf("padded array set in Qsquare \n");

  int rank =1;
  int *dims;
  dims = (int*) malloc(rank*sizeof(int)); 
  dims[0]=(2*k-1); //dims[1]=nz; dims[2]=2*k-1;
  int howmany = nx*nz;
  //inembed=onembed=dims

  int istride = 1;
  int ostride = 1;
  int idist = (2*k-1);
  int odist = idist;
  // printf("variables set ...\n");
  /*
 clock_t start,end;
 double elapsed;
 start=clock();
  */
  fftw_plan plan_forward,plan_backward;
  plan_forward = fftw_plan_many_dft(rank,dims,howmany,FQ1,NULL,istride,idist,FQ1,NULL,ostride,odist,FFTW_FORWARD,FFTW_ESTIMATE);
  plan_backward = fftw_plan_many_dft(rank,dims,howmany,FQ1,NULL,istride,idist,FQ2,NULL,ostride,odist,FFTW_BACKWARD,FFTW_ESTIMATE);

  fftw_execute(plan_forward);
  //printf("plan executed \n");

 for (l=0;l<2*k-1;l++)
    {//printf("l= %d\n",l);
      for (m=0;m<nx;m++)
	{
	  for(n=0;n<nz;n++)
	    {//printf("%d ",l+(k)*(m*nz+n));
	    FQ1[l+(2*k-1)*(m*nz+n)] = FQ1[l+(2*k-1)*(m*nz+n)]*FQ1[l+(2*k-1)*(m*nz+n)]/(2*k-1); //printf("set \n");
	    }
	}
    }
 /*
 end=clock();
 elapsed = ((double)(end-start))/CLOCKS_PER_SEC;
 printf(" time taken to set up plans in Qsquare: %lf \n",elapsed);
 */
 fftw_execute(plan_backward);

 //cleaning
 fftw_destroy_plan(plan_forward);
 fftw_destroy_plan(plan_backward);
 free(FQ1);


}
//constructs a symbol that has conjugate symmetry from the first half of the coefficients
//again size of input FQ: nx*nz*k , output FQ2:nx*nz*(2k-1)
void Qsym(complex const *FQ,complex *FQ2,int nx,int nz,int k){
  int m,n,l;
            
 for (l=0;l<k;l++)
    {//printf("l= %d\n",l);
      for (m=0;m<nx;m++)
	{
	  for(n=0;n<nz;n++)
	    {
	      FQ2[l+(k-1)+(2*k-1)*(m*nz+n)]=FQ[l+(k)*(m*nz+n)];
	      if( l !=0)
	      FQ2[(k-1)-l+(2*k-1)*(m*nz+n)]=conj(FQ[l+(k)*(m*nz+n)]);
	    }
	}
    }
}


//Shifts the coefficients by shift entries, 
//size of input FQ: nx*nz*k size of output:FQ2: nx*nz*(2k-1)
//the zero shift leaves FQ in the center but pads with (k-1)/2 zeros on each side
void Qshift(complex const *FQ, complex *FQ2, int nx, int nz, int k, int shift){

  //make sure FQ2 is zeroed before filling
 

  int m,n,l;
  /*
  for (m=0;m<(2*k-1)*nx*nz;m++)
    { FQ2[m]=0;}
  */
  memset(FQ2,0,(2*k-1)*nx*nz*sizeof(complex));
  if (shift >(k-1)/2 || shift <-(k-1)/2)
    printf("Illegal Shift,returning zeros\n");
else
  {
    // printf("setting stuff\n");
 for (l=0;l<k;l++)
    {//printf("l= %d\n",l);
      for (m=0;m<nx;m++)
	{
	  for(n=0;n<nz;n++)
	    {

	      FQ2[(k-1)/2+shift+ l+(2*k-1)*(m*nz+n)] = FQ[l+k*(m*nz+n)];
	    }
	}
    }
  }
  
}





