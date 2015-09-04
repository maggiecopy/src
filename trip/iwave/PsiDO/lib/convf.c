#include "conv.h"
#include <string.h>
//same as conv.c but all operations are done in float precision,
//please consult conv.c for documentaion
/*
void conv(complex const *x,complex const *y,complex *c,int nx,int ny)
     {

     complex *X,*Y;
     fftw_complex *FX,*FY,*FC;
     fftw_plan planx,plany,planc;

     X =(complex *) calloc(nx+ny-1,sizeof(complex));
     Y= (complex *) calloc(nx+ny-1,sizeof(complex));
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
     free(X);
     free(Y);
     fftw_free(FX);
     fftw_free(FY);
     fftw_free(FC);
     } */

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


void Qsquaref(fftwf_complex const *FQ,fftwf_complex *FQ2,int nx,int nz,int k){

 fftwf_complex *FQ1;
  FQ1 = (fftwf_complex *) fftwf_malloc((2*k-1)*nx*nz*sizeof(fftwf_complex));

  int m,n,l;
  // printf(" k = %d  \n",k);
  for (l=0;l<k;l++)
    {//printf("l= %d\n",l);
      for (m=0;m<nx;m++)
	{
	  for(n=0;n<nz;n++)
	    {//printf("%d ",l+(k)*(m*nz+n));
	    FQ1[l+(2*k-1)*(m*nz+n)] = FQ[l+k*(m*nz+n)]; 
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

  fftwf_plan plan_forward,plan_backward;
  plan_forward = fftwf_plan_many_dft(rank,dims,howmany,FQ1,NULL,istride,idist,FQ1,NULL,ostride,odist,FFTW_FORWARD,FFTW_ESTIMATE);
  plan_backward = fftwf_plan_many_dft(rank,dims,howmany,FQ1,NULL,istride,idist,FQ2,NULL,ostride,odist,FFTW_BACKWARD,FFTW_ESTIMATE);

  fftwf_execute(plan_forward);
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

 fftwf_execute(plan_backward);

 //cleaning
 fftwf_destroy_plan(plan_forward);
 fftwf_destroy_plan(plan_backward);
 fftwf_free(FQ1);


}

void Qsymf(fftwf_complex const *FQ,fftwf_complex *FQ2,int nx,int nz,int k){
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



void Qshiftf(fftwf_complex const *FQ, fftwf_complex *FQ2, int nx, int nz, int k, int shift){


  //make sure FQ2 is zeroed before filling
  

  int m,n,l;
  /*
  for (m=0;m<(2*k-1)*nx*nz;m++)
    { FQ2[m]=0;}
  */
  memset(FQ2,0,(2*k-1)*nx*nz*sizeof(fftwf_complex));
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




