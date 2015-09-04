using namespace std;

#include <stdio.h>
#include <string.h>
#include <iostream>
#include "objective_real.h"
//#include "PsiDOf.h"


//int objective_real(complex const *U,float const *coeff_c,int nxc,int nzc,int k,double x0, int nx , double dx, double z0, int nz, double dz, double m_ord, float Obj, float *Grad);


extern "C" { 

  void bicubic_(float& , float& , int& ,
		float&, float& , int&,  
		float*,
		float& , float& , int& , 
		float& , float&, int&, 
		float*,
		int&, float*,
		int&, int&);
}

extern "C" { 

  void bicubicadj_(float& , float& , int& ,
		   float&, float& , int&, 
		   float*,
		   float& , float& , int& , 
		   float& , float&, int&, 
		   float*,
		   int&, float*,
		   int&, int&);
}

extern "C"
{

  void PsiDOf_real(fftwf_complex  const*, fftwf_complex  const*, fftwf_complex *, int, int, int, double, double, double);
}

extern "C"
{
  void PsiDOf_l(fftwf_complex const  *,fftwf_complex const *,fftwf_complex * ,int , int , double ,double ,double ,int );
}

void objective_real_p(fftwf_complex const *Np_b, float const *coeff_c, int nxc, int nzc, int k, int p, float x0, int nx, float dx, float z0, int nz, float dz, float m_ord, float *Obj, float *Grad,int flag_grad,float lam){

  float xf = x0 + (nx-1)*dx;
  float zf = z0 + (nz-1)*dz;

  float x0c = x0;
  float z0c = z0;

  float dxc = (xf-x0)/(float)(nxc-1);
  float dzc = (zf-z0)/(float)(nzc-1);

  //int wl = nxc*nzc+max(nxc,nzc)+max(nzc,nx)+max(4*nzc-8,4*nxc-3);
  int wl = nxc*nzc+max(nxc,nzc)+max(nxc,nz)+max(4*nxc-8,4*nzc-3);
  //int wl2 = max(nx,nz)+max(nxc,nzc)+ max(7*nxc-8+nx,7*nzc-8+nz)+nxc*nz;
  int wl2 = max(nx,nz)+max(nxc,nzc)+ max(7*nxc-8+nx,7*nzc-8+nz)+nzc*nx;
  float * work = (float *)malloc(wl *sizeof(float));
  float * work2 = (float *)malloc(wl2 *sizeof(float));

  int ierr=0;
  int iend=1;

  fftwf_complex * FQ =  (fftwf_complex *)malloc(nx*nz*k*p*sizeof(fftwf_complex));

  float * FQtemp_real = (float *)malloc(nx*nz*sizeof(float));
  float * FQtemp_imag = (float *)malloc(nx*nz*sizeof(float));
  float * coeff_temp = (float *)malloc(nxc*nzc*sizeof(float));

  int l,m,n,ip,ip2;

  for (ip=0;ip<p;ip++)
    {

      for (l=0;l<k;l++)
	{ 
	  if(l==0)
	    {
	      memcpy(coeff_temp,&coeff_c[ip*nxc*nzc*(2*k-1)+0],nxc*nzc*sizeof(float));

	      //printf("obj fct calculation 0 \n");
	      bicubic_(z0c,dzc,nzc,x0c,dxc,nxc,coeff_temp,z0,dz,nz,x0,dx,nx,FQtemp_real,iend,work,wl,ierr);
	      // bicubic_(x0c,dxc,nxc,z0c,dzc,nzc,coeff_temp,x0,dx,nx,z0,dz,nz,FQtemp_real,iend,work,wl,ierr);
	      // for (m=0;m<nxc*nzc;m++){ printf("%f",coeff_temp[m]);printf("    ");} 
	      for (m=0;m<nx;m++)
		{
		  for (n=0;n<nz;n++)
		    {
		      FQ[ip*nx*nz*k+k*(m*nz+n)][0]=FQtemp_real[(m*nz+n)]; 
		      FQ[ip*nx*nz*k+k*(m*nz+n)][1]=0.f; 
		      // FQ[l+k*(m*nz+n)]=coeff_temp[l+k*(m*nz+n)]+I*0;

		    }
		}
	    }
	  else
	    {// printf("obj fct calculation 1 \n");
	      memcpy(coeff_temp,&coeff_c[ip*nxc*nzc*(2*k-1)+nxc*nzc*l],nxc*nzc*sizeof(float)); 

	      bicubic_(z0c,dzc,nzc,x0c,dxc,nxc,coeff_temp,z0,dz,nz,x0,dx,nx,FQtemp_real,iend,work,wl,ierr);
	      //bicubic_(x0c,dxc,nxc,z0c,dzc,nzc,coeff_temp,x0,dx,nx,z0,dz,nz,FQtemp_real,iend,work,wl,ierr);

	      memcpy(coeff_temp,&coeff_c[ip*nxc*nzc*(2*k-1)+nxc*nzc*k+nxc*nzc*(l-1)],nxc*nzc*sizeof(float));
	      bicubic_(z0c,dzc,nzc,x0c,dxc,nxc,coeff_temp,z0,dz,nz,x0,dx,nx,FQtemp_imag,iend,work,wl,ierr);
	      // bicubic_(x0c,dxc,nxc,z0c,dzc,nzc,coeff_temp,x0,dx,nx,z0,dz,nz,FQtemp_imag,iend,work,wl,ierr);

	      for (m=0;m<nx;m++)
		{
		  for (n=0;n<nz;n++)
		    {
		      FQ[ip*nx*nz*k+l+k*(m*nz+n)][0]=FQtemp_real[(m*nz+n)];
		      FQ[ip*nx*nz*k+l+k*(m*nz+n)][1]=FQtemp_imag[(m*nz+n)];
		    }
		}
	    }
	}
    }
  /*
    for (l=0;l<k;l++)
    {  printf("l= %d \n",l);
    for (m=0;m<nx;m++)
    {
    for (n=0;n<nz;n++)
    {
    printf("%3.3f ",creal(FQ[l+k*(m*nz+n)]));
    printf("+I %3.3f",cimag(FQ[l+k*(m*nz+n)]));
    printf("   ");
	
    }
    printf("\n");
    }
    }
  */
  
  
  fftwf_complex * QU = (fftwf_complex *)calloc(nx*nz*p*p,sizeof(fftwf_complex));
  fftwf_complex * FQ_temp=(fftwf_complex *)calloc(nx*nz*k,sizeof(fftwf_complex));
  fftwf_complex * U_temp = (fftwf_complex *)calloc(nx*nz,sizeof(fftwf_complex));
  fftwf_complex * QU_temp = (fftwf_complex *)calloc(nx*nz,sizeof(fftwf_complex));
 
  for(ip=0;ip<p;ip++)
    {
      //      for (l=0;l<k;l++)
      //        { 
      // 	 for (m=0;m<nx;m++)
      // 	   {
      // 	     for (n=0;n<nz;n++)
      // 	       {
      // 		 FQ_temp[l+k*(m*nz+n)]=FQ[(ip-1)*nx*nz*(2*k-1)+l+k*(m*nz+n)];
      // 		 U_temp[l+k*(m*nz+n)]=FQ[(ip-1)*nx*nz*(2*k-1)+l+k*(m*nz+n)
      // 		   }
      // 	   }
      //        }
      memcpy(FQ_temp,&FQ[(ip)*nx*nz*k],nx*nz*k*sizeof(fftwf_complex));

      for(ip2=0;ip2<p;ip2++)
	{	  
          memcpy(U_temp,&Np_b[ip2*(p+1)*nx*nz+(ip+1)*nx*nz],nx*nz*sizeof(fftwf_complex));

	  PsiDOf_real(U_temp,FQ_temp,QU_temp,nx,nz,k,(double)dx,(double)dz,(double)(ip+1)*m_ord);

	  memcpy(&QU[ip2*p*nx*nz+(ip)*nx*nz],QU_temp,nx*nz*sizeof(fftwf_complex));
	}
    
    }

  Obj[0]=0;
  float alpha;

  fftwf_complex *J=(fftwf_complex*)calloc(nx*nz*p,sizeof(fftwf_complex));
  for(ip2=0;ip2<p;ip2++)
    {
      memcpy(&J[ip2*nx*nz],&Np_b[ip2*(p+1)*nx*nz],nx*nz*sizeof(fftwf_complex));

      for(ip=0;ip<p;ip++)
	{
	  for (m=0;m<nx;m++)
	    {
	      for (n=0;n<nz;n++)
		{ 
		  J[ip2*nx*nz+m*nz+n][0]=J[ip2*nx*nz+m*nz+n][0]-QU[ip2*p*nx*nz+ip*nx*nz+m*nz+n][0];
		  J[ip2*nx*nz+m*nz+n][1]=J[ip2*nx*nz+m*nz+n][1]-QU[ip2*p*nx*nz+ip*nx*nz+m*nz+n][1];
		}
	    }
	}
    }
  for(ip2=0;ip2<p;ip2++)
    {
      for (m=0;m<nx;m++)
	{
	  for (n=0;n<nz;n++)
	    { 
//	      alpha = cabs(J[ip2*nx*nz+m*nz+n]);
	      Obj[0] = Obj[0] + J[ip2*nx*nz+m*nz+n][0]*J[ip2*nx*nz+m*nz+n][0]
                              + J[ip2*nx*nz+m*nz+n][1]*J[ip2*nx*nz+m*nz+n][1]; //alpha*alpha;
	      // printf("Obj= %f \n",Obj[0]);
	    }
	}
    }
  for (m=0;m<nxc*nzc*(2*k-1);m++)
    {
      Obj[0]=Obj[0]+lam*coeff_c[m]*coeff_c[m];
    }

  //printf(" Obj in objective = %f \n", Obj[0]);
  
  if(flag_grad)
    {// printf("In Grad");
      fftwf_complex * ones = ( fftwf_complex * ) malloc(nx*nz*sizeof(fftwf_complex));
      fftwf_complex * QUtemp = (fftwf_complex *) malloc(nx*nz*sizeof(fftwf_complex));
      fftwf_complex * QUtemp2 = (fftwf_complex *) malloc(nx*nz*sizeof(fftwf_complex));
      float *Grad2 = (float *)malloc(nxc*nzc*(2*k-1)*p*sizeof(float));
  
      for (m=0;m<nx;m++)
	{
	  for (n=0;n<nz;n++)
	    { 
	      ones[m*nz+n][0]=1.f; 
	      ones[m*nz+n][1]=0.f; 
	    }
	}
      for(ip2=0;ip2<p;ip2++)
	{
	  for(ip=0;ip<p;ip++)
	    { 
   
	      memcpy(U_temp,&Np_b[ip2*(p+1)*nx*nz+(ip+1)*nx*nz],nx*nz*sizeof(fftwf_complex));
      
	      for(l=0;l<k;l++)
		{ 

		  if (l==0)
		    {
	  
		      PsiDOf_l(U_temp ,ones,QUtemp ,nx, nz, dx,dz,(ip+1)*m_ord,2*l);
		      for (m=0;m<nx;m++)
			{
			  for (n=0;n<nz;n++)
			    {
			      FQtemp_real[m*nz+n]=2*dxc/dx*dzc/dz*
                                     ( QUtemp[m*nz+n][0]*J[ip2*nx*nz+m*nz+n][0]
                                      +QUtemp[m*nz+n][1]*J[ip2*nx*nz+m*nz+n][1]);
//			      FQtemp_real[m*nz+n]=2*dxc/dx*dzc/dz*creal(QUtemp[m*nz+n]*(creal(J[ip2*nx*nz+m*nz+n])-I*cimag(J[ip2*nx*nz+m*nz+n])));
			    }
			}
		      bicubicadj_(z0,dz,nz,x0,dx,nx,FQtemp_real,z0c,dzc,nzc,x0c,dxc,nxc,coeff_temp,iend,work2,wl2,ierr);

		      // bicubicadj_(x0,dx,nx,z0,dz,nz,FQtemp_real,x0c,dxc,nxc,z0c,dzc,nzc,coeff_temp,iend,work2,wl2,ierr);
		      memcpy(&Grad2[ip*nxc*nzc*(2*k-1)+0],coeff_temp,nxc*nzc*sizeof(float));
		      //need to scale by dx/dxc*dz/dzc
		      //printf("  %d",ip*nxc*nzc*(2*k-1)+0);
		    }

		  else
		    {
	 
	  
		      PsiDOf_l(U_temp,ones,QUtemp,nx,nz,dx,dz,(ip+1)*m_ord,2*l);
		      PsiDOf_l(U_temp,ones,QUtemp2,nx,nz,dx,dz,(ip+1)*m_ord,-2*l);
		      for (m=0;m<nx;m++)
			{
			  for (n=0;n<nz;n++)
			    {
			      FQtemp_real[m*nz+n]=2*dxc/dx*dzc/dz*
                                     ((QUtemp[m*nz+n][0]+QUtemp2[m*nz+n][0])*J[ip2*nx*nz+m*nz+n][0]
                                     +(QUtemp[m*nz+n][1]+QUtemp2[m*nz+n][1])*J[ip2*nx*nz+m*nz+n][1]);
//			      FQtemp_real[m*nz+n]=2*dxc/dx*dzc/dz*creal((QUtemp[m*nz+n]+QUtemp2[m*nz+n])*(creal(J[ip2*nx*nz+m*nz+n])-I*cimag(J[ip2*nx*nz+m*nz+n])));

			      FQtemp_imag[m*nz+n]=2*dxc/dx*dzc/dz*
                                     ((QUtemp2[m*nz+n][1]-QUtemp[m*nz+n][1])*J[ip2*nx*nz+m*nz+n][0]
                                     -(QUtemp2[m*nz+n][0]-QUtemp[m*nz+n][0])*J[ip2*nx*nz+m*nz+n][1]);
//			      FQtemp_imag[m*nz+n]=2*dxc/dx*dzc/dz*cimag((QUtemp2[m*nz+n]-QUtemp[m*nz+n])*(creal(J[ip2*nx*nz+m*nz+n])-I*cimag(J[ip2*nx*nz+m*nz+n])));
			      //FQtemp_imag[m*nz+n]=cimag(QUtemp[m*nz+n]);
			    }
			}
		      bicubicadj_(z0,dz,nz,x0,dx,nx,FQtemp_real,z0c,dzc,nzc,x0c,dxc,nxc,coeff_temp,iend,work2,wl2,ierr);
	  
		      //bicubicadj_(x0,dx,nx,z0,dz,nz,FQtemp_real,x0c,dxc,nxc,z0c,dzc,nzc,coeff_temp,iend,work2,wl2,ierr);
		      memcpy(&Grad2[ip*nxc*nzc*(2*k-1)+nxc*nzc*l],coeff_temp,nxc*nzc*sizeof(float));
		      // printf("  %d",ip*nxc*nzc*(2*k-1)+nxc*nzc*l);
		      bicubicadj_(z0,dz,nz,x0,dx,nx,FQtemp_imag,z0c,dzc,nzc,x0c,dxc,nxc,coeff_temp,iend,work2,wl2,ierr);
		      // bicubicadj_(x0,dx,nx,z0,dz,nz,FQtemp_imag,x0c,dxc,nxc,z0c,dzc,nzc,coeff_temp,iend,work2,wl2,ierr);
		      memcpy(&Grad2[ip*nxc*nzc*(2*k-1)+nxc*nzc*k+nxc*nzc*(l-1)],coeff_temp,nxc*nzc*sizeof(float));
		      // printf("  %d\n",ip*nxc*nzc*(2*k-1)+nxc*nzc*k+nxc*nzc*(l-1));
		    }
		}
	    }
	   if (ip2==0)
	     {
	      memcpy(&Grad[0],&Grad2[0],nxc*nzc*(2*k-1)*p*sizeof(float));
		    
	     }	 
	      else
		{
		  for (m=0;m<p*nxc*nzc*(2*k-1);m++)
		    {
		      Grad[m]=Grad[m]+Grad2[m];
		    }
		}
	}	    
	
	  //the minus sign fixes the sign left out in the gradient
	  for (m=0;m<p*nxc*nzc*(2*k-1);m++)
	    {
	      Grad[m]=-Grad[m]+2*lam*coeff_c[m];
	    }

		 
	
	  /*	 
  
	  for (m=0;m<nxc;m++)
	  {
	  for (n=0;n<nzc;n++)
	  {
	  printf("%3.3f ",coeff_temp[(m*nzc+n)]);
	  //printf("+I %3.3f",cimag(QU[(m*nz+n)]));
	  printf("   ");
	
	  }
	  printf("\n");
	  }
	  */
  
	  //printf("Gradient [0] = %6.6f \n",coeff_temp[2000]);
	  //printf(" dxc/dx = %1.6f \n",dxc/dx);
	  free(QUtemp);free(QUtemp2);free(ones);free(Grad2);
	}



      free(FQtemp_real);free(FQtemp_imag);free(coeff_temp);free(work);free(work2);free(QU); free(FQ);free(QU_temp);free(J);
    }

      
