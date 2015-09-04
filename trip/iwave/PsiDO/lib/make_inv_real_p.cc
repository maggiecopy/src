using namespace std;

#include <iostream>
#include <stdio.h>
#include <string.h>
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

void make_inv_real_p(fftwf_complex const *Np_b, fftwf_complex * inv, float const *coeff_c, int nxc, int nzc, int k,int p, float x0, int nx, float dx, float z0, int nz, float dz, float m_ord){

  float xf = x0 + (nx-1)*dx;
  float zf = z0 + (nz-1)*dz;

  float x0c = x0;
  float z0c = z0;

  float dxc = (xf-x0)/(float)(nxc-1);
  float dzc = (zf-z0)/(float)(nzc-1);

  // int wl = nxc*nzc+max(nxc,nzc)+max(nzc,nx)+max(4*nzc-8,4*nxc-3);
 int wl = nxc*nzc+max(nxc,nzc)+max(nxc,nz)+max(4*nxc-8,4*nzc-3);
  int wl2 = max(nx,nz)+max(nxc,nzc)+ max(7*nxc-8+nx,7*nzc-8+nz)+nxc*nz;

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
  
  
  // fftwf_complex * QU = (fftwf_complex *)calloc(nx*nz,sizeof(fftwf_complex));

 // PsiDOf_real(U,FQ,inv,nx,nz,k,(double)dx,(double)dz,(double)m_ord);
  
  fftwf_complex * QU = (fftwf_complex *)calloc(nx*nz*p*p,sizeof(fftwf_complex));
  fftwf_complex * FQ_temp=(fftwf_complex *)calloc(nx*nz*k,sizeof(fftwf_complex));
  fftwf_complex * U_temp = (fftwf_complex *)calloc(nx*nz,sizeof(fftwf_complex));
  fftwf_complex * QU_temp = (fftwf_complex *)calloc(nx*nz,sizeof(fftwf_complex));

  for(ip=0;ip<p;ip++)
    {

      memcpy(FQ_temp,&FQ[(ip)*nx*nz*k],nx*nz*k*sizeof(fftwf_complex));

      for(ip2=0;ip2<p;ip2++)
	{	  
          memcpy(U_temp,&Np_b[ip2*(p+1)*nx*nz+(ip)*nx*nz],nx*nz*sizeof(fftwf_complex));

	  PsiDOf_real(U_temp,FQ_temp,QU_temp,nx,nz,k,(double)dx,(double)dz,(double)(ip+1)*m_ord);

	  memcpy(&QU[ip2*p*nx*nz+(ip)*nx*nz],QU_temp,nx*nz*sizeof(fftwf_complex));
	}
    
    }
  fftwf_complex *J=(fftwf_complex*)calloc(nx*nz*p,sizeof(fftwf_complex));//zeros the bytes so OK
  for(ip2=0;ip2<p;ip2++)
    {
      //memcpy(&J[ip2*nx*nz],&Np_b[ip2*(p+1)*nx*nz],nx*nz*sizeof(fftwf_complex));

      for(ip=0;ip<p;ip++)
	{
	  for (m=0;m<nx;m++)
	    {
	      for (n=0;n<nz;n++)
		{ 
		  J[ip2*nx*nz+m*nz+n][0]=J[ip2*nx*nz+m*nz+n][0]+QU[ip2*p*nx*nz+ip*nx*nz+m*nz+n][0];
		  J[ip2*nx*nz+m*nz+n][1]=J[ip2*nx*nz+m*nz+n][1]+QU[ip2*p*nx*nz+ip*nx*nz+m*nz+n][1];
		}
	    }
	}
    }
  memcpy(inv,J,nx*nz*p*sizeof(fftwf_complex));

  free(FQ);free(FQtemp_real);free(FQtemp_imag);free(work);free(work2);free(QU);free(FQ_temp);free(U_temp);free(QU_temp);free(J);
}
