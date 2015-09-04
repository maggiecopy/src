#include "PsiDO.h"
#include <time.h>
#include <string.h>
/*
Applies the action of a pseudodifferential operator ona  function U
parametres: U- input function
            FQ- Fourier coefficients of the symbol
	    QU- output
	    nx-size of first dimesion in U
	    nz-size of second dimension in U
	    k-number of fourier coefficients of the symbol
	    dx-spacing in x
	    dz-spacing in z
	    m_ord- the order of the psido
Note that all vectors are long vectors in row major form

*/

void PsiDO(fftwf_complex const  *U ,fftwf_complex const *FQ, fftwf_complex *QU ,int nx, int nz, int k, double dx,double dz,double m_ord){

  const double pi = 3.14159265358979323846;

  fftw_complex *FU,*FR,*R;
  fftw_plan plan_forward,plan_backward;
  //reinterpret_cast<fftw_complex*>(U);

  //in= (fftw_complex *) fftw_malloc(nx*nz*sizeof(fftw_complex));
  FU=(fftw_complex *) fftw_malloc(nx*nz*sizeof(fftw_complex));
  FR= (fftw_complex *) fftw_malloc(nx*nz*sizeof(fftw_complex));
  R=(fftw_complex *) fftw_malloc(nx*nz*sizeof(fftw_complex));
  /*
 clock_t start,end;
 double elapsed;
 start=clock();
  */
  plan_forward = fftw_plan_dft_2d(nx,nz,U,FU,FFTW_FORWARD,FFTW_ESTIMATE);
  plan_backward= fftw_plan_dft_2d(nx,nz,FR,R,FFTW_BACKWARD,FFTW_ESTIMATE);
  /*
  end=clock();
 elapsed = ((double)(end-start))/CLOCKS_PER_SEC;
 printf(" time taken to set up plans in PsiDO: %lf \n",elapsed);
  */
  /* reinterpret_cast <fftw_complex*>(U) may work here try it */

  int m,n,l,p,r;


  fftw_execute(plan_forward);

 

  double X = (nx-1)*dx;
  double Z= (nz-1)*dz;
  double dxi = 1.0/X;
  double deta =1.0/Z;
  double xi,eta,omega,theta,alpha;
  /*
  for(m=0;m<nx;m++)
    {
      for(n=0;n<nz;n++)
	{
	  QU[m*nz+n] = 0;
	}
	}*/
  //Initialize QU to zero
  memset(QU,0,nx*nz*sizeof(fftwf_complex));

  


  
  for (l=-(k-1)/2;l<=(k-1)/2;l++)
    {
      /*
      for(p = -(nx-1)/2;p<= (nx-1)/2;p++)
	{

	  xi = 2*pi*dxi*p;
	for(r = - (nz-1)/2; r <= (nz-1)/2; r++)
	  {
	    eta = 2*pi*deta*r;
	    omega = pow((xi*xi+eta*eta),(0.5));
	    
	    if (p !=0 || r !=0)
	      {
		if( p >=0 && r >=0)
		        FR[p*nz+r] = pow(omega,(m_ord-l))*cpow((xi+I*eta),l)*FU[p*nz+r];	
		else if (p >=0 && r <0)
		  
			FR[p*nz+r+nz] = pow(omega,(m_ord-l))*cpow((xi+I*eta),l)*FU[p*nz+r+nz];
	     
		else if (p <0 && r >=0)
		  
			FR[(p+nx)*nz+r]= pow(omega,(m_ord-l))*cpow((xi+I*eta),l)*FU[(p+nx)*nz+r];
		  
		else
		  
			FR[(p+nx)*nz+r+nz]= pow(omega,(m_ord-l))*cpow((xi+I*eta),l)*FU[(p+nx)*nz+r+nz];
       
	      }
	 
	    else
	      {
		if(m_ord ==0)
		FR[p*nz+r]= FU[p*nz+r];//the singularity
		else
		  FR[p*nz+r]=0;
	      }
	    
	  }
	}
      
      */
      if ( m_ord == 0)
	{//printf("entering m_ord =0 \n");
 for(p = 1;p<= (nx-1)/2;p++)
	{

	  xi = 2*pi*dxi*p;
	  for(r = 1; r <= (nz-1)/2; r++)
	    {
	    eta = 2*pi*deta*r;
	    //omega = pow((xi*xi+eta*eta),(0.5));
	    theta = atan(eta/xi);
	   
	    FR[p*nz+r] = (cos(l*theta)+I*sin(l*theta))*FU[p*nz+r];
	    
	    FR[p*nz-r+nz] = (cos(l*theta)-I*sin(l*theta))*FU[p*nz-r+nz];
	   
	    FR[(-p+nx)*nz+r]= (cos(l*(pi-theta))+I*sin(l*(pi-theta)))*FU[(-p+nx)*nz+r];
	    
	    FR[(-p+nx)*nz-r+nz]= (cos(l*(pi+theta))+I*sin(l*(pi+theta)))*FU[(-p+nx)*nz-r+nz];
	  }
	}
      for (r = 1;r<=(nz-1)/2;r++)
	{
	  // eta = 2*pi*deta*r;
	  
	  FR[r] =(cos(l*pi/2)+I*sin(l*pi/2))*FU[r];
	  FR[-r+nz] = (cos(l*pi/2)-I*sin(l*pi/2))*FU[-r+nz];
	}

      for (p=1;p<=(nx-1)/2;p++)
	{
	  //xi = 2*pi*dxi*p;
	  FR[p*nz]= (1)*FU[p*nz];
	  FR[(-p+nx)*nz]=(cos(l*pi))*FU[(-p+nx)*nz];
	}
      //printf("exiting m_ord = 0\n");
	}

      
    
      else
	{
	  //printf("entering m_ord !=0 \n");
      for(p = 1;p<= (nx-1)/2;p++)
	{

	  xi = 2*pi*dxi*p;
	  for(r = 1; r <= (nz-1)/2; r++)
	    {
	    eta = 2*pi*deta*r;
	    alpha = pow((xi*xi+eta*eta),(0.5*m_ord));
	    theta = atan(eta/xi);
	   
	    FR[p*nz+r] = alpha*(cos(l*theta)+I*sin(l*theta))*FU[p*nz+r];
	    
	    FR[p*nz-r+nz] = alpha*(cos(l*theta)-I*sin(l*theta))*FU[p*nz-r+nz];
	   
	    FR[(-p+nx)*nz+r]=alpha* (cos(l*(pi-theta))+I*sin(l*(pi-theta)))*FU[(-p+nx)*nz+r];
	    
	    FR[(-p+nx)*nz-r+nz]= alpha*(cos(l*(pi+theta))+I*sin(l*(pi+theta)))*FU[(-p+nx)*nz-r+nz];
	  }
	}
      for (r = 1;r<=(nz-1)/2;r++)
	{
	  eta = 2*pi*deta*r;
	  
	  FR[r] =pow(eta,(m_ord))*(cos(l*pi/2)+I*sin(l*pi/2))*FU[r];
	  FR[-r+nz] = pow(eta,(m_ord))*(cos(l*pi/2)-I*sin(l*pi/2))*FU[-r+nz];
	}
      for (p=1;p<=(nx-1)/2;p++)
	{
	  xi = 2*pi*dxi*p;
	  FR[p*nz]= pow(xi,m_ord)*(1)*FU[p*nz];
	  FR[(-p+nx)*nz]=pow(xi,m_ord)*(cos(l*pi))*FU[(-p+nx)*nz];
	}
      
      //printf("exiting  m_ord !=0\n");
	}
	
      
	if(m_ord ==0)
	  FR[0]=FU[0];
	else 
	  FR[0]=0;
      
      
 
      fftw_execute(plan_backward);

      //printf("fft inv ok ..\n");

      for (m=0;m<nx;m++)
	{
	  for(n=0;n<nz;n++)
	    {
	      QU[m*nz+n]=QU[m*nz+n]+ FQ[(k-1)/2+l+k*(m*nz+n)]*R[m*nz+n]/(nx*nz);
	     
	    }
	}

      //printf("setting QU in psido ok ..\n");

    }
  
  //printf("freeing memory...");
  //free memory 
  fftw_destroy_plan(plan_forward);
  fftw_destroy_plan(plan_backward);

  fftw_free(FU);
  fftw_free(FR);
  fftw_free(R);


  //printf("memory freed...");


}


