#include "conv.h"

int main(int argc, char **argv) {

  complex *x,*y,*c;
  int nx=3;
  int ny=4;
  int nc;

  x = (complex *)malloc(nx*sizeof(complex));
  y = (complex *)malloc(ny*sizeof(complex));


  x[0]=0.5;x[1]=0;x[2]=0.5;
  y[0]=1;y[1]=2;y[2]=2+I*1;y[3]=3;

  nc = 2*nx-1;
  c = (complex *)malloc(nc*sizeof(complex));

  conv(x,x,c,nx,nx);
  printf("conv(x,x)\n");
  int m;
  for (m=0;m<nc;m++){
    printf("%f ",creal(c[m]));printf("+I* %f\n",cimag(c[m]));
  }
  
  nc = nx+ny-1;
  c = (complex *)malloc(nc*sizeof(complex));
  conv(x,y,c,nx,ny);
  printf("conv(x,y)\n");
  
  for (m=0;m<nc;m++){
    printf("%f ",creal(c[m]));printf("+I* %f\n",cimag(c[m]));
  }

  nc = 2*ny-1;
  c = (complex *)malloc(nc*sizeof(complex));
  conv(y,y,c,ny,ny);
  printf("conv(y,y) \n");
  
  for (m=0;m<nc;m++){
    printf("%f ",creal(c[m]));printf("+I* %f\n",cimag(c[m]));
  }
  
  free(x);free(y);free(c);

}
