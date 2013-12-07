/* Local prediction filter for complex numbers (n-dimensional). */
/*
  Copyright (C) 2006 University of Texas at Austin
   
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
   
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
   
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>

#include "cmultidivn.h"
#include "cweight2.h"

int main(int argc, char* argv[])
{
    bool verb;
    int m[1], rect[1];
    int i1,iw,niter,nshift,n1,n2,n12,ns,i3,n3,i,ii1,is;
    sf_complex *d1,**d,*sh, *f;
    float mean,*mk;
    sf_file dat, flt, mask, pre, zshift, zdata;

    sf_init(argc,argv);

    dat = sf_input("in");
    mask = sf_input("mask");
    flt = sf_output("out");
    if (!sf_getint("rect",rect)) rect[0]=10;
	/* smoothing radius (in space X) */

    if (!sf_getint("ns",&ns)) ns=1;
	/* shifts of both sides npef=2*ns+1*/

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */
    if (SF_COMPLEX != sf_gettype(dat)) sf_error("Need complex input");


    /*Predicted Data*/
    if (NULL != sf_getstring("pred")) {
	pre = sf_output("pred");
    } else {
	pre = NULL;
    }

    /*Shifts zeroed */
    if (NULL != sf_getstring("zshift")) {
	zshift = sf_output("zshift");
        sf_shiftdim(dat, zshift, 2);
        sf_putint(zshift,"n2",2*ns);
        sf_putfloat(zshift,"d2",1);
        sf_putfloat(zshift,"o2",-ns);
    } else {
	zshift = NULL;
    }

    /*Data zeroed */
    if (NULL != sf_getstring("zdata")) {
	zdata = sf_output("zdata");
    } else {
	zdata = NULL;
    }
    
    if (!sf_histint(dat,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(dat,"n2",&n2)) sf_error("No n2= in input");
    n12=n1*n2;
    n3 = sf_leftsize(dat,2);



    nshift=n1*2*ns;
    
    sf_shiftdim(dat, flt, 2);
    sf_putint(flt,"n2",2*ns);
    sf_putfloat(flt,"d2",1);
    sf_putfloat(flt,"o2",-ns);

    d  = sf_complexalloc2(n1,n2);
    d1 = sf_complexalloc(n1);
    sh = sf_complexalloc(nshift);
    mk = sf_floatalloc(n1);
    f  = sf_complexalloc(nshift);
    //f = sf_complexalloc(n12);
    //g = sf_complexalloc(nd);
    m[0]=n1;

    for (i3=0; i3<n3; i3++) { // y

         sf_complexread(d[0],n12,dat);
         sf_floatread(mk,n1,mask);

         for (iw=0; iw<n2; iw++) { // w
	      sf_warning("iw = %d of %d", iw, n2);
              for (i1=0;i1<n1;i1++){
                   d1[i1] = d[iw][i1];
                   for (is=0; is < (2*ns); is++) {
                        sh[is*n1+i1]=sf_cmplx(0,0);
                   }
              }
              // shifts
              for (is=0;is<ns;is++) {//forward
                   for (i1=0;i1<ns-is;i1++) {
                        sh[is*n1+i1]=sf_cmplx(0,0);
                   }
                   for (i1=ns-is;i1<n1;i1++) {
                        sh[is*n1+i1]=d1[i1-ns+is];
                   }
                   
              }
              for (is=0;is<ns;is++) { // backward
                   for (i1=0;i1<n1-is-1;i1++) {
                        sh[n1*ns+is*n1+i1]=d1[i1+is+1];
                   }
                   for (i1=n1-is-1;i1<n1;i1++) {
                        sh[n1*ns+is*n1+i1]=sf_cmplx(0,0);
                   }
                   
              }                   
              // shift end
              //Shifts end
              for (i1=0;i1<n1;i1++) { //Zero mask
                   if (mk[i1] == 0.) { 
                       for (ii1=i1-ns;ii1<i1+ns;ii1++){
                            if(ii1 >= 0 && ii1< n1){
                                 d1[ii1] = sf_cmplx(0,0);
                                 for (is=0; is < (2*ns); is++) {
                                      sh[is*n1+ii1]=sf_cmplx(0,0);
                                 }
                            }
                       }
                   } 
                   
              }         
              if (zshift) {
                  sf_complexwrite(sh,nshift,zshift);
              }
              if (zdata) {
                  sf_complexwrite(d1,n1,zdata);
              }
              
             cmultidivn_init(2*ns, 1, n1, m, rect, sh, false);    
             mean = 0.;
             for(i = 0; i < nshift; i++) {
#ifdef SF_HAS_COMPLEX_H
	         mean += crealf(conjf(sh[i])*sh[i]);
#else
	         mean += crealf(sf_cmul(conjf(sh[i]),sh[i]));
#endif
             }
             if (mean == 0.) {
                 for(i=0; i < nshift; i++) {
                     f[i]=sf_cmplx(0,0); ///////////////////////////////////////////////all zeros
                 }
	         sf_complexwrite(f,nshift,flt);  
	         continue;
             }

             mean = sqrtf (nshift/mean);
    
            for(i=0; i < nshift; i++) {
#ifdef SF_HAS_COMPLEX_H
	       sh[i] *= mean;
#else
 	       sh[i] = sf_crmul(sh[i],mean);
#endif
            }
            for(i=0; i < n1; i++) {
#ifdef SF_HAS_COMPLEX_H
	        d1[i] *= mean;
#else
	        d1[i] = sf_crmul(d1[i],mean);
#endif
            }
            cmultidivn (d1,f,niter);  
            sf_complexwrite(f,nshift,flt);

            if (pre) {
	        for(i=0; i < nshift; i++) {
#ifdef SF_HAS_COMPLEX_H
	           sh[i] /= mean;
#else
	           sh[i] = sf_crmul(sh[i],1./mean);
#endif
	        }
	
	       cweight2_lop(false,false,nshift,n1,f,d1);
	       sf_complexwrite(d1,n1,pre);
            }

         } // w end
    } // y end



    exit(0);
}