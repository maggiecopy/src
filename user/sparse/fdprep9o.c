/* 2D 9-point finite difference scheme */
/* C-H. Jo, C. Shin, J.H. Suh, 1996, An optimal 9-point, finite 
   difference, frequency-space, 2-D scalar wave extrapolator, Geophysics, 
   61, 529-537. */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
#include <umfpack.h>

#include "fdprep9o.h"

void fdprep9o(const double omega,
	      const float a0, const float f0,
	      const int n1, const int n2,
	      const float d1, const float d2,
	      float **v,
	      const int npml,
	      const int pad1, const int pad2,
	      SuiteSparse_long n, SuiteSparse_long nz,
	      SuiteSparse_long *Ti, SuiteSparse_long *Tj,
	      double* Tx, double *Tz)
/*< discretization >*/
{
    int i, j, index;
    double eta1, eta2, c;
    double *g1, *g2, **pad;
    double complex *s1, *s2, neib, cent;
    SuiteSparse_long count;
    
    /* prepare PML */
    eta1 = (double) npml*d1;
    eta2 = (double) npml*d2;
    
    c = 2.*SF_PI*a0*f0;

    g1 = (double*) sf_alloc(pad1,sizeof(double));
    for (i=0; i < pad1; i++) {
	if (i < npml) {
	    g1[i] = c*pow(((npml-i)*d1/eta1),2.);
	} else if (i >= pad1-npml) {
	    g1[i] = c*pow(((i-(pad1-npml-1))*d1/eta1),2.);
	} else {
	    g1[i] = 0.;
	}
    }

    s1 = (double complex*) sf_alloc(pad1,sizeof(double complex));
    for (i=0; i < pad1; i++) {
	s1[i] = 1.-I*g1[i]/omega;
    }

    g2 = (double*) sf_alloc(pad2,sizeof(double));
    for (j=0; j < pad2; j++) {
	if (j < npml) {
	    g2[j] = c*pow(((npml-j)*d2/eta2),2.);
	} else if (j >= pad2-npml) {
	    g2[j] = c*pow(((j-(pad2-npml-1))*d2/eta2),2.);
	} else {
	    g2[j] = 0.;
	}
    }

    s2 = (double complex*) sf_alloc(pad2,sizeof(double complex));
    for (j=0; j < pad2; j++) {
	s2[j] = 1.-I*g2[j]/omega;
    }

    /* extend model */
    pad = (double**) sf_alloc(pad2,sizeof(double*));
    pad[0] = (double*) sf_alloc(pad1*pad2,sizeof(double));
    for (j=1; j < pad2; j++) {
	pad[j] = pad[0]+j*pad1;
    }

    for (j=0; j < npml; j++) {
	for (i=0; i < npml; i++) {
	    pad[j][i] = v[0][0];
	}
	for (i=npml; i < pad1-npml; i++) {
	    pad[j][i] = v[0][i-npml];
	}
	for (i=pad1-npml; i < pad1; i++) {
	    pad[j][i] = v[0][n1-1];
	}
    }
    for (j=npml; j < pad2-npml; j++) {
	for (i=0; i < npml; i++) {
	    pad[j][i] = v[j-npml][0];
	}
	for (i=npml; i < pad1-npml; i++) {
	    pad[j][i] = v[j-npml][i-npml];
	}
	for (i=pad1-npml; i < pad1; i++) {
	    pad[j][i] = v[j-npml][n1-1];
	}
    }
    for (j=pad2-npml; j < pad2; j++) {
	for (i=0; i < npml; i++) {
	    pad[j][i] = v[n2-1][0];
	}
	for (i=npml; i < pad1-npml; i++) {
	    pad[j][i] = v[n2-1][i-npml];
	}
	for (i=pad1-npml; i < pad1; i++) {
	    pad[j][i] = v[n2-1][n1-1];
	}
    }

    /* assemble matrix in triplet form */
    count = 0;
    for (j=1; j < pad2-1; j++) {
	for (i=1; i < pad1-1; i++) {
	    index = (j-1)*(pad1-2)+(i-1);
	    
	    cent = 0.+I*0.;

	    /* left */
	    neib = 0.5461*(s2[j]/s1[i]+s2[j]/s1[i-1])/(2.*d1*d1);
	    cent += -neib;
	    neib += 0.09381*pow(omega/pad[j][i-1],2.)*(s1[i-1]*s2[j]);

	    if (i != 1) {
		Ti[count] = index;
		Tj[count] = index-1;
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* right */
	    neib = 0.5461*(s2[j]/s1[i]+s2[j]/s1[i+1])/(2.*d1*d1);
	    cent += -neib;
	    neib += 0.09381*pow(omega/pad[j][i+1],2.)*(s1[i+1]*s2[j]);

	    if (i != pad1-2) {
		Ti[count] = index;
		Tj[count] = index+1;
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* down */
	    neib = 0.5461*(s1[i]/s2[j]+s1[i]/s2[j-1])/(2.*d2*d2);
	    cent += -neib;
	    neib += 0.09381*pow(omega/pad[j-1][i],2.)*(s1[i]*s2[j-1]);

	    if (j != 1) {
		Ti[count] = index;
		Tj[count] = index-(pad1-2);
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* up */
	    neib = 0.5461*(s1[i]/s2[j]+s1[i]/s2[j+1])/(2.*d2*d2);
	    cent += -neib;
	    neib += 0.09381*pow(omega/pad[j+1][i],2.)*(s1[i]*s2[j+1]);

	    if (j != pad2-2) {
		Ti[count] = index;
		Tj[count] = index+(pad1-2);
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* left down */
	    neib = 0.4539*(s2[j]/s1[i]+s2[j-1]/s1[i-1])/(2.*(d1*d1+d2*d2));
	    cent += -neib;
	    neib += -0.00001*pow(omega/pad[j-1][i-1],2.)*(s1[i-1]*s2[j-1]);

	    if (i != 1 && j != 1) {
		Ti[count] = index;
		Tj[count] = index-1-(pad1-2);
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* right up */
	    neib = 0.4539*(s2[j]/s1[i]+s2[j+1]/s1[i+1])/(2.*(d1*d1+d2*d2));
	    cent += -neib;
	    neib += -0.00001*pow(omega/pad[j+1][i+1],2.)*(s1[i+1]*s2[j+1]);

	    if (i != pad1-2 && j != pad2-2) {
		Ti[count] = index;
		Tj[count] = index+1+(pad1-2);
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }
	    
	    /* left up */
	    neib = 0.4539*(s1[i]/s2[j]+s1[i-1]/s2[j+1])/(2.*(d1*d1+d2*d2));
	    cent += -neib;
	    neib += -0.00001*pow(omega/pad[j+1][i-1],2.)*(s1[i-1]*s2[j+1]);

	    if (i != 1 && j != pad2-2) {
		Ti[count] = index;
		Tj[count] = index-1+(pad1-2);
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* right down */
	    neib = 0.4539*(s1[i]/s2[j]+s1[i+1]/s2[j-1])/(2.*(d1*d1+d2*d2));
	    cent += -neib;
	    neib += -0.00001*pow(omega/pad[j-1][i+1],2.)*(s1[i+1]*s2[j-1]);

	    if (i != pad1-2 && j != 1) {
		Ti[count] = index;
		Tj[count] = index+1-(pad1-2);
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }
	    
	    /* center */
	    cent += 0.6248*pow(omega/pad[j][i],2.)*(s1[i]*s2[j]);
	    
	    Ti[count] = index;
	    Tj[count] = index;
	    Tx[count] = creal(cent);
	    Tz[count] = cimag(cent);
	    
	    count++;
	}
    }    
}

void fdpad9o(const int npml,
	     const int pad1, const int pad2,
	     sf_complex **dat,
	     double *Bx, double *Bz)
/*< pad >*/
{
    int i, j;

    for (j=1; j < pad2-1; j++) {
	for (i=1; i < pad1-1; i++) {
	    if (i < npml || i >= pad1-npml || 
		j < npml || j >= pad2-npml) {
		Bx[(j-1)*(pad1-2)+(i-1)] = 0.;
		Bz[(j-1)*(pad1-2)+(i-1)] = 0.;
	    } else {
		Bx[(j-1)*(pad1-2)+(i-1)] = creal(dat[j-npml][i-npml]);
		Bz[(j-1)*(pad1-2)+(i-1)] = cimag(dat[j-npml][i-npml]);
	    }
	}
    }
}

void fdcut9o(const int npml,
	     const int pad1, const int pad2,
	     sf_complex **dat,
	     double *Xx, double *Xz)
/*< cut >*/
{
    int i, j;

    for (j=npml; j < pad2-npml; j++) {
	for (i=npml; i < pad1-npml; i++) {
	    dat[j-npml][i-npml] = sf_cmplx((float) Xx[(j-1)*(pad1-2)+(i-1)], 
					   (float) Xz[(j-1)*(pad1-2)+(i-1)]);
	}	    
    }
}