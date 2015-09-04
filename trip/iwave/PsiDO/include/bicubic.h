#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <complex.h>
#include "fftw3.h"

extern  void bicubic_(float & o0in, float & d0in, int & n0in,
			 float & o1in, float & d1in, int & n1in, 
			 float * datain,
			 float & o0out, float & d0out, int & n0out, 
			 float & o1out, float & d1out, int & n1out, 
			 float * dataout,
			 int & iend, float * work,
			 int & wl, int & ier);

extern  void bicubicadj_(float & o0in, float & d0in, int & n0in,
			    float & o1in, float & d1in, int & n1in, 
			    float * datain,
			    float & o0out, float & d0out, int & n0out, 
			    float & o1out, float & d1out, int & n1out, 
			    float * dataout,
			    int & iend, float * work,
			    int & wl, int & ier);
