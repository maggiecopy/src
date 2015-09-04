#include "gridbandpassop.hh"

static void reverse (int n1, float* trace) {
    int i1;
    float t;

    for (i1=0; i1 < n1/2; i1++) { 
        t=trace[i1];
        trace[i1]=trace[n1-1-i1];
        trace[n1-1-i1]=t;
    }
}

butter butter_init(bool low     /* low-pass (or high-pass) */, 
		   float cutoff /* cut off frequency */, 
		   int nn       /* number of poles */)
/*< initialize >*/
{
    int j;
    float arg, ss, sinw, cosw, fact;
    butter bw;

    arg = 2.*M_PI*cutoff;
    sinw = sinf(arg);
    cosw = cosf(arg);

    bw = (butter) usermalloc_ (sizeof(butter));
    bw->nn = nn;
    bw->low = low;

    size_t i2;
    bw->den = (float**) usermalloc_ ((nn+1)/2*sizeof(float*));
    bw->den[0] = (float*) usermalloc_ (2*(nn+1)/2);
    for (i2=1; i2 < (nn+1)/2; i2++) {
	bw->den[i2] = bw->den[0]+i2*2;
    }

    if (nn%2) {
	if (low) {
	    fact = (1.+cosw)/sinw;
	    bw->den[nn/2][0] = 1./(1.+fact);
	    bw->den[nn/2][1] = 1.-fact;
	} else {
	    fact = sinw/(1.+cosw);
	    bw->den[nn/2][0] = 1./(fact+1.);
	    bw->den[nn/2][1] = fact-1.;
	}
    }

    fact = low? sinf(0.5*arg): cosf(0.5*arg);
    fact *= fact;
    
    for (j=0; j < nn/2; j++) {
	ss = sinf(M_PI*(2*j+1)/(2*nn))*sinw;
	bw->den[j][0] = fact/(1.+ss);
	bw->den[j][1] = (1.-ss)/fact;
    }
    bw->mid = -2.*cosw/fact;

    return bw;
}

void butter_apply (const butter bw, int nx, float *x /* data [nx] */)
/*< filter the data (in place) >*/
{
    int ix, j, nn;
    float d0, d1, d2, x0, x1, x2, y0, y1, y2;

    d1 = bw->mid;
    nn = bw->nn;

    if (nn%2) {
	d0 = bw->den[nn/2][0];
	d2 = bw->den[nn/2][1];
	x0 = y1 = 0.;
	for (ix=0; ix < nx; ix++) { 
	    x1 = x0; x0 = x[ix];
	    y0 = (bw->low)? 
		(x0 + x1 - d2 * y1)*d0:
		(x0 - x1 - d2 * y1)*d0;
	    x[ix] = y1 = y0;
	}
    }

    for (j=0; j < nn/2; j++) {
	d0 = bw->den[j][0];
	d2 = bw->den[j][1];
	x1 = x0 = y1 = y2 = 0.;
	for (ix=0; ix < nx; ix++) { 
	    x2 = x1; x1 = x0; x0 = x[ix];
	    y0 = (bw->low)? 
		(x0 + 2*x1 + x2 - d1 * y1 - d2 * y2)*d0:
		(x0 - 2*x1 + x2 - d1 * y1 - d2 * y2)*d0;
	    y2 = y1; x[ix] = y1 = y0;
	}
    }
}

void butter_close(butter bw)
/*< Free allocated storage >*/
{
    userfree_(bw->den[0]);
    userfree_(bw->den);
    userfree_(bw);
}


namespace TSOpt {

  using RVL::ScalarFieldTraits;
  using RVL::SpaceTest;
  using RVL::Operator;
  using RVL::LinearOp;
  using RVL::Space;
  using RVL::ProductSpace;
  using RVL::Vector;
  using RVL::Components;
  using RVL::ProtectedDivision;
  using RVL::RnArray;
  using RVL::RVLScale;
  using RVL::BinaryLocalFunctionObject;
  using RVL::RVLException;
  using RVL::ContentPackage;
  using RVL::LocalDataContainer;
  using RVL::MPISerialFunctionObject;
  using RVL::MPISerialFunctionObjectRedn;

  void GridScaleFO::operator()(LocalDataContainer<ireal> & x,
			      LocalDataContainer<ireal> const & y) {
    try {
      // cerr<<"GridWindowFO::operator() begin\n";
      ContentPackage< ireal, RARR > const & gy =
	dynamic_cast<ContentPackage< ireal, RARR > const &>(y);
            
      ContentPackage< ireal, RARR > & gx =
	dynamic_cast<ContentPackage< ireal, RARR > &>(x);
            
      //ContentPackage< ireal, RARR > const & gscale =
      //dynamic_cast<ContentPackage< ireal, RARR > const &>(vscale);
      // precondition - metadata are same dimn
      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();
      int dimx; int dimy;
      ra_ndim(&rax,&dimx);
      ra_ndim(&ray,&dimy);
      if (dimx != dimy) {
	RVLException e;
	e<<"Error: GridScaleFO::operator()\n";
	e<<"arguments have different dims:\n";
	e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
	throw e;
      }
      // compute grid params
      IPNT gsx; IPNT gex;
      IPNT gsy; IPNT gey;
      IPNT s; IPNT e;
      ra_a_gse(&rax,gsx,gex);
      ra_a_gse(&ray,gsy,gey);
      // calculate grid overlap
      for (int ii=0;ii<dimx;ii++)  {
	s[ii]=max(gsy[ii],gsx[ii]);
	e[ii]=min(gey[ii],gex[ii]);
      }
      IPNT i;
#if RARR_MAX_NDIM > 0
      if (dimx==1) {
#pragma ivdep
        for (i[0]=s[0];i[0]<=e[0];i[0]++) {
          rax._s1[i[0]]*=ray._s1[i[0]];
        }
      }
#endif
#if RARR_MAX_NDIM > 1
      if (dimx==2) {
	for (i[1]=s[1];i[1]<=e[1];i[1]++) {
#pragma ivdep
          for (i[0]=s[0];i[0]<=e[0];i[0]++) {
              rax._s2[i[1]][i[0]]*=ray._s2[i[1]][i[0]];
          }
	}
      }
#endif
#if RARR_MAX_NDIM > 2
      if (dimx==3) {
	for (i[2]=s[2];i[2]<=e[2];i[2]++) {
	  for (i[1]=s[1];i[1]<=e[1];i[1]++) {
#pragma ivdep
	    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
                rax._s3[i[2]][i[1]][i[0]]*=ray._s3[i[2]][i[1]][i[0]];
	    }
	  }
	}
      }
#endif
      if (dimx<1 || dimx>3) {
	RVLException e;
	e<<"Error: GridScaleFO::operator()\n";
	e<<"dim = "<<dimx<<" outside of admissible set {1, 2, 3}\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"\nError: GridScaleFO::operator()\n";
      e<<"at least one arg is not ContentPackage<ireal,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridScaleFO::operator()\n";
      throw e;
    }
  }

  void GridBandPassFO::operator()(LocalDataContainer<ireal> & x,
			      LocalDataContainer<ireal> const & y) {
    try {
      // cerr<<"GridWindowFO::operator() begin\n";
      ContentPackage< ireal, RARR > const & gy =
	dynamic_cast<ContentPackage< ireal, RARR > const &>(y);
            
      ContentPackage< ireal, RARR > & gx =
	dynamic_cast<ContentPackage< ireal, RARR > &>(x);
            
      //ContentPackage< ireal, RARR > const & gscale =
      //dynamic_cast<ContentPackage< ireal, RARR > const &>(vscale);
      // precondition - metadata are same dimn
      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();
      int dimx; int dimy;
      ra_ndim(&rax,&dimx);
      ra_ndim(&ray,&dimy);
      if (dimx != dimy) {
	RVLException e;
	e<<"Error: GridBandPassFO::operator()\n";
	e<<"arguments have different dims:\n";
	e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
	throw e;
      }
      // compute grid params
      IPNT gsx; IPNT gex;
      IPNT gsy; IPNT gey;
      IPNT s; IPNT e;
      ra_a_gse(&rax,gsx,gex);
      ra_a_gse(&ray,gsy,gey);
      // calculate grid overlap
      for (int ii=0;ii<dimx;ii++)  {
	s[ii]=max(gsy[ii],gsx[ii]);
	e[ii]=min(gey[ii],gex[ii]);
      }
      int nz=e[0]-s[0]+1;
      float * indata = (float *)usermalloc_(nz*sizeof(float));      
      IPNT i;
      RPNT fac;
      RASN(fac,RPNT_1);
#if RARR_MAX_NDIM > 0
      if (dimx==1) {
#pragma ivdep
        for (i[0]=s[0];i[0]<=e[0];i[0]++) {
          indata[i[0]-s[0]]=ray._s1[i[0]];
        }
        if (blo != NULL) {
          butter_apply (blo, nz, indata);

          if (!phase) {
            reverse (nz, indata);
            butter_apply (blo, nz, indata);
            reverse (nz, indata);
          }
        }

        if (bhi != NULL) {
          butter_apply (bhi, nz, indata);
            
          if (!phase) {
            reverse (nz, indata);
            butter_apply (bhi, nz, indata);    
            reverse (nz, indata);                   
          }
        }
#pragma ivdep
        for (i[0]=s[0];i[0]<=e[0];i[0]++) {
          rax._s1[i[0]]=indata[i[0]-s[0]];
        }
      }
#endif
#if RARR_MAX_NDIM > 1
      if (dimx==2) {
	for (i[1]=s[1];i[1]<=e[1];i[1]++) {
#pragma ivdep
          for (i[0]=s[0];i[0]<=e[0];i[0]++) {
              indata[i[0]-s[0]]=ray._s2[i[1]][i[0]];
          }
          if (blo != NULL) {
            butter_apply (blo, nz, indata);

            if (!phase) {
              reverse (nz, indata);
              butter_apply (blo, nz, indata);
              reverse (nz, indata);
            }
          }

          if (bhi != NULL) {
            butter_apply (bhi, nz, indata);
              
            if (!phase) {
              reverse (nz, indata);
              butter_apply (bhi, nz, indata);    
              reverse (nz, indata);                   
            }
          }
#pragma ivdep
          for (i[0]=s[0];i[0]<=e[0];i[0]++) {
              rax._s2[i[1]][i[0]]=indata[i[0]-s[0]];
          }
	}
      }
#endif
#if RARR_MAX_NDIM > 2
      if (dimx==3) {
	for (i[2]=s[2];i[2]<=e[2];i[2]++) {
	  for (i[1]=s[1];i[1]<=e[1];i[1]++) {
#pragma ivdep
	    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
                indata[i[0]-s[0]]=ray._s3[i[2]][i[1]][i[0]];
	    }
            if (blo != NULL) {
              butter_apply (blo, nz, indata);
            
              if (!phase) {
                reverse (nz, indata);
                butter_apply (blo, nz, indata);
                reverse (nz, indata);
              }
            }
            
            if (bhi != NULL) {
              butter_apply (bhi, nz, indata);
                
              if (!phase) {
                reverse (nz, indata);
                butter_apply (bhi, nz, indata); 
                reverse (nz, indata);		
              }
            }
#pragma ivdep
	    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
		rax._s3[i[2]][i[1]][i[0]]=indata[i[0]-s[0]];
	    }
	  }
	}
      }
#endif
      if (dimx<1 || dimx>3) {
	RVLException e;
	e<<"Error: GridBandPassFO::operator()\n";
	e<<"dim = "<<dimx<<" outside of admissible set {1, 2, 3}\n";
	throw e;
      }
      userfree_(indata);
    }
    catch (bad_cast) {
      RVLException e;
      e<<"\nError: GridBandPassFO::operator()\n";
      e<<"at least one arg is not ContentPackage<ireal,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridBandPassFO::operator()\n";
      throw e;
    }
  }
    
   
#ifdef IWAVE_USE_MPI
  typedef MPIGridSpace myGridSpace;
#else
  typedef GridSpace myGridSpace;
#endif
    
  GridBandPassOp::GridBandPassOp(
                  Space<ireal> const & _dom,
                  Vector<ireal> const & _vscale,
                  bool _appsc,
                  ireal _flo, ireal _fhi, bool _phase,
                  int _nplo, int _nphi)
    : dom(_dom), vscale(_vscale), flo(_flo), fhi(_fhi),
      phase(_phase), nplo(_nplo), nphi(_nphi), appsc(_appsc) {
    try {
  
      blo=NULL, bhi=NULL;
      const float eps=0.0001;
      IASN(n_arr,IPNT_0);
      RASN(d_arr,RPNT_0);
      // branch on product structure - unfortunately required
      ProductSpace<ireal> const * pdom = dynamic_cast<ProductSpace<ireal> const *>(&dom);
      ProductSpace<ireal> const * prng = dynamic_cast<ProductSpace<ireal> const *>(&(vscale.getSpace()));
      if (pdom && prng) {
        if (pdom->getSize() != prng->getSize()) {
          RVLException e;
          e<<"Error GridBandPassOp constructor\n";
          e<<"  domain and scale vector are product spaces with different numbers of components\n";
          throw e;
        }
	myGridSpace const & gref = dynamic_cast<myGridSpace const &>((*pdom)[0]);
	if (retrieveGlobalRank()==0) {
	  grid const & g = gref.getGrid();
          get_d(d_arr,g);
          get_n(n_arr,g);
          if (flo < 0.) {
            RVLException e;
            e<<"Error: GridBandPassOp::GridBandPassOp\n";
            e<<"  Negative input low freq flo=" << flo << "\n";
            throw e;
          }          
          flo = flo * d_arr[0];
          fhi = fhi * d_arr[0];
          if (fhi < flo) {
            RVLException e;
            e<<"Error: GridBandPassOp::GridBandPassOp\n";
            e<<"  Need flo < fhi, got flo=" << flo << ", fhi=" << fhi << "\n";
            throw e;
          } 
          if (fhi > 0.5f){
            fhi = 0.5f;
          }
          if (nplo < 1)           nplo=1;
          if (nplo > 1 && !phase) nplo = nplo/2;
          if (nphi < 1)           nphi = 1;
          if (nphi > 1 && !phase) nphi = nphi/2; 
          if (flo > eps)     blo = butter_init(false, flo, nplo);
          if (fhi < 0.5-eps) bhi = butter_init(true,  fhi, nphi);
	}
      }
      else {
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &> (dom);
	myGridSpace const & grng = dynamic_cast<myGridSpace const &>(vscale.getSpace());
	if (retrieveGlobalRank()==0) {
	  if (compatible_grid(gdom.getGrid(),grng.getGrid())) {
	    RVLException e;
	    e<<"Error: GridBandPassOp constructor\n";
	    e<<"  domain, scale vector defined on incompatible grids\n";
	    e<<"  domain:\n";
	    for (int i=0;i<gdom.getGrid().gdim;i++) 
	      e<<"    axis "<<i<<" d="<<gdom.getGrid().axes[i].d<<" o="<<gdom.getGrid().axes[i].o<<"\n";
	    e<<"  range:\n";
	    for (int i=0;i<grng.getGrid().gdim;i++) 
	      e<<"    axis "<<i<<" d="<<grng.getGrid().axes[i].d<<" o="<<grng.getGrid().axes[i].o<<"\n";
	    throw e;
	  }
	  grid const & g = gdom.getGrid();
          get_d(d_arr,g);
          get_n(n_arr,g);
          if (flo < 0.) {
            RVLException e;
            e<<"Error: GridBandPassOp::GridBandPassOp\n";
            e<<"  Negative input low freq flo=" << flo << "\n";
            throw e;
          }          
          flo = flo * d_arr[0];
          fhi = fhi * d_arr[0];
          if (fhi < flo) {
            RVLException e;
            e<<"Error: GridBandPassOp::GridBandPassOp\n";
            e<<"  Need flo < fhi, got flo=" << flo << ", fhi=" << fhi << "\n";
            throw e;
          } 
          if (fhi > 0.5f){
            fhi = 0.5f;
          }
          if (nplo < 1)           nplo=1;
          if (nplo > 1 && !phase) nplo = nplo/2;
          if (nphi < 1)           nphi = 1;
          if (nphi > 1 && !phase) nphi = nphi/2; 
 
          if (flo > eps)     blo = butter_init(false, flo, nplo);
          if (fhi < 0.5-eps) bhi = butter_init(true,  fhi, nphi);
	}
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridWindowOp constructor\n";
      e<<"  domain is neither product nor a GridSpace,\n";
      e<<"  or some component is not a GridSpace\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridMaskOp constructor\n";
      throw e;
    }
  }
    
  void GridBandPassOp::apply(Vector<ireal> const & x,
			 Vector<ireal> & y) const {
    try {
      GridBandPassFO op(blo,bhi,phase);
      MPISerialFunctionObject<ireal> mpiop(op);
      y.zero();
      y.eval(mpiop,x);
      if(appsc){
        GridScaleFO sop('m');
        MPISerialFunctionObject<ireal> mpisop(sop);
        y.eval(mpisop,vscale);
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridBandPassOp::apply\n";
      throw e;
    }
  }
    
  void GridBandPassOp::applyAdj(Vector<ireal> const & x,
			      Vector<ireal> & y) const {
    try {
      GridBandPassFO op(blo,bhi,phase);  
      MPISerialFunctionObject<ireal> mpiop(op);
      y.zero();
      y.eval(mpiop,x);
      if (appsc){
        GridScaleFO sop('m');
        MPISerialFunctionObject<ireal> mpisop(sop);
        y.eval(mpisop,vscale);
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridBandPassOp::applyDeriv\n";
      throw e;
    }
  }
    
  GridBandPassOp::~GridBandPassOp(){
    if (blo != NULL) butter_close(blo);
    if (bhi != NULL) butter_close(bhi);
  }

  ostream & GridBandPassOp::write(ostream & str) const {
    if (!retrieveGlobalRank()) {
      str<<"GridBandPassOp\n";
    }
    return str;
  }
    
}
