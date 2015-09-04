#ifndef __RVL_GRID_FLOAT_CUBIC
#define __RVL_GRID_FLOAT_CUBIC

#include "griddata.hh"

namespace RVLGrid {

  using RVL::BinaryLocalFunctionObject;
  using RVL::RVLException;

  /** CubicAdjInterp and CubicFwdInterp have the SAME order of
      arguments in their constructor: ie. for interpolation
      from R^{n1} to R^{n2}, BOTH constructors look like
      (classname)(n1,...,n2,...) 
  */
  class CubicFwdInterp: public BinaryLocalFunctionObject<float> {
  private:
    // Data for the deprecated constructor
    int nin, nout;
    float oin, oout, din, dout;
    int wl;
    float * work;
  protected:
    
    void cubic_interp(float & oin, float & din, float const * datain, int & nin,
		      float & oout, float & dout, float * dataout, int & nout, 
		      int & iend, float * work, int & wl, int & ier);

  public:

    CubicFwdInterp(): nin(0), nout(0), 
		      oin(0.0), oout(0.0), 
		      din(0.0), dout(0.0),
		      wl(0), work(NULL) {}
    // Deprecated constructor.
    CubicFwdInterp(int _nin, float _din, float _oin,
		   int _nout, float _dout, float _oout)
      :nin(_nin), nout(_nout), oin(_oin), oout(_oout), 
       din(_din), dout(_dout) 
    {
      wl = 4*nin;
      work = new float[wl];
    }

    CubicFwdInterp(const CubicFwdInterp & f)
      : nin(f.nin), nout(f.nout), 
	oin(f.oin), oout(f.oout), 
	din(f.din), dout(f.dout),
	wl(f.wl), work(NULL) {
      if( wl > 0 )
	work = new float[wl];
    }
    ~CubicFwdInterp() {
      if( work != NULL ) delete [] work;
    }
  
    void operator() (LocalDataContainer<float> & dcout,
		     LocalDataContainer<float> const & dcin);
    string getName() const { return "CubicFwdInterp"; }
  };

  /** CubicAdjInterp and CubicFwdInterp have the SAME order of
      arguments in their constructor: ie. for interpolation
      from R^{n1} to R^{n2}, BOTH constructors look like
      (classname)(n1,...,n2,...) 
  */
  class CubicAdjInterp: public BinaryLocalFunctionObject<float> {
  private:
    // Data for the deprecated constructor
    int nin, nout;
    float oin, oout, din, dout;
    int wl;
    float * work;
  protected:

    void cubicadj_interp(float & oin, float & din, float const * datain, int & nin,
			 float & oout, float & dout, float * dataout, int & nout, 
			 int & iend, float * work, int & wl, int & ier);

  public:
    CubicAdjInterp(): nin(0), nout(0), 
		      oin(0.0), oout(0.0), 
		      din(0.0), dout(0.0),
		      wl(0), work(NULL) {}
    // Deprecated constructor.
    CubicAdjInterp(int _nin, float _din, float _oin,
		   int _nout, float _dout, float _oout)
      :nin(_nin), nout(_nout), oin(_oin), oout(_oout), 
       din(_din), dout(_dout) 
    {
      wl = 8*nin+8*nout;
      work = new float[wl];
    }
    CubicAdjInterp(const CubicAdjInterp & f)
      : nin(f.nin), nout(f.nout), 
	oin(f.oin), oout(f.oout), 
	din(f.din), dout(f.dout),
	wl(f.wl), work(NULL) {
      if( wl > 0 )
	work = new float[wl];
    }
    ~CubicAdjInterp() {
      if( work != NULL ) delete [] work;
    }
      
    void operator() (LocalDataContainer<float> & dcout,
		     LocalDataContainer<float> const & dcin);
    string getName() const { return "CubicAdjInterp"; }
  };
  
  class BiCubicFwdInterp: public BinaryLocalFunctionObject<float> {

  protected:

    void bicubic_interp(float & o0in, float & d0in, int & n0in,
			float & o1in, float & d1in, int & n1in, 
			float const * datain,
			float & o0out, float & d0out, int & n0out, 
			float & o1out, float & d1out, int & n1out, 
			float * dataout,
			int & iend, float * work,
			int & wl, int & ier);

  public:

    BiCubicFwdInterp() {}
    BiCubicFwdInterp(const BiCubicFwdInterp &) {}
    ~BiCubicFwdInterp() {}
  
    void operator() (LocalDataContainer<float> & dout,
		     LocalDataContainer<float> const & din);
    string getName() const { return "BiCubicFwdInterp"; }
  };
    
  class BiCubicAdjInterp: public BinaryLocalFunctionObject<float> {
  private:
    int rn;   // if set, range treated as Rn for adjoint defn
  protected:
    void bicubicadj_interp(float & o0in, float & d0in, int & n0in,
			   float & o1in, float & d1in, int & n1in, 
			   float const * datain,
			   float & o0out, float & d0out, int & n0out, 
			   float & o1out, float & d1out, int & n1out, 
			   float * dataout,
			   int & iend, float * work,
			   int & wl, int & ier);
    
  public:

    BiCubicAdjInterp(int _rn=0): rn(_rn) {}
    BiCubicAdjInterp(const BiCubicAdjInterp & b): rn(b.rn) {}
    ~BiCubicAdjInterp() {}
  
    void operator() (LocalDataContainer<float> & dout,
		     LocalDataContainer<float> const & din);
    string getName() const { return "BiCubicAdjInterp"; }
  };

  class TriCubicFwdInterp: public BinaryLocalFunctionObject<float> {
  protected:
    void cubic_interp(float & oin, float & din, float const * datain, int & nin,
		      float & oout, float & dout, float * dataout, int & nout, 
		      int & iend, float * work, int & wl, int & ier);
    
    void bicubic_interp(float & o0in, float & d0in, int & n0in,
			float & o1in, float & d1in, int & n1in, 
			float const * datain,
			float & o0out, float & d0out, int & n0out, 
			float & o1out, float & d1out, int & n1out, 
			float * dataout,
			int & iend, float * work,
			int & wl, int & ier);
  public:

    TriCubicFwdInterp() {}
    TriCubicFwdInterp(const TriCubicFwdInterp &) {}
    ~TriCubicFwdInterp() {}
  
    void operator() (LocalDataContainer<float> & dout,
		     LocalDataContainer<float> const & din);
    string getName() const { return "TriCubicFwdInterp"; }
  };
      
  class TriCubicAdjInterp: public BinaryLocalFunctionObject<float> {
  protected:
    void cubicadj_interp(float & oin, float & din, float const * datain, int & nin,
			 float & oout, float & dout, float * dataout, int & nout, 
			 int & iend, float * work, int & wl, int & ier);

    void bicubicadj_interp(float & o0in, float & d0in, int & n0in,
			   float & o1in, float & d1in, int & n1in, 
			   float const * datain,
			   float & o0out, float & d0out, int & n0out, 
			   float & o1out, float & d1out, int & n1out, 
			   float * dataout,
			   int & iend, float * work,
			   int & wl, int & ier);
  public:

    TriCubicAdjInterp() {}
    TriCubicAdjInterp(const TriCubicAdjInterp &) {}
    ~TriCubicAdjInterp() {}
  
    void operator() (LocalDataContainer<float> & dout,
		     LocalDataContainer<float> const & din);
    string getName() const { return "TriCubicAdjInterp"; }
  };

} 
      
#endif








  
