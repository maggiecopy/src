#ifndef __TSOPT_GRIDBP_OPS__
#define __TSOPT_GRIDBP_OPS__

#include "rn.hh"
#include "op.hh"
#include "productspace.hh"
#include "mpiserialfo.hh"

#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif
// required for old NCAR-fft-based helmholtz
//#include <f2c.h>

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
struct Butter {
    bool low;
    int nn;
    float **den, mid;
};

typedef struct Butter *butter;


namespace TSOpt {

  using RVL::BinaryLocalFunctionObject;
  using RVL::RVLException;
  using RVL::ContentPackage;
  using RVL::LocalDataContainer;

    
  /** function object for BandPass a 1D grid function, in the form of
      a CP with Grid metadata. Scalar arrays ffset between input and
      output grids computed by calling unit
  */
  class GridScaleFO: public BinaryLocalFunctionObject<ireal> {
        
  private:
        
    char method;

    GridScaleFO();
        
  public:
        
    GridScaleFO(char _method)
    : method(_method){
    }
        
    GridScaleFO(GridScaleFO const & f)
    : method(f.method){
    }
        
    using RVL::LocalEvaluation<ireal>::operator();
    void operator()(LocalDataContainer<ireal> &,
		    LocalDataContainer<ireal> const &);
        
    string getName() const { string tmp = "GridScaleFO"; return tmp; }
        
  };
    
  class GridBandPassFO: public BinaryLocalFunctionObject<ireal> {
        
  private:
        
    butter blo, bhi;
    bool  phase;

    GridBandPassFO();
        
  public:
        
    GridBandPassFO(butter _blo, butter _bhi, bool _phase=false)
      : blo(_blo), bhi(_bhi), phase(_phase){
    }
        
    GridBandPassFO(GridBandPassFO const & f)
    : blo(f.blo), bhi(f.bhi),phase(f.phase){
    }
        
    using RVL::LocalEvaluation<ireal>::operator();
    void operator()(LocalDataContainer<ireal> &,
		    LocalDataContainer<ireal> const &);
        
    string getName() const { string tmp = "GridBandPassFO"; return tmp; }
        
  };
    
  /** BandPass operator for grid objects. Apply method outputs
      masked version of background Vector data member: thus
     
      \f$ y = x outside of mask, or 0 inside mask\f$
     
      Derivative and adjoint derivative are implement standard linear injection and
      projection operators.
  */
  class GridBandPassOp: public LinearOp<ireal> {
        
  private:
        
    Space<ireal> const & dom;
    Vector<ireal> const & vscale;
    ireal flo;
    ireal fhi;
    bool  phase;
    bool  appsc;
    int   nplo;
    int   nphi;

    butter blo, bhi;
    IPNT n_arr;
    RPNT d_arr;

    GridBandPassOp();
        
  protected:
        
    void apply(Vector<ireal> const &,
	       Vector<ireal> &) const;
    void applyAdj(Vector<ireal> const &,
		  Vector<ireal> &) const;
        
    Operator<ireal> * clone() const { return new GridBandPassOp(*this); }
        
  public:
        
    /** main constructor: takes Grid defining increment window, and width
	parameter defining zone of smooth taper - same on all sides.
    */
    GridBandPassOp(Space<ireal> const & _dom,
                   Vector<ireal> const & _vscale,
                   bool appsc=false,
	           ireal _flo=0.f, ireal _fhi=0.5f, bool phase=false, 
                   int _nplo=6, int _nphi=6);
        
    /** Copy constructor - memberwise */
    GridBandPassOp(GridBandPassOp const & op)
    : dom(op.dom), vscale(op.vscale), flo(op.flo), fhi(op.fhi),
      phase(op.phase), nplo(op.nplo), nphi(op.nphi), appsc(op.appsc){
    }
        
    ~GridBandPassOp();
        
    Space<ireal> const & getDomain() const { return dom; }
    Space<ireal> const & getRange() const { return dom; }
        
    ostream & write(ostream & str) const;
  };

    

    
}
#endif
