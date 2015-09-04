#ifndef __PSIDO_OP
#define __PSIDO_OP

#define DEFAULT_SNAPS 10

#include "parser.h"
#include "op.hh"
#include "linop.hh"
#include "ocdc.hh"
#include "gridpp.hh"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#endif
#include "traceio.h"

/////////  for PsiDO 
#include <fftw3.h>
#include "lbfgs.h"
#include "objective_real.h"    

typedef struct {
    fftwf_complex  *U;
    fftwf_complex  *QUtest;
    float *coeff_c;

    float * Obj;
    float *Grad;

    IPNT   nc;       // coarse grid info
    IPNT   n;        // image grid dim
    IPNT   d;        // image grid size
    IPNT   o;        // image grid origin

    float m_ord;     // order of psido
    int k;           // 2*k-1 #of angles
    float lam;       // L2 regularization weight
    FILE * stream;
} PsiDO_PARS;


namespace TSOpt {

  //using namespace RVLAlg;
  using RVL::DataContainer;
  using RVL::ProductDataContainer;
  using RVL::StdProductDataContainer;
  using RVL::Space;
  using RVL::parse;
  using RVL::valparse;
  using RVL::SpaceDCF;
  using RVL::ProductSpace;
  using RVL::ConstContainer;
  using RVL::STRING_PAIR;
  using RVL::Vector;
  using RVL::Components;
  using RVL::FunctionObject;
  using RVL::LinearOp;
  using RVL::Writeable;
  using RVL::AssignParams;
  
  class PsiDiffOp: public LinearOp<ireal>  {
      
  private:

    Space<ireal> const & dom;
    Vector<ireal> & coeff;
    int flag;                   // flag=0 use supplied coeff, otherwise compute  
    mutable FILE * stream;              /* output stream            */
    PARARRAY * pars;            /* parameter array ref copy */

    // verbosity control
    int dump_steps;
    int dump_pars;
    int dump_term;
    grid g;
    grid gc;
    string datatype;
    string cdatatype;
    // other verbosity control handled within iwave code
    
    // verbose announcements
    ostream & announce;

    PsiDiffOp();
      
  protected:
      
    void apply(const Vector<ireal> & x, 
	       Vector<ireal> & y) const;

    void applyAdj(const Vector<ireal> & x, 
		  Vector<ireal> & y) const;
    
  public:
    
    PsiDiffOp(Space<ireal> const & _dom,
              Vector<ireal> & _coeff,
              PARARRAY _pars, FILE * _stream,
              int _flag=0, 
              string _datatype="nontype",
              string _cdatatype="coeff",
	      ostream & _announce=cerr);
    
    PsiDiffOp(PsiDiffOp const & A);

    ~PsiDiffOp(){ps_delete(&pars);}

    LinearOp<float> * clone() const { return new PsiDiffOp(*this); }

    const Space<ireal> & getDomain() const { return dom; }
    const Space<ireal> & getRange() const { return dom; }

    // added 23.06.10 to facilitate using source as variable
    // without admitting that it's part of domain
    PARARRAY & getPar() { return *pars;}
    PARARRAY const & getPar() const { return *pars;}

    ostream & write(ostream & str) const {
      str<<"PsiDiffOp object\n";
      return str;
    }
  };
}

#endif
