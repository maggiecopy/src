#ifndef __RVLGRID_DERIV
#define __RVLGRID_DERIV

#include "gridspace.H"
#include "linop.H"

using namespace RVL;

namespace RVLGrid {

  /* finite difference approximation to indx partial deriv
     on regular grid LDC. Central difference except at ends.
     In each copy of the indx axis,

     out(0)   = (in(1)-in(0))/d;
     out(n-1) = (in(n-1)-in(n-2))/d
     out(i)   = (in(i+1)-in(i-1))/(2d), i=i,...,n-2

  */

  class DerivFO: public BinaryLocalFunctionObject<float> {

  private:

    int indx;

  public:

    DerivFO(int _indx=0): indx(_indx) {}
    DerivFO(const DerivFO & d): indx(d.indx) {}
    ~DerivFO() {}

    void operator()(LocalDataContainer<float> & y,    // output
		    LocalDataContainer<float> const & x);   // input

    bool readsData(int i) { if (i==1) return true; return false; }
    bool writesData(int i) { return !readsData(i); }

    string getName() const { string tmp = "DerivFO"; return tmp; }

  };

  /* finite difference approximation to adjoint of indx partial deriv
     on regular grid LDC. Central difference except at ends.
     In each copy of the indx axis,

     out(0)   = -in(1)/2d-in(0))/d;
     out(1)   =  in(0)/d - in(2)/2d
     out(n-2) =  in(n-3)/2d - in(n-1)/d
     out(n-1) = in(n-1)/2d+in(n-2)/d
     out(i)   = (-in(i+1)+in(i-1))/(2d), i=2,...,n-3

  */

  class AdjDerivFO: public BinaryLocalFunctionObject<float> {

  private:

    int indx;

  public:

    AdjDerivFO(int _indx=0): indx(_indx) {}
    AdjDerivFO(const AdjDerivFO & d): indx(d.indx) {}
    ~AdjDerivFO() {}

    void operator()(LocalDataContainer<float> & y,    // output
		    LocalDataContainer<float> const & x);   // input

    string getName() const { string tmp = "AdjDerivFO"; return tmp; }

  };

  class DerivOp: public LinearOpPair<float> {

  private:

    const GridSpace<float> dom;
    int indx;

    DerivOp();

  protected:

    virtual LinearOpPair<float> * clone() const { return new DerivOp(*this); }

    void apply(const Vector<float> & x,
	       Vector<float> & y) const;

    void applyAdj(const Vector<float> & x,
		  Vector<float> & y) const;

  public:

    DerivOp(Grid<float> & g, int _indx);
    DerivOp(const DerivOp & d): dom(d.dom), indx(d.indx) {}

    ~DerivOp() {}

    const Space<float> & getDomain() const { return dom; }
    const Space<float> & getRange() const { return dom; }

    void write(RVLException & e) const;
    ostream & write(ostream & str) const;

  };
}

#endif
