#ifndef __RVLGRID__STACKOP__
#define __RVLGRID__STACKOP__

#include "gridspace.H"
#include "gridstack.H"
#include "linop.H"
#include "op.H"

using namespace RVL;

namespace RVLGrid {

  vector<int> intseq(int n);

  class SpreadOp: public LinearOpPairFO<float> {

  private: 

    GridSpace<float> rng;
    Grid<float> domgrid;
    GridSpace<float> dom;
    Spread spr;
    Stack stk;

    SpreadOp();

  protected:

    LinearOpPair<float> * clone() const { return new SpreadOp(*this); }

  public:

    SpreadOp(const SpreadOp & sop) 
      : LinearOpPairFO<float>(*this), 
	rng(sop.rng),
	domgrid(sop.domgrid),
	dom(sop.dom), 
	spr(), 
	stk() {}

    SpreadOp(const GridSpace<float> & sp, int def) 
      : LinearOpPairFO<float>(dom,rng,spr,stk),
	rng(sp), 
	domgrid(rng.getGrid(),intseq(rng.getGrid().getNaxes()-def)),
	dom(domgrid), 
	spr(), 
	stk() {} 
    
    ~SpreadOp() {}

  };

  class NLSpreadOp: public LNLOperator<float> {

    friend class OperatorEvaluation<float>;

  private:

    Vector<float> cmax;
    Vector<float> cmin;
    const Vector<float> & cminref;
    const Vector<float> & cmaxref;

    NLSpreadOp();

  protected:

    Operator<float> * clone() const {
      return new NLSpreadOp(*this);
    }

  public:

    NLSpreadOp(SpreadOp & LL, float _cmin, float _cmax);
    NLSpreadOp(const NLSpreadOp & op);
    ~NLSpreadOp() {}
    virtual float getMaxStep(const Vector<float> & x,
			     const Vector<float> & dx) const;
  };
}

#endif


