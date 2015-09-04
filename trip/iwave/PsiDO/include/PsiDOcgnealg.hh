// PsiDOcgnealg.H
// created by Yin Huang 07/13/15

// Much code is shamelessly stolen from the umin cgne.H
// originally written by WWS, with his permission

/*************************************************************************

Copyright Rice University, 2004.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/


#ifndef __PsiDO_CGNE_H
#define __PsiDO_CGNE_H

#include "alg.hh"
#include "terminator.hh"
#include "linop.hh"
#include "table.hh"

using namespace RVLAlg;

namespace RVLUmin {

  using namespace RVL;
  using namespace RVLAlg;    

  // forward declaration
  template<typename Scalar>
  class PsiCGNEAlg;

  /** Preconditioned conjugate gradient iteration for the normal
      equations with PsiDO as a preconditioner.
      
      On construction, internal workspace allocated and initialized.
      Each step updates internal state of CGNEStep object. Since
      solution vector, residual norm, and normal residual norm are
      stored as mutable references to external objects, these external
      objects are updated as well.

      IMPORTANT NOTE: this version of the algorithm assumes that the
      solution vector reference (internal data member x) refers to a
      zero vector on initialization. To accommodate nontrivial initial
      guess, <i>modify the right-hand-side vector</i> (argument _b)
      externally.

      Solution vector (x), iteration count, residual norm, and
      gradient norm are all references to external objects, which may
      be monitored by appropriate terminators to build a LoopAlg out
      of this Algorithm.

      See CGNEAlg for description of a fully functional algorithm
      class, combining this step with a Terminator to make a LoopAlg.
  */
  template<typename Scalar>
  class PsiCGNEStep : public Algorithm {

    friend class PsiCGNEAlg<Scalar>;

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:

    PsiCGNEStep(LinearOp<Scalar> const & _A,
	      LinearOp<Scalar> const & _M,
	      Vector<Scalar> & _x,
	      Vector<Scalar> const & _b,
	      atype & _rnorm, 
	      atype & _nrnorm)
      : A(_A), M(_M), x(_x), b(_b), rnorm(_rnorm), nrnorm(_nrnorm), 
	r(A.getRange()), ng(A.getDomain()), g(A.getDomain()),
	q(A.getRange()), p(A.getDomain()),
        Aimg(A.getRange()), preng(A.getDomain()) { 
      // sanity tests - cannot sensibly test M for SPD, but 
      // at least get domain right
      if ((M.getDomain() != A.getDomain()) ||
	  (M.getRange()  != A.getDomain())) {
	RVLException e;
	e<<"Error PsiCGNEStep constructor\n";
	e<<"  preconditioning operator does not have domain, range\n";
	e<<"  same as domain of system operator\n";
	e<<"  preconditioning operator:\n";
	M.write(e);
	e<<"  system operator:\n";
	A.write(e);
	throw e;
      }
      // NOTE: initial x assumed to be zero vector
      r.copy(b);                   // r0=b-Ax0 with x0=0
      rnorm=r.norm();          
      A.applyAdjOp(r,ng);          // ng = A^T r0
      //modify according to PsiDO 
      //M.applyOp(ng,g);             // g  = M A^T r0
      A.applyOp(ng,Aimg);     
      A.applyAdjOp(Aimg,g);
      M.applyOp(ng,g);             // g  = M A^T r0
      p.copy(g);                   // p  = g
      gamma = abs(g.inner(ng));
      nrnorm=sqrt(gamma);
    }
      
    /**
       Run a single step of the conjugate gradients for the normal equations
    */
    void run() {
      try {
	A.applyOp(p,q);
	atype qtq = q.normsq();
	atype absalpha;
	if (ProtectedDivision<atype>(gamma,qtq,absalpha)) {
	  RVLException e;
	  e<<"Error: PsiCGNEStep::run() from ProtectedDivision: alpha\n";
	  throw e;
	}
	// can use q for workspace since it is not needed again until
	// reinitialized

	Scalar alpha=absalpha;
	x.linComb(alpha,p);
	r.linComb(-alpha,q);

        preng.copy(ng);
	A.applyAdjOp(r,ng);
        A.applyOp(ng,Aimg);
        A.applyAdjOp(Aimg,g);
	M.applyOp(ng,g);

	atype newgamma = abs(g.inner(ng));
        atype tmp = newgamma - abs(g.inner(preng));        
	atype absbeta;
	if (ProtectedDivision<atype>(tmp,gamma,absbeta)) {
	  RVLException e;
	  e<<"Error: PsiCGNEStep::run() from ProtectedDivision: beta\n";
	  throw e;
	}

	Scalar beta = absbeta;
	p.linComb(ScalarFieldTraits<Scalar>::One(),g,beta);
	gamma=newgamma;
	rnorm=r.norm();
	nrnorm=sqrt(gamma);
      }
      catch (RVLException & e) {
	e<<"\ncalled from PsiCGNEStep::run()\n";
	throw e;
      }
     
    }

    ~PsiCGNEStep() {}

  protected:

    // references to external objects
    LinearOp<Scalar> const & A;
    LinearOp<Scalar> const & M;
    Vector<Scalar> & x;
    Vector<Scalar> const & b;
    atype & rnorm;
    atype & nrnorm;

  private:

    // need four work vectors and one scalar as persistent object data
    Vector<Scalar> r;    // residual
    Vector<Scalar> ng;   // normal residual
    Vector<Scalar> g;    // gradient = preconditioned normal residual
    Vector<Scalar> q;    // image of search direction
    Vector<Scalar> p;    // search direction
    atype gamma;         // gradient norm, inner product defined by inverse of M
    // working vectors for PsiDO
    Vector<Scalar> Aimg;
    Vector<Scalar> preng;
  };
  

  template<typename Scalar>
  class PsiScaleCGNEStep : public Algorithm {

    friend class PsiCGNEAlg<Scalar>;

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:

    PsiScaleCGNEStep(LinearOp<Scalar> const & _A,
	      LinearOp<Scalar> const & _M,
	      LinearOp<Scalar> const & _S,
	      LinearOp<Scalar> const & _S2,
	      Vector<Scalar> & _x,
	      Vector<Scalar> const & _b,
	      atype & _rnorm, 
	      atype & _nrnorm)
      : A(_A), M(_M), S(_S), S2(_S2), x(_x), b(_b), rnorm(_rnorm), nrnorm(_nrnorm), 
	r(A.getRange()), ng(A.getDomain()), g(A.getDomain()),
        sng(A.getDomain()), sg(A.getDomain()),
	q(A.getRange()), p(A.getDomain()), swg(A.getDomain()),
        Aimg(A.getRange()), preng(A.getDomain()) { 
      // sanity tests - cannot sensibly test M for SPD, but 
      // at least get domain right
      if ((M.getDomain() != A.getDomain()) ||
	  (M.getRange()  != A.getDomain())) {
	RVLException e;
	e<<"Error PsiScaleCGNEStep constructor\n";
	e<<"  preconditioning operator does not have domain, range\n";
	e<<"  same as domain of system operator\n";
	e<<"  preconditioning operator:\n";
	M.write(e);
	e<<"  system operator:\n";
	A.write(e);
	throw e;
      }
      // NOTE: initial x assumed to be zero vector
      r.copy(b);                   // r0=b-Ax0 with x0=0
      rnorm=r.norm();          
      A.applyAdjOp(r,ng);          // ng = A^T r0
      //modify according to PsiDO 
      //M.applyOp(ng,g);             // g  = M A^T r0
      S.applyOp(ng,sng);             // sng= S ng
      A.applyOp(sng,Aimg);     
      A.applyAdjOp(Aimg,sg);
      S.applyOp(sg,swg);
      M.applyOp(sng,swg);
      S2.applyOp(swg,g);
      //S.applyOp(sg,g);
      //M.applyOp(sng,g);             // g  = M A^T r0
      p.copy(g);                   // p  = g
      gamma = abs(g.inner(ng));
      nrnorm=sqrt(gamma);
    }
      
    /**
       Run a single step of the conjugate gradients for the normal equations
    */
    void run() {
      try {
	A.applyOp(p,q);
	atype qtq = q.normsq();
	atype absalpha;
	if (ProtectedDivision<atype>(gamma,qtq,absalpha)) {
	  RVLException e;
	  e<<"Error: PsiScaleCGNEStep::run() from ProtectedDivision: alpha\n";
	  throw e;
	}
	// can use q for workspace since it is not needed again until
	// reinitialized

	Scalar alpha=absalpha;
	x.linComb(alpha,p);
	r.linComb(-alpha,q);

        preng.copy(ng);
	A.applyAdjOp(r,ng);
        S.applyOp(ng,sng);             // sng= S ng
        A.applyOp(sng,Aimg);     
        A.applyAdjOp(Aimg,sg);
        S.applyOp(sg,swg);
        M.applyOp(sng,swg);
        S2.applyOp(swg,g);
        //S.applyOp(sg,g);
        //M.applyOp(sng,g);             // g  = M A^T r0

	atype newgamma = abs(g.inner(ng));
        atype tmp = newgamma - abs(g.inner(preng));        
	atype absbeta;
	if (ProtectedDivision<atype>(tmp,gamma,absbeta)) {
	  RVLException e;
	  e<<"Error: PsiScaleCGNEStep::run() from ProtectedDivision: beta\n";
	  throw e;
	}

	Scalar beta = absbeta;
	p.linComb(ScalarFieldTraits<Scalar>::One(),g,beta);
	gamma=newgamma;
	rnorm=r.norm();
	nrnorm=sqrt(gamma);
      }
      catch (RVLException & e) {
	e<<"\ncalled from PsiScaleCGNEStep::run()\n";
	throw e;
      }
     
    }

    ~PsiScaleCGNEStep() {}

  protected:

    // references to external objects
    LinearOp<Scalar> const & A;
    LinearOp<Scalar> const & M;
    LinearOp<Scalar> const & S;
    LinearOp<Scalar> const & S2;
    Vector<Scalar> & x;
    Vector<Scalar> const & b;
    atype & rnorm;
    atype & nrnorm;

  private:

    // need four work vectors and one scalar as persistent object data
    Vector<Scalar> r;    // residual
    Vector<Scalar> ng;   // normal residual
    Vector<Scalar> g;    // gradient = preconditioned normal residual
    Vector<Scalar> q;    // image of search direction
    Vector<Scalar> p;    // search direction
    atype gamma;         // gradient norm, inner product defined by inverse of M
    // working vectors for PsiDO
    Vector<Scalar> Aimg;
    Vector<Scalar> preng;
    Vector<Scalar> sng;   // normal residual
    Vector<Scalar> sg;    // gradient = preconditioned normal residual
    Vector<Scalar> swg;   
  };



  /** Conjugate gradient algorithm - efficient implementation for
      normal equations
      \f[ A^{\prime} A x = A^{\prime} b\f]
      for solving the linear least squares problem
      \f[ \min_{x} \vert A x - b \vert^2 \f].

      This is Algorithm CGLS as stated in Paige and
      Saunders, ACM TOMS vol. 8 pp. 43-72 1982 (see p. 57). We use
      variable names aping Paige and Saunder's notation insofar as
      possible. 

      Structure and function: Combines CGNEStep with a Terminator
      which displays iteration count, residual norm, and normal
      residual norm on output stream (constructor argument _str), and
      terminates if iteration count exceeds max or residual norm or
      normal residual norm fall below threshhold (default =
      10*sqrt(macheps)). Also terminates if the length of the solution
      vector exceeds a specified bound (maxstep argument to
      constructor). In this latter case, the computed step is
      projected onto the ball of radius maxstep centered at the
      initial estimate. This maximum step limit and projection turns
      the algorithm into an approximate trust region subproblem
      solver, similar to Steihaug-Toint. The default choice of maxstep
      is the max Scalar, which effectively turns off the trust region
      feature.

      Usage: construct CGNEAlg object by supplying appropriate
      arguments to constructor. On return from constructor, solution
      vector initialized to zero, residual norm to norm of RHS, and
      normal residual norm to norm of image of RHS under adjoint of
      operator. Then call run() method. Progress of iteration written
      on output unit. On return from run(), solution vector stores
      final estimate of solution, and residual norm and normal
      residual norm scalars have corresponding values.

      Typical Use: see <a href="../../testsrc/testcgne.cc">
      functional test source</a>.

      IMPORTANT NOTE: This class is also an RVLAlg::Terminator
      subclass. Its query() method returns true if the trust region
      constraint was binding (raw LS solution larger than trust
      radius), else false.

      IMPORTANT NOTE: The solution vector and residual and normal
      residual scalars are external objects, for which this algorithm
      stores mutable references. These objects are updated by
      constructing a CGNEAlg object, and by calling its run() method.

      IMPORTANT NOTE: this version of the algorithm initializes the
      solution vector to zero. To accommodate nontrivial initial
      guess, <i>modify the right-hand-side vector</i> (argument _rhs)
      externally.

      IMPORTANT NOTE: in order that this algorithm function properly
      for complex scalar types, a careful distinction is maintained
      between the main template parameter (Scalar) type and its
      absolute value type. All of the scalars appearing in the
      algorithm are actually of the latter type.

      See constructor documentation for description of parameters.


  */

  template<typename Scalar>
  class PsiCGNEAlg: public Algorithm, public Terminator {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:


    /** Constructor - preconditioned algorithm

	@param x - mutable reference to solution vector (external),
	initialized to zero vector on construction, estimated solution
	on return from CGNEAlg::run().

	@param A - const reference to LinearOp (external) defining
	problem

	@param M - const reference to LinearOp (external) preconditioner 
	= inverse of operator defining inner product in domain space. Only this
	operator is necessary, not the inner product operator itself. M is 
	assumed to be SPD, else no guarantees about the behaviour of this algorithm

	@param rhs - const reference to RHS or target vector
	(external)

	@param _rnorm - mutable reference to residual norm scalar
	(external), initialized to norm of RHS on construction, norm
	of estimated residual at solution on return from
	CGNEAlg::run()

	@param _nrnorm - mutable reference to normal residual (least
	squares gradient) norm scalar (external), initialized to morm
	of image of RHS under adjoint of problem LinearOp on
	construction, norm of estimated normal residual at solution on
	return from CGNEAlg::run()

	@param _rtol - stopping threshold for residual norm, default
	value = 100.0*macheps

	@param _nrtol - stopping threshold for normal residual norm,
	default value = 100.0*macheps

	@param _maxcount - max number of iterations permitted, default
	value = 10

	@param _maxstep - max permitted step length (trust radius),
	default value = max absval scalar (which makes the trust
	region feature inactive)

	@param _str - output stream

     */
    PsiCGNEAlg(RVL::Vector<Scalar> & x, 
	    LinearOp<Scalar> const & A, 
	    LinearOp<Scalar> const & M, 
	    Vector<Scalar> const & rhs, 
	    atype & _rnorm,
	    atype & _nrnorm,
	    atype _rtol = 100.0*numeric_limits<atype>::epsilon(),
	    atype _nrtol = 100.0*numeric_limits<atype>::epsilon(),
	    int _maxcount = 10,
	    atype _maxstep = numeric_limits<atype>::max(),
	    ostream & _str = cout)
    : rnorm(_rnorm), 
      nrnorm(_nrnorm), 
      rtol(_rtol), 
      nrtol(_nrtol), 
      maxstep(_maxstep), 
      maxcount(_maxcount), 
      count(0), 
      proj(false), 
      str(_str)
    { 
      try {
	x.zero(); 
	step = new PsiCGNEStep<Scalar>(A,M,x,rhs,rnorm,nrnorm);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEAlg constructor (preconditioned)\n";
	throw e;
      }
    }
    PsiCGNEAlg(RVL::Vector<Scalar> & x, 
	    LinearOp<Scalar> const & A, 
	    LinearOp<Scalar> const & M, 
	    LinearOp<Scalar> const & S, 
	    LinearOp<Scalar> const & S2, 
	    Vector<Scalar> const & rhs, 
	    atype & _rnorm,
	    atype & _nrnorm,
	    atype _rtol = 100.0*numeric_limits<atype>::epsilon(),
	    atype _nrtol = 100.0*numeric_limits<atype>::epsilon(),
	    int _maxcount = 10,
	    atype _maxstep = numeric_limits<atype>::max(),
	    ostream & _str = cout)
    : rnorm(_rnorm), 
      nrnorm(_nrnorm), 
      rtol(_rtol), 
      nrtol(_nrtol), 
      maxstep(_maxstep), 
      maxcount(_maxcount), 
      count(0), 
      proj(false), 
      str(_str)
    { 
      try {
	x.zero(); 
	step = new PsiScaleCGNEStep<Scalar>(A,M,S,S2,x,rhs,rnorm,nrnorm);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEAlg constructor (preconditioned)\n";
	throw e;
      }
    }

    ~PsiCGNEAlg() {
      delete step;
    }

    bool query() { return proj; }

    void run() { 
      try {
	//	cerr<<"cg::run\n";
	// access to internals
	PsiCGNEStep<Scalar> * pcg = dynamic_cast<PsiCGNEStep<Scalar> *>(step);
	PsiScaleCGNEStep<Scalar> * pscg = dynamic_cast<PsiScaleCGNEStep<Scalar> *>(step);
	if ((!pscg && !pcg) || (pscg && pcg)) {
	  RVLException e;
	  e<<"Error: PsiCGNEAlg::run\n";
	  e<<"  unable to determine whether step data member\n";
	  e<<"  is scaled or not. PANIC!\n";
	  throw e;
	}
	
	// terminator for CGNE iteration
	vector<string> names(2);
	vector<atype *> nums(2);
	vector<atype> tols(2);
	atype rnorm0=rnorm;
	atype nrnorm0=nrnorm;
          
	names[0]="Residual Norm"; nums[0]=&rnorm; tols[0]=rtol*rnorm0;
	names[1]="Gradient Norm"; nums[1]=&nrnorm; tols[1]=nrtol*nrnorm0;
	if (pcg) str<<"========================== BEGIN PsiCGNE =========================\n";
	else 	str<<"========================== BEGIN PsiScaleCGNE =========================\n";

	VectorCountingThresholdIterationTable<atype> stop1(maxcount,names,nums,tols,str);
	stop1.init();
	// terminator for Trust Region test and projection
	//      BallProjTerminator<Scalar> stop2(x,maxstep,str);

	Terminator * stop2 = NULL;
	if (pcg)  stop2 = new BallProjTerminator<Scalar>(pcg->x,maxstep,str); 
	if (pscg) stop2 = new BallProjTerminator<Scalar>(pscg->x,maxstep,str); 
	// terminate if either
	OrTerminator stop(stop1,*stop2);
	// loop
	LoopAlg doit(*step,stop);
	//	cerr<<"stop1="<<stop1.query()<<" stop2="<<stop2->query()<<endl;
	//	cerr<<"rtol="<<rtol<<" nrtol="<<nrtol<<" maxcount="<<maxcount<<endl;
	doit.run();
	// must recompute residual if scaling occured 
	//	cerr<<"stop1="<<stop1.query()<<" stop2="<<stop2->query()<<" xnorm="<<cg->x.norm()<<" maxstep="<<maxstep<<"\n";
	proj = stop2->query();
	if (proj) {
	  if (pcg) {
	    Vector<Scalar> temp((pcg->A).getRange());
	    (pcg->A).applyOp((pcg->x),temp);
	    temp.linComb(-1.0,(pcg->b));
	    rnorm=temp.norm();
	    Vector<Scalar> temp1((pcg->A).getDomain());
	    Vector<Scalar> temp2((pcg->A).getDomain());
	    (pcg->A).applyAdjOp(temp,temp1);
	    (pcg->M).applyOp(temp1,temp2);
	    atype tmpgamma = abs(temp1.inner(temp2));
	    nrnorm=sqrt(tmpgamma);
	  }
	  if (pscg) {
//	    Vector<Scalar> temp((pcg->A).getRange());
//	    (pcg->A).applyOp((pcg->x),temp);
//	    temp.linComb(-1.0,(pcg->b));
//	    rnorm=temp.norm();
//	    Vector<Scalar> temp1((pcg->A).getDomain());
//	    Vector<Scalar> temp2((pcg->A).getDomain());
//	    (pcg->A).applyAdjOp(temp,temp1);
//	    (pcg->M).applyOp(temp1,temp2);
//	    atype tmpgamma = abs(temp1.inner(temp2));
//	    nrnorm=sqrt(tmpgamma);
	  }
	}

	count = stop1.getCount();
	if (pcg) str<<"=========================== END PsiCGNE ==========================\n";
	else str<<"=========================== END PsiScaleCGNE ==========================\n";

          // display results
          str<<"\n ****************** PsiCGNE summary *****************  "<<endl;
          str<<"initial residual norm      = "<<rnorm0<<endl;
          str<<"residual norm              = "<<rnorm<<endl;
          str<<"residual redn              = "<<rnorm/rnorm0<<endl;
          str<<"initial gradient norm      = "<<nrnorm0<<endl;
          str<<"gradient norm              = "<<nrnorm<<endl;
          str<<"gradient redn              = "<<nrnorm/nrnorm0<<endl;
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGNEAlg::run\n";
	throw e;
      }
    }

    int getCount() const { return count; }

  private:

    atype & rnorm;                 // residual norm
    atype & nrnorm;                // gradient norm
    atype rtol;                    // tolerance residual norm
    atype nrtol;                   // tolerance gradient norm 
    atype maxstep;                 // upper bound for net step x-x0
    int maxcount;                  // upper bound for iteration count
    int count;                     // actual iteration count
    mutable bool proj;             // whether step is projected onto TR boundary
    ostream & str;                 // stream for report output
    Algorithm * step;              // CGNE or PsiCGNE step

    // disable default, copy constructors
    PsiCGNEAlg();
    PsiCGNEAlg(PsiCGNEAlg<Scalar> const &);

  };

  /** data class for PsiCGNE policy
  */
  template<typename Scalar>
  class PsiCGNEPolicyData {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
  
  public:

    atype rtol;
    atype nrtol;
    atype Delta;
    int maxcount;
    bool verbose;

    PsiCGNEPolicyData(atype _rtol = numeric_limits<atype>::max(),
		   atype _nrtol = numeric_limits<atype>::max(),
		   atype _Delta = numeric_limits<atype>::max(),
		   int _maxcount = 0,
		   bool _verbose = false)
      : rtol(_rtol), nrtol(_nrtol), Delta(_Delta), maxcount(_maxcount), verbose(_verbose) {}
      
    PsiCGNEPolicyData(PsiCGNEPolicyData<Scalar> const & a) 
      : rtol(a.rtol), nrtol(a.nrtol), Delta(a.Delta), maxcount(a.maxcount), verbose(a.verbose) {}

    ostream & write(ostream & str) const {
      str<<"\n";
      str<<"==============================================\n";
      str<<"PsiCGNEPolicyData: \n";
      str<<"rtol      = "<<rtol<<"\n";
      str<<"nrtol     = "<<nrtol<<"\n";
      str<<"Delta     = "<<Delta<<"\n";
      str<<"maxcount  = "<<maxcount<<"\n";
      str<<"verbose   = "<<verbose<<"\n";
      str<<"==============================================\n";
      return str;
    }
  };

  /** policy class for creation of CGNEAlg in trust region solver and 
      any other algorithm needing a least squares solver component - build
      method creates CGNEAlg with these attributes:

      rtol     = residual threshhold for convergence
      nrtol    = normal residual (LS gradient) tolerance for convergence
      Delta    = trust radius - truncate iteration when reached
      maxcount = max number of iterations permitted

      Default values set to cause immediate return from CGNEAlg::run.

      Other attributes are arguments of build method.

      Conforms to specifications described in TRGNAlg docs.

      Usage: use as policy type, i.e. class name as template parameter
      to TRGNAlg constructor. After construction of TRGNAlg, which IS
      a CGNEPolicy by inheritance, call the assign method on it to set
      rtol, nrtol, and maxcount, BEFORE calling TRGNAlg::run() - else
      you will get immediate termination, as intended.
  */

  template<typename Scalar> 
  class PsiCGNEPolicy {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:
      
    PsiCGNEAlg<Scalar> * build(Vector<Scalar> & x,
                            LinearOp<Scalar> const & A,
                            LinearOp<Scalar> const & M,
                            Vector<Scalar> const & d,
                            atype & rnorm,
                            atype & nrnorm,
                            ostream & str) const {
    try {
        if (verbose)
            return new PsiCGNEAlg<Scalar>(x,A,M,d,rnorm,nrnorm,rtol,nrtol,maxcount,Delta,str);
        else
            return new PsiCGNEAlg<Scalar>(x,A,M,d,rnorm,nrnorm,rtol,nrtol,maxcount,Delta,nullstr);
        }
        catch (RVLException & e) {
            e<<"\ncalled from PsiCGNEPolicy::build\n";
            e<<"inputs: \n";
            e<<"**** x:\n";
            x.write(e);
            e<<"**** A:\n";
            A.write(e);
            e<<"**** d:\n";
            d.write(e);
            throw e;
        }
    }
      
    /** post-construction initialization
	@param _rtol - residual norm stopping threshhold
	@param _nrtol - normal residual (LS gradient) norm stopping threshhold
	@param _maxcount - max number of permitted iterations
    */
    void assign(atype _rtol, atype _nrtol, atype _Delta, int _maxcount, bool _verbose) {
      rtol=_rtol; nrtol=_nrtol; Delta=_Delta; maxcount=_maxcount; verbose=_verbose;
    }

    /** parameter table overload */
    void assign(Table const & t) {
      rtol=getValueFromTable<atype>(t,"CGNE_ResTol");
      nrtol=getValueFromTable<atype>(t,"CGNE_GradTol"); 
      Delta=getValueFromTable<atype>(t,"TR_Delta");
      maxcount=getValueFromTable<int>(t,"CGNE_MaxItn"); 
      verbose=getValueFromTable<bool>(t,"CGNE_Verbose");
    }

    /** data struct overload */
    void assign(PsiCGNEPolicyData<Scalar> const & s) {
      rtol=s.rtol;
      nrtol=s.nrtol;
      Delta=s.Delta;
      maxcount=s.maxcount;
      verbose=s.verbose;
    }

    /** only Delta need be changed repeatedly, as opposed
	to set post-construction. Simplest way to do this - make
	it public
    */
    mutable atype Delta;

    /** main constructor - acts as default. Default values of
	parameters set to result in immediate return, no
	iteration. Note that policy design requires that default
	construction must be valid, and all run-time instance data be
	initiated post-construction, in this case by the assign
	function, to be called by drivers of user classes (subclassed
	from this and with this as template param).
    */
    PsiCGNEPolicy(atype _rtol = numeric_limits<atype>::max(),
	       atype _nrtol = numeric_limits<atype>::max(),
	       atype _Delta = numeric_limits<atype>::max(),
	       int _maxcount = 0,
	       bool _verbose = true)
      : Delta(_Delta), rtol(_rtol), nrtol(_nrtol), maxcount(_maxcount), verbose(_verbose), nullstr(0) {}

    PsiCGNEPolicy(PsiCGNEPolicy<Scalar> const & p)
      : 	Delta(p.Delta), 
		rtol(p.rtol), 
		nrtol(p.nrtol), 
		maxcount(p.maxcount), 
		verbose(p.verbose), 
		nullstr(0) {}
      
  private:
    mutable atype rtol;
    mutable atype nrtol;
    mutable int maxcount;
    mutable bool verbose;
    mutable std::ostream nullstr;
  };
}

#endif









