// PsiDOsdalg.H
// created by Yin Huang 08/05/15

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


#ifndef __PsiDO_SDNE_H
#define __PsiDO_SDNE_H

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
  class PsiSDNEAlg;

  template<typename Scalar>
  class SDNEStep : public Algorithm {

    friend class PsiSDNEAlg<Scalar>;

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:

    SDNEStep(LinearOp<Scalar> const & _A,
	      Vector<Scalar> & _x,
	      Vector<Scalar> const & _b,
	      atype & _rnorm, 
	      atype & _nrnorm)
      : A(_A), x(_x), b(_b), rnorm(_rnorm), nrnorm(_nrnorm), 
	r(A.getRange()), g(A.getDomain()),
	q(A.getRange()) {
      // NOTE: initial x assumed to be zero vector
      r.copy(b);                   // r0=b-Ax0 with x0=0
      rnorm=r.norm();          
      A.applyAdjOp(r,g);          // ng = A^T r0
      gamma = abs(g.inner(g));
      nrnorm=sqrt(gamma);
    }
      
    /**
       Run a single step of the conjugate gradients for the normal equations
    */
    void run() {
      try {
	A.applyOp(g,q);
	atype qtq = q.normsq();
	atype absalpha;
	if (ProtectedDivision<atype>(gamma,qtq,absalpha)) {
	  RVLException e;
	  e<<"Error: SDNEStep::run() from ProtectedDivision: alpha\n";
	  throw e;
	}
	// can use q for workspace since it is not needed again until
	// reinitialized

	Scalar alpha=absalpha;
	x.linComb(alpha,g);
	r.linComb(-alpha,q);

        A.applyAdjOp(r,g);
	gamma = abs(g.inner(g));
	rnorm=r.norm();
	nrnorm=sqrt(gamma);
      }
      catch (RVLException & e) {
	e<<"\ncalled from SDNEStep::run()\n";
	throw e;
      }
     
    }

    ~SDNEStep() {}

  protected:

    // references to external objects
    LinearOp<Scalar> const & A;
    Vector<Scalar> & x;
    Vector<Scalar> const & b;
    atype & rnorm;
    atype & nrnorm;

  private:

    // need four work vectors and one scalar as persistent object data
    Vector<Scalar> r;    // residual
    Vector<Scalar> g;    // minus gradient A^T(b-Ax)
    Vector<Scalar> q;    // image of search direction
    atype gamma;         
  };

  /** Preconditioned steepest descent method for the normal
      equations with PsiDO as a preconditioner.
  */
  template<typename Scalar>
  class PsiSDNEStep : public Algorithm {

    friend class PsiSDNEAlg<Scalar>;

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:

    PsiSDNEStep(LinearOp<Scalar> const & _A,
	      LinearOp<Scalar> const & _M,
	      Vector<Scalar> & _x,
	      Vector<Scalar> const & _b,
	      atype & _rnorm, 
	      atype & _nrnorm)
      : A(_A), M(_M), x(_x), b(_b), rnorm(_rnorm), nrnorm(_nrnorm), 
	r(A.getRange()), ng(A.getDomain()), g(A.getDomain()),
	q(A.getRange()), Aimg(A.getRange()) { 
      // sanity tests - cannot sensibly test M for SPD, but 
      // at least get domain right
      if ((M.getDomain() != A.getDomain()) ||
	  (M.getRange()  != A.getDomain())) {
	RVLException e;
	e<<"Error PsiSDNEStep constructor\n";
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
      gamma = abs(g.inner(ng));
      nrnorm=sqrt(gamma);
    }
      
    /**
       Run a single step of the conjugate gradients for the normal equations
    */
    void run() {
      try {
	A.applyOp(g,q);
	atype qtq = q.normsq();
	atype absalpha;
	if (ProtectedDivision<atype>(gamma,qtq,absalpha)) {
	  RVLException e;
	  e<<"Error: SDNEStep::run() from ProtectedDivision: alpha\n";
	  throw e;
	}
	// can use q for workspace since it is not needed again until
	// reinitialized

	Scalar alpha=absalpha;
	x.linComb(alpha,g);
	r.linComb(-alpha,q);

        // need to call A 3 times because of PsiDO
	A.applyAdjOp(r,ng);
        A.applyOp(ng,Aimg);
        A.applyAdjOp(Aimg,g);
	M.applyOp(ng,g);

	gamma = abs(ng.inner(g));
	rnorm=r.norm();
	nrnorm=sqrt(gamma);
      }
      catch (RVLException & e) {
	e<<"\ncalled from PsiSDNEStep::run()\n";
	throw e;
      }
     
    }

    ~PsiSDNEStep() {}

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
    atype gamma;         // gradient norm, inner product defined by inverse of M
    // working vectors for PsiDO
    Vector<Scalar> Aimg;
  };

/** Steepest descent algorithm - efficient implementation for
      normal equations
      \f[ A^{\prime} A x = A^{\prime} b\f]
      for solving the linear least squares problem
      \f[ \min_{x} \vert A x - b \vert^2 \f].

      This is Steepest descent Algorithm: 

      Step 1:
        \f $ r_0=b-Ax_0 $;
        \f $ rnorm=\|r\| $;
        \f $ g_0=A^Tr_0 $;
        \f $ nrnorm=\|g\| $;

      Step 2:
        for \f$ k=0,1,\cdot $
            \f$ q_k=Ag_k $; 
            \f$ \alpha=\frac{\|g_k\|^2}{\|q_k\|^2}$;
            \f$ x_{k+1}=x_k+\alpha g_k $;
            \f$ r_{k+1}=r_k-\alpha q_k $;
            \f$ g_{k+1}=A^Tr_{k+1} $;

      IMPORTANT NOTE: This class is also an RVLAlg::Terminator
      subclass. Its query() method returns true if the trust region
      constraint was binding (raw LS solution larger than trust
      radius), else false.

      IMPORTANT NOTE: The solution vector and residual and normal
      residual scalars are external objects, for which this algorithm
      stores mutable references. These objects are updated by
      constructing a SDNEAlg object, and by calling its run() method.

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
  class PsiSDNEAlg: public Algorithm, public Terminator {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:


    /** Constructor - preconditioned algorithm

	@param x - mutable reference to solution vector (external),
	initialized to zero vector on construction, estimated solution
	on return from SDNEAlg::run().

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
	SDNEAlg::run()

	@param _nrnorm - mutable reference to normal residual (least
	squares gradient) norm scalar (external), initialized to morm
	of image of RHS under adjoint of problem LinearOp on
	construction, norm of estimated normal residual at solution on
	return from SDNEAlg::run()

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
    PsiSDNEAlg(RVL::Vector<Scalar> & x, 
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
	step = new PsiSDNEStep<Scalar>(A,M,x,rhs,rnorm,nrnorm);
      }
      catch (RVLException & e) {
	e<<"\ncalled from SDNEAlg constructor (preconditioned)\n";
	throw e;
      }
    }
    PsiSDNEAlg(RVL::Vector<Scalar> & x, 
	    LinearOp<Scalar> const & A, 
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
	step = new SDNEStep<Scalar>(A,x,rhs,rnorm,nrnorm);
      }
      catch (RVLException & e) {
	e<<"\ncalled from SDNEAlg constructor (preconditioned)\n";
	throw e;
      }
    }

    ~PsiSDNEAlg() {
      delete step;
    }

    bool query() { return proj; }

    void run() { 
      try {
	//	cerr<<"cg::run\n";
	// access to internals
	PsiSDNEStep<Scalar> * psd = dynamic_cast<PsiSDNEStep<Scalar> *>(step);
	SDNEStep<Scalar>    *  sd = dynamic_cast<   SDNEStep<Scalar> *>(step);
	if ((!sd && !psd) || (sd && psd)) {
	  RVLException e;
	  e<<"Error: PsiSDNEAlg::run\n";
	  e<<"  unable to determine whether step data member\n";
	  e<<"  is preconditioned or not. PANIC!\n";
	  throw e;
	}
	
	// terminator for SDNE iteration
	vector<string> names(2);
	vector<atype *> nums(2);
	vector<atype> tols(2);
	atype rnorm0=rnorm;
	atype nrnorm0=nrnorm;
          
	names[0]="Residual Norm"; nums[0]=&rnorm; tols[0]=rtol*rnorm0;
	names[1]="Gradient Norm"; nums[1]=&nrnorm; tols[1]=nrtol*nrnorm0;
	if (psd) str<<"========================== BEGIN PsiSDNE =========================\n";
	else 	 str<<"========================== BEGIN SDNE =========================\n";

	VectorCountingThresholdIterationTable<atype> stop1(maxcount,names,nums,tols,str);
	stop1.init();
	// terminator for Trust Region test and projection
	//      BallProjTerminator<Scalar> stop2(x,maxstep,str);

	Terminator * stop2 = NULL;
	if (psd) stop2 = new BallProjTerminator<Scalar>(psd->x,maxstep,str); 
	if (sd)  stop2 = new BallProjTerminator<Scalar>( sd->x,maxstep,str); 
	// terminate if either
	OrTerminator stop(stop1,*stop2);
	// loop
	LoopAlg doit(*step,stop);
	doit.run();
	// must recompute residual if scaling occured 
	proj = stop2->query();
	if (proj) {
	  if (psd) {
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
	  if (sd) {
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
	if (psd) str<<"=========================== END PsiSDNE ==========================\n";
	else     str<<"=========================== END SDNE ==========================\n";

          // display results
          str<<"\n ****************** PsiSDNE summary *****************  "<<endl;
          str<<"initial residual norm      = "<<rnorm0<<endl;
          str<<"residual norm              = "<<rnorm<<endl;
          str<<"residual redn              = "<<rnorm/rnorm0<<endl;
          str<<"initial gradient norm      = "<<nrnorm0<<endl;
          str<<"gradient norm              = "<<nrnorm<<endl;
          str<<"gradient redn              = "<<nrnorm/nrnorm0<<endl;
      }
      catch (RVLException & e) {
	e<<"\ncalled from SDNEAlg::run\n";
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
    Algorithm * step;              // SDNE or PsiSDNE step

    // disable default, copy constructors
    PsiSDNEAlg();
    PsiSDNEAlg(PsiSDNEAlg<Scalar> const &);

  };

  /** data class for PsiSDNE policy
  */
  template<typename Scalar>
  class PsiSDNEPolicyData {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
  
  public:

    atype rtol;
    atype nrtol;
    atype Delta;
    int maxcount;
    bool verbose;

    PsiSDNEPolicyData(atype _rtol = numeric_limits<atype>::max(),
		   atype _nrtol = numeric_limits<atype>::max(),
		   atype _Delta = numeric_limits<atype>::max(),
		   int _maxcount = 0,
		   bool _verbose = false)
      : rtol(_rtol), nrtol(_nrtol), Delta(_Delta), maxcount(_maxcount), verbose(_verbose) {}
      
    PsiSDNEPolicyData(PsiSDNEPolicyData<Scalar> const & a) 
      : rtol(a.rtol), nrtol(a.nrtol), Delta(a.Delta), maxcount(a.maxcount), verbose(a.verbose) {}

    ostream & write(ostream & str) const {
      str<<"\n";
      str<<"==============================================\n";
      str<<"PsiSDNEPolicyData: \n";
      str<<"rtol      = "<<rtol<<"\n";
      str<<"nrtol     = "<<nrtol<<"\n";
      str<<"Delta     = "<<Delta<<"\n";
      str<<"maxcount  = "<<maxcount<<"\n";
      str<<"verbose   = "<<verbose<<"\n";
      str<<"==============================================\n";
      return str;
    }
  };

  /** policy class for creation of SDNEAlg in trust region solver and 
      any other algorithm needing a least squares solver component - build
      method creates SDNEAlg with these attributes:

      rtol     = residual threshhold for convergence
      nrtol    = normal residual (LS gradient) tolerance for convergence
      Delta    = trust radius - truncate iteration when reached
      maxcount = max number of iterations permitted

      Default values set to cause immediate return from SDNEAlg::run.

      Other attributes are arguments of build method.

      Conforms to specifications described in TRGNAlg docs.

      Usage: use as policy type, i.e. class name as template parameter
      to TRGNAlg constructor. After construction of TRGNAlg, which IS
      a SDNEPolicy by inheritance, call the assign method on it to set
      rtol, nrtol, and maxcount, BEFORE calling TRGNAlg::run() - else
      you will get immediate termination, as intended.
  */

  template<typename Scalar> 
  class PsiSDNEPolicy {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:
      
    PsiSDNEAlg<Scalar> * build(Vector<Scalar> & x,
                            LinearOp<Scalar> const & A,
                            LinearOp<Scalar> const & M,
                            Vector<Scalar> const & d,
                            atype & rnorm,
                            atype & nrnorm,
                            ostream & str) const {
    try {
        if (verbose)
            return new PsiSDNEAlg<Scalar>(x,A,M,d,rnorm,nrnorm,rtol,nrtol,maxcount,Delta,str);
        else
            return new PsiSDNEAlg<Scalar>(x,A,M,d,rnorm,nrnorm,rtol,nrtol,maxcount,Delta,nullstr);
        }
        catch (RVLException & e) {
            e<<"\ncalled from PsiSDNEPolicy::build\n";
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
      rtol=getValueFromTable<atype>(t,"SDNE_ResTol");
      nrtol=getValueFromTable<atype>(t,"SDNE_GradTol"); 
      Delta=getValueFromTable<atype>(t,"TR_Delta");
      maxcount=getValueFromTable<int>(t,"SDNE_MaxItn"); 
      verbose=getValueFromTable<bool>(t,"SDNE_Verbose");
    }

    /** data struct overload */
    void assign(PsiSDNEPolicyData<Scalar> const & s) {
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
    PsiSDNEPolicy(atype _rtol = numeric_limits<atype>::max(),
	       atype _nrtol = numeric_limits<atype>::max(),
	       atype _Delta = numeric_limits<atype>::max(),
	       int _maxcount = 0,
	       bool _verbose = true)
      : Delta(_Delta), rtol(_rtol), nrtol(_nrtol), maxcount(_maxcount), verbose(_verbose), nullstr(0) {}

    PsiSDNEPolicy(PsiSDNEPolicy<Scalar> const & p)
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









