#include "acd_defn.hh"
#include "grid.h"
#include "iwsim.hh"
#include "iwop.hh"
#include "gridpp.hh"
#include "functions.hh"
#include "PsiDOcgnealg.hh"
#include "psido_selfdoc.h"
#include "psidiffop.hh"
#include "adjtest.hh"
#include "segyops.hh"

//#define GTEST_VERBOSE
IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csqext",    0, true,  true },
  {"data",   1, false, true },
  {"source", 1, true,  false},
  {"",       0, false, false}
};

using RVL::parse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::Operator;
using RVL::OperatorEvaluation;
using RVL::LinearOp;
using RVL::LinearOpFO;
using RVL::OpComp;
using RVL::SymmetricBilinearOp;
using RVL::AssignFilename;
using RVL::AssignParams;
using RVL::RVLRandomize;
using RVL::AdjointTest;
using TSOpt::IWaveEnvironment;
using TSOpt::IWaveTree;
using TSOpt::IWaveSampler;
using TSOpt::IWaveSim;
using TSOpt::TASK_RELN;
using TSOpt::IOTask;
using TSOpt::IWaveOp;
using TSOpt::SEGYLinMute;
#ifdef IWAVE_USE_MPI
using TSOpt::MPISEGYSpace;
#else
using TSOpt::SEGYSpace;
#endif

using TSOpt::PsiDiffOp;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
typedef TSOpt::MPIGridSpace gsp;
#else
using TSOpt::GridSpace;
typedef TSOpt::GridSpace gsp;
#endif

using RVLUmin::PsiCGNEAlg;

int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

    PARARRAY * pars = NULL;
    FILE * stream = NULL;
    IWaveEnvironment(argc, argv, 0, &pars, &stream);
    
    // the Op
    IWaveOp op(*pars,stream);
    
    Vector<ireal> m(op.getDomain());
    Vector<ireal> dm(op.getDomain());
    Vector<ireal> dd(op.getRange());
    Vector<ireal> mdd(op.getRange());

    AssignFilename mfn(valparse<std::string>(*pars,"rcsq"));
    m.eval(mfn);

    AssignFilename dmfn(valparse<std::string>(*pars,"icsq"));
    dm.eval(dmfn);
    dm.zero();

    AssignFilename ddfn(valparse<std::string>(*pars,"data"));
    dd.eval(ddfn);

    float rtol=valparse<float>(*pars,"ResidualTol",100.0*numeric_limits<float>::epsilon());
    float nrtol=valparse<float>(*pars,"GradientTol",100.0*numeric_limits<float>::epsilon());
    int maxcount=valparse<int>(*pars,"MaxIter",10);
    float maxstep=valparse<float>(*pars,"MaxStep",numeric_limits<float>::max());

    RVLRandomize<float> rnd(getpid(),-1.0,1.0);
    
    /* output stream */
    std::stringstream res;
    res<<scientific;
    
    res<<endl<<"*******************************************************"<<endl;
    res<<"* Acoustic Constant Density Linearized Inversion via";
    res<<"* Conjugate Gradient Algorithm for Normal Eqns"<<endl;
    res<<"* max iterations       = "<<maxcount<<endl;
    res<<"* residual tolerance   = "<<rtol<<endl;
    res<<"* normal res tolerance = "<<nrtol<<endl;
    res<<"* trust radius         = "<<maxstep<<endl;
    res<<"*******************************************************"<<endl;
    
        string coeff   = valparse<string>(*pars,"coeff");
        string csqext   = valparse<string>(*pars,"csqext");

        gsp csp(coeff, "coeff", true
#ifdef IWAVE_USE_MPI 
           , retrieveGlobalComm()
#endif
           );

        Vector<float> incoeff(csp);
        AssignFilename afincoeff(coeff);
        incoeff.eval(afincoeff);
        /* PsiDO op */
        PsiDiffOp psido(op.getDomain(),incoeff,*pars,stream,1);


    /* create CGNE object */
    float rnorm;
    float nrnorm;
    OperatorEvaluation<ireal> opeval(op,m);
    //    AdjointTest<float>(opeval.getDeriv(),rnd,cerr);
    PsiCGNEAlg<float> alg(dm,opeval.getDeriv(),psido, dd,
		       rnorm, nrnorm, rtol, nrtol, maxcount, maxstep, res);
    float nrnorm0=nrnorm;
    float rnorm0=rnorm;
    
    alg.run();
    
    // display results
    res<<"\n ******* summary ********  "<<endl;
    res<<"initial residual norm      = "<<rnorm0<<endl;
    res<<"residual norm              = "<<rnorm<<endl;
    res<<"residual redn              = "<<rnorm/rnorm0<<endl;
    res<<"initial gradient norm      = "<<nrnorm0<<endl;
    res<<"gradient norm              = "<<nrnorm<<endl;
    res<<"gradient redn              = "<<nrnorm/nrnorm0<<endl;
    
    std::string dataest = valparse<std::string>(*pars,"dataest","");
    std::string datares = valparse<std::string>(*pars,"datares","");
    if (dataest.size()>0) {
      Vector<float> est(op.getRange());
      AssignFilename estfn(dataest);
      est.eval(estfn);
      opeval.getDeriv().applyOp(dm,est);
      if (datares.size()>0) {
	Vector<float> dres(op.getRange());
	AssignFilename resfn(datares);
	dres.eval(resfn);
	dres.copy(mdd);
	dres.linComb(-1.0f,est);
      } 
    }

    if (retrieveRank() == 0) {
      std::string outfile = valparse<std::string>(*pars,"outfile","");
      if (outfile.size()>0) {
	ofstream outf(outfile.c_str());
	outf<<res.str();
	outf.close();
      }
      else {
	cout<<res.str();
      }
    }
    
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
  
}
