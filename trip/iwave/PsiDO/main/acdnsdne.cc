#include "acd_defn.hh"
#include "grid.h"
#include "iwsim.hh"
#include "iwop.hh"
#include "gridpp.hh"
#include "gridops.hh"
#include "functions.hh"
#include "PsiDOsdnealg.hh"
#include "adjtest.hh"
#include "segyops.hh"

//#define GTEST_VERBOSE
IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, true,  true },
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
using TSOpt::GridMaskOp;
using TSOpt::IWaveEnvironment;
using TSOpt::IWaveTree;
using TSOpt::IWaveSampler;
using TSOpt::IWaveSim;
using TSOpt::TASK_RELN;
using TSOpt::IOTask;
using TSOpt::IWaveOp;
using TSOpt::SEGYLinMute;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif

using RVLUmin::PsiSDNEAlg;

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

    AssignFilename mfn(valparse<std::string>(*pars,"rcsq"));
    m.eval(mfn);

    AssignFilename dmfn(valparse<std::string>(*pars,"icsq"));
    dm.eval(dmfn);
    dm.zero();

    AssignFilename ddfn(valparse<std::string>(*pars,"data"));
    dd.eval(ddfn);

    RPNT swind,ewind,width;
    RASN(swind,RPNT_0);
    RASN(ewind,RPNT_0);
    swind[0]=valparse<float>(*pars,"sww0",0.0f);
    swind[1]=valparse<float>(*pars,"sww1",0.0f);
    swind[2]=valparse<float>(*pars,"sww2",0.0f);
    ewind[0]=valparse<float>(*pars,"eww0",0.0f);
    ewind[1]=valparse<float>(*pars,"eww1",0.0f);
    ewind[2]=valparse<float>(*pars,"eww2",0.0f);
    width[0]=valparse<float>(*pars,"width0",0.0f);
    width[1]=valparse<float>(*pars,"width1",0.0f);
    width[2]=valparse<float>(*pars,"width2",0.0f);

    // need to read in model space for bg input to GridWindowOp
    Vector<ireal> m_in(op.getDomain());
    AssignFilename minfn(valparse<std::string>(*pars,"csqext"));
    Components<ireal> cmin(m_in);
    cmin[0].eval(minfn);
    GridMaskOp mop(op.getDomain(),m_in,swind,ewind);

    OperatorEvaluation<float> mopeval(mop,m_in);
    LinearOp<float> const & lmop=mopeval.getDeriv();

    float rtol=valparse<float>(*pars,"ResidualTol",100.0*numeric_limits<float>::epsilon());
    float nrtol=valparse<float>(*pars,"GradientTol",100.0*numeric_limits<float>::epsilon());
    int maxcount=valparse<int>(*pars,"MaxIter",10);
    float maxstep=valparse<float>(*pars,"MaxStep",numeric_limits<float>::max());

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
    
    /* create CGNE object */
    float rnorm;
    float nrnorm;
    OperatorEvaluation<ireal> opeval(op,m);
    int appw=valparse<int>(*pars,"appwind",1);
    if(appw){
      CompLinearOp<float> cop(lmop,opeval.getDeriv());
      PsiSDNEAlg<float> alg(dm,cop,dd,
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
    }
    else {
      PsiSDNEAlg<float> alg(dm,opeval.getDeriv(),dd,
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
    }



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
	dres.copy(dd);
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
