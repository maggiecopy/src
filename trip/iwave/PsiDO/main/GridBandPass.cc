#include "grid.h"
#include "functions.hh"
#include "adjtest.hh"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif
#include "gridops.hh"
#include "gridbandpassop.hh"

using RVL::parse;
using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::Operator;
using RVL::OperatorEvaluation;
using RVL::LinearOp;
using RVL::LinearOpFO;
using RVL::AssignFilename;
using RVL::AssignParams;
using RVL::RVLRandomize;
using RVL::AdjointTest;
using TSOpt::GridWindowOp;
using TSOpt::GridBandPassOp;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
typedef TSOpt::MPIGridSpace gsp;
#else
using TSOpt::GridSpace;
typedef TSOpt::GridSpace gsp;
#endif


int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
    storeGlobalComm(MPI_COMM_WORLD);
#endif

    PARARRAY * pars = pars = ps_new();

    if (ps_createargs(pars,argc-1,&(argv[1]))) {
      RVLException e;
      e<<"ERROR: GridDerivOp from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++)
        e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }

    string inp = valparse<string>(*pars,"in");
    string outp = valparse<string>(*pars,"out");

    gsp sp(inp,"notype",true
#ifdef IWAVE_USE_MPI 
           , retrieveGlobalComm()
#endif
           ); 

    Vector<ireal> min(sp);
    Vector<ireal> mout(sp);

    AssignFilename minfn(inp);
    Components<ireal> cmin(min);
    cmin[0].eval(minfn);
      
    AssignFilename moutfn(outp);
    Components<ireal> cmout(mout);
    cmout[0].eval(moutfn);

    Vector<ireal> mvscale(sp);
    string vscale = valparse<string>(*pars,"vscale");
    if (vscale.size()>0) {
      AssignFilename mvfn(vscale);
      mvscale.eval(mvfn);
    }
    else mvscale.zero();

    GridBandPassOp wop(sp,mvscale,
                       valparse<float>(*pars,"flo",0.f),
                       valparse<float>(*pars,"fhi",0.05f),
                       valparse<bool>(*pars,"phase",false));

    RVLRandomize<float> rnd(getpid(),-1.0,1.0);
    
    //AdjointTest<float>(wop,rnd,cerr);
    
    if (valparse<bool>(*pars,"adjoint")) {
      wop.applyAdjOp(min,mout);
    }
    else {
      wop.applyOp(min,mout);
    }
    ps_delete(&pars);    
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
