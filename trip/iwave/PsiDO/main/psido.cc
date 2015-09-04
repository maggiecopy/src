#include <stdio.h>
#include <cstring>
#include <fftw3.h>
#include "lbfgs.h"
#include "objective_real.h"
#include <time.h>
#include <iostream>
#include "parser.h"
#include "psido_selfdoc.h"
#include "psidiffop.hh"

#include "grid.h"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif
#include "gridops.hh"
#include "iwenv.hh"
#include "functions.hh"
#include "op.hh"
#include "acd_defn.hh"

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, true,  true },
  {"data",   1, false, true },
  {"source", 1, true,  false},
  {"",       0, false, false}
};


using RVL::RVLException;
using RVL::valparse;
using RVL::AssignFilename;
using TSOpt::IWaveEnvironment;
using TSOpt::PsiDiffOp;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
typedef TSOpt::MPIGridSpace gsp;
#else
using TSOpt::GridSpace;
typedef TSOpt::GridSpace gsp;
#endif


int main(int argc, char **argv){
    try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

    PARARRAY * pars = NULL;
    FILE * stream = NULL;
    IWaveEnvironment(argc, argv, 0, &pars, &stream);

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }
#ifdef IWAVE_USE_MPI
    if (retrieveGroupID() == MPI_UNDEFINED) {
      fprintf(stream,"NOTE: finalize MPI, cleanup, exit\n");
      fflush(stream);
    }
    else {
#endif

        string coeff   = valparse<string>(*pars,"coeff");
        string inmig   = valparse<string>(*pars,"inmig");
        string inremig = valparse<string>(*pars,"inremig");
        string outig   = valparse<string>(*pars,"outig");

        gsp sp(inmig, "csq", true
#ifdef IWAVE_USE_MPI 
           , retrieveGlobalComm()
#endif
           );
        gsp csp(coeff, "coeff", true
#ifdef IWAVE_USE_MPI 
           , retrieveGlobalComm()
#endif
           );

        cerr << "in PsiDO\n";
        Vector<float> incoeff(csp);
        AssignFilename afincoeff(coeff);
        incoeff.eval(afincoeff);
        /* PsiDO op */
        PsiDiffOp psido(sp,incoeff,*pars,stream,1);

        Vector<float> invec(sp);
        Vector<float> inrmvec(sp);
        Vector<float> outvec(sp);
        AssignFilename afin(inmig);
        AssignFilename afinrm(inremig);
        AssignFilename afout(outig);
        invec.eval(afin);
        inrmvec.eval(afinrm);
        outvec.eval(afout);
        outvec.copy(inrmvec);
        
        fprintf(stream, "inrmvec.norm()=%f\n", inrmvec.norm());
        fprintf(stream, "outvec.norm()=%f\n",outvec.norm());
#ifdef IWAVE_USE_MPI
        MPI_Barrier(retrieveGlobalComm());
#endif
        fprintf(stream, "outvec.norm()=%f\n",outvec.norm());

        psido.applyOp(invec,outvec);

        fprintf(stream, "outvec.norm()=%f\n",outvec.norm());
#ifdef IWAVE_USE_MPI
    }
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
