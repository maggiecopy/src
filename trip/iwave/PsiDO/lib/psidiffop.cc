#include "psidiffop.hh"

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n1,
    const lbfgsfloatval_t step
    ){
    PsiDO_PARS *psido_pars = (PsiDO_PARS *)instance;
    int len = psido_pars->nc[0]
             *psido_pars->nc[1]
             *(2*psido_pars->k-1);
    int m;

    lbfgsfloatval_t fx = 0.0;

    for(m=0;m<len;m++)
    {
      psido_pars->coeff_c[m]=x[m];
    }

    //fprintf(psido_pars->stream, "before fx = %f\n", fx);
    objective_real(psido_pars->U,     psido_pars->QUtest, psido_pars->coeff_c,
                   psido_pars->nc[1], psido_pars->nc[0],  psido_pars->k,
                   psido_pars->o[1],  psido_pars->n[1],   psido_pars->d[1],
                   psido_pars->o[0],  psido_pars->n[0],   psido_pars->d[0],
                   psido_pars->m_ord, psido_pars->Obj,    psido_pars->Grad,
                   1, psido_pars->lam);

    fx=psido_pars->Obj[0];
    //fprintf(psido_pars->stream, "after fx = %f\n", fx);

    for(m=0;m<len;m++)
    {
      g[m]=psido_pars->Grad[m];
    }

  return fx;
}


static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n1,
    int k1,
    int ls
    )
{
    //fprintf(psido_pars->stream, "after fx = %f\n", fx);
      PsiDO_PARS *psido_pars = (PsiDO_PARS *)instance;
      fprintf(psido_pars->stream,"Iteration %d:\n", k1);
      fprintf(psido_pars->stream,"  fx = %f, x[0] = %f, x[center] = %f\n", fx, x[0], x[(psido_pars->nc[1]+1)*psido_pars->nc[0]/2+(psido_pars->nc[0]+1)/2]);
      fprintf(psido_pars->stream,"  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
      fprintf(psido_pars->stream,"\n");
      return 0;
}



namespace TSOpt {

#ifdef IWAVE_USE_MPI
  typedef MPIGridSpace myGridSpace;
#else
  typedef GridSpace myGridSpace;
#endif

  PsiDiffOp::PsiDiffOp(Space<ireal> const & _dom,
                       Vector<ireal> & _coeff,
                       PARARRAY _pars, FILE * _stream,
                       int _flag,
                       string _datatype,
                       string _cdatatype,
                       ostream & _announce)
    :dom(_dom), coeff(_coeff), flag(_flag), 
     stream(_stream), pars(NULL), announce(_announce),
     datatype(_datatype),cdatatype(_cdatatype),
     dump_steps(0), dump_pars(0), dump_term(0) {
  try{ 
    if(flag==0){
      cerr << "PsiDiffOp::PsiDiffOp  WARNing!!!!!!\n";
      cerr << "    default value of flag=0,\n";
      cerr << "    input coefficient is used to compute the PsiDO.\n";
      cerr << "    give flag=1 or any nonzero ints to trigger the inversion.\n";
    }
    int err; 
    // copy input par array to data member 
    if ((err=ps_copy(&pars,_pars))) {
      RVLException e;
      e<<"Error: PsiDiffOp constructor from ps_copy, err="<<err<<"\n";
      throw e;
    }
    // set dump controls
    ps_flint(*pars,"dump_steps",&dump_steps);
    ps_flint(*pars,"dump_pars",&dump_pars);
    ps_flint(*pars,"dump_term",&dump_term);

    // see what we've got
    if (dump_pars) {
      fprintf(stream,"PARS IN PSIDIFFOP CONSTRUCTOR\n");
      ps_printall(*pars,stream);
    }

    ProductSpace<ireal> const * pdom = dynamic_cast<ProductSpace<ireal> const *>(&dom);
    if (pdom){
      myGridSpace const & gdom = dynamic_cast<myGridSpace const &>((*pdom)[0]);
      if (!(&gdom)) {
        RVLException e;
        e<<"Error: PsiDiffOp::PsiDiffOp\n";
        e<<"  input space is not a GridSpace\n";	
        e<<"  description:\n";
        dom.write(e);	
        throw e;
      }
#ifdef IWAVE_USE_MPI
      if (retrieveGlobalRank() == 0) {
#endif
        g  = gdom.getGrid(); 
#ifdef IWAVE_USE_MPI
      }
#endif
    }
    else {
      myGridSpace const & gdom = dynamic_cast<myGridSpace const &> (dom);
#ifdef IWAVE_USE_MPI
      if (retrieveGlobalRank() == 0) {
#endif
        g  = gdom.getGrid(); 
#ifdef IWAVE_USE_MPI
      }
#endif
    }
    myGridSpace const * gcdom = dynamic_cast<myGridSpace const *>(&(coeff.getSpace()));
    if (!(&gcdom)) {
      RVLException e;
      e<<"Error: PsiDiffOp::PsiDiffOp\n";
      e<<"  input coefficient vector is not in a GridSpace\n";	
      e<<"  description:\n";
      coeff.write(e);	
      throw e;
    }
#ifdef IWAVE_USE_MPI
    if (retrieveGlobalRank() == 0) {
#endif
      gc = gcdom->getGrid();
#ifdef IWAVE_USE_MPI
    }
    if (MPI_Bcast(&g, sizeof(g), MPI_BYTE,0,retrieveGlobalComm()) || 
        MPI_Bcast(&gc,sizeof(gc),MPI_BYTE,0,retrieveGlobalComm()) ) {
      RVLException e;
      e<<"Error: PsiDiffOp, rank="<<retrieveGlobalRank()<<"\n";
      e<<"  failed to bcast grid\n";
      throw e;
    }
#endif
    }
    catch (bad_cast) {
      RVLException e;
      e<<"\nError: PsiDiffOp::PsiDiffOp\n";
      e<<"  either domain or range is neither product nor a GridSpace,\n";
      e<<"  or some component is not a GridSpace\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from PsiDiffOp::PsiDiffOp\n";
      throw e;
    }
  }

  PsiDiffOp::PsiDiffOp(PsiDiffOp const & A)
   :dom(A.dom), coeff(A.coeff), flag(A.flag), stream(A.stream), pars(NULL),
    announce(A.announce), datatype(A.datatype), cdatatype(A.cdatatype),
    dump_steps(A.dump_steps), dump_pars(A.dump_pars),
    dump_term(A.dump_term){
  
    int err;
    // copy input par array to data member 
    if ((err=ps_copy(&pars,*(A.pars)))) {
      RVLException e;
      e<<"Error: PsiDiffOp constructor from ps_copy, err="<<err<<"\n";
      throw e;
    }
  
    // see what we've got
    if (dump_pars) {
      fprintf(stream,"PARS IN PSIDIFFOP CONSTRUCTOR\n");
      ps_printall(*pars,stream);
    }

  }

  void PsiDiffOp::apply(const Vector<ireal> & x,
                        Vector<ireal> & y) const {
    try {
      SpaceTest(this->getDomain(),x,"TSOpt::PsiDiffOp::apply (dom)");
      SpaceTest(this->getRange(),y,"TSOpt::PsiDiffOp::apply (rng)");

      PARARRAY * locpars = ps_new();
      PARARRAY * xpars = ps_new();
      PARARRAY * ypars = ps_new();
      PARARRAY * cpars = ps_new();
      ps_copy(&locpars,*pars);

      AssignParams xap(*xpars,datatype,stream);
      AssignParams yap(*ypars,datatype,stream);
      AssignParams cap(*cpars,cdatatype,stream);
      x.eval(xap);  y.eval(yap); coeff.eval(cap); 

      if (dump_pars) {
        fprintf(stream,"PARAM ARRAY CREATED IN PSIDIFFOP::APPLY\n");
        ps_printall(*locpars,stream);
        fflush(stream);
      }


      if (dump_steps) {
	  fprintf(stream,"PsiDiffOp::APPLY -> STEP\n");
	  fflush(stream);
      }

      string xname, yname, cname;
      parse(*xpars,datatype,xname);
      parse(*ypars,datatype,yname);
      parse(*cpars,cdatatype,cname);
      if (dump_term) {
        fprintf(stream, "input vector x file name = %s\n", xname.c_str());
        fprintf(stream, "input vector y file name = %s\n", yname.c_str());
        fprintf(stream, "input vector coeff file name = %s\n", cname.c_str());
      }

      IPNT gsc_arr;
      IPNT nc_arr;
      IPNT gs_arr;
      IPNT n_arr;
      RPNT d_arr;
      RPNT o_arr;
      if (g.dim != 2) {
        RVLException e;
        e<<"Error: PsiDiffOp::apply\n";
        e<<"  current implementation is 2D only\n";
        throw e;
      }

      get_gs(gsc_arr,gc);
      get_n(nc_arr,gc);
      get_gs(gs_arr,g);
      get_d(d_arr,g);
      get_n(n_arr,g);
      get_o(o_arr,g);
      if (dump_term) {
         for (int i=0; i<g.dim; i++){
           fprintf(stream, "  n_arr[%d]=%d\n ", i, n_arr[i]);
           fprintf(stream, "  d_arr[%d]=%f\n ", i, d_arr[i]);
           fprintf(stream, "  o_arr[%d]=%f\n ", i, o_arr[i]);
         }
         for (int i=0; i<gc.dim; i++){
           fprintf(stream, "  nc_arr[%d]=%d\n ", i, nc_arr[i]);
         }
      }  

      int first;
      int last;
      int panelnum=1;
      int panelindex=0;
      // compute istart
      IPNT istart;
      IPNT istop;
      get_gs(istart,g);
      get_ge(istop,g);

      // nontrivial extended axes
      if (g.gdim > g.dim) {
        for (int i=g.dim; i<g.gdim; i++) {
          panelnum *= istop[i]-istart[i]+1;
        }
      }
      else {
        panelnum=1;
        panelindex=0;
      }
      calc_group(&first, &last, panelnum);
      for (panelindex=first; panelindex<=last; panelindex++) { 
        announce<<"    panel="<<panelindex<<" in range=["<<first<<", "<<last<<"], rank="<<retrieveGlobalRank()<<endl;

        int lenc=nc_arr[0]*nc_arr[1]*nc_arr[2];
        int len = n_arr[0]* n_arr[1];
        PsiDO_PARS *psido_pars;
        psido_pars = (PsiDO_PARS *)usermalloc_(sizeof(PsiDO_PARS));
        psido_pars->nc[0]  = nc_arr[0];
        psido_pars->nc[1]  = nc_arr[1];
        psido_pars->n[0]   = n_arr[0];
        psido_pars->n[1]   = n_arr[1];
        psido_pars->d[0]   = d_arr[0];
        psido_pars->d[1]   = d_arr[1];
        psido_pars->o[0]   = o_arr[0];
        psido_pars->o[1]   = o_arr[1];
        psido_pars->stream = stream;

        psido_pars->U      = (fftwf_complex *)fftwf_malloc(len*sizeof(fftwf_complex));
        psido_pars->QUtest = (fftwf_complex *)fftwf_malloc(len*sizeof(fftwf_complex));
    
        psido_pars->coeff_c = (ireal *)malloc(lenc*sizeof(ireal));
        psido_pars->Grad    = (ireal *)malloc(lenc*sizeof(ireal));
        psido_pars->Obj     = (ireal *)malloc(sizeof(ireal));
#pragma ivdep
        for (int m=0; m<lenc; m++){
          psido_pars->coeff_c[m]=0.f;
          psido_pars->Grad[m]=0.f;
        }
        psido_pars->Obj[0]=0.f;
    
        psido_pars->m_ord = valparse<ireal>(*pars,"m_ord",-1);
        psido_pars->k     = (nc_arr[2]+1)/2;
        psido_pars->lam   = 0.f;
    
        int ret = 0;
        lbfgsfloatval_t fx;
        lbfgsfloatval_t *x = lbfgs_malloc(lenc);
        lbfgs_parameter_t param;

        if (x == NULL) {
          RVLException e;
          e<<"Error: PsiDiffOp::apply\n";
          e<<"  Failed to allocate a memory block for variables of lbfgs.\n";
          throw e;
        }

        float * Ureal=(float *)malloc(len*sizeof(float));
#pragma ivdep
        for (int m=0; m<len; m++){
          Ureal[m]=0.f;
        }
//        int err=rsfread(Ureal,gs_arr,n_arr,"../rmcsqext.rsf",1.0,stream,panelindex);
        int err=rsfread(Ureal,gs_arr,n_arr,yname.c_str(),1.0,stream,panelindex);
        fprintf(stream, "  panel=%d, filename=%s\n ", panelindex, yname.c_str());
	if (err) {
	  RVLException e;
	  e<<"Error: Error: PsiDiffOp::apply, rsfread input re-migrated image\n";
	  e<<"  failed to read panel "<<panelindex<<" from "<<yname<<"\n";
	  throw e;
	}
        float * QUreal=(float *)malloc(len*sizeof(float));
#pragma ivdep
        for (int m=0; m<len; m++){
          QUreal[m]=0.f;
        }
        err=rsfread(QUreal,gs_arr,n_arr,xname.c_str(),1.0,stream,panelindex);
	if (err) {
	  RVLException e;
	  e<<"Error: Error: PsiDiffOp::apply, rsfread input migrated image\n";
	  e<<"  failed to read panel "<<panelindex<<" from "<<xname<<"\n";
	  throw e;
	}

        //cerr << "QUreal[825400] = " << QUreal[825400] << endl;
        //cerr << " Ureal[825400] = " <<  Ureal[825400] << endl;
        float xr,zr;
        float max_U =0;
        float n_QU=0;

        for (int m=0;m<n_arr[1];m++){
          for (int n=0;n<n_arr[0];n++){
            xr=o_arr[1]+m*d_arr[1];
            zr=o_arr[0]+n*d_arr[0];
               
            psido_pars->U[m*n_arr[0]+n][0]=Ureal[m*n_arr[0]+n];
            psido_pars->U[m*n_arr[0]+n][1]=0.f; ;
            if(Ureal[m*n_arr[0]+n] > max_U)
              max_U = Ureal[m*n_arr[0]+n];

            psido_pars->QUtest[m*n_arr[0]+n][0]=QUreal[m*n_arr[0]+n];
            psido_pars->QUtest[m*n_arr[0]+n][1]=0.f;
            n_QU = n_QU + (psido_pars->QUtest[m*n_arr[0]+n][0]*psido_pars->QUtest[m*n_arr[0]+n][0] 
                        +  psido_pars->QUtest[m*n_arr[0]+n][1]*psido_pars->QUtest[m*n_arr[0]+n][1]);
          }
        }

        n_QU=pow(n_QU,0.5);
        fprintf(stream, "  panel=%d, n_QU=%f\n ", panelindex, n_QU);
        fprintf(stream, "  panel=%d, max_U=%f\n ", panelindex, max_U);
        //cerr << " n_QU = " << n_QU << endl;
#pragma ivdep
        for(int m=0;m<len;m++){
            psido_pars->U[m][0]     =psido_pars->U[m][0]/n_QU;
            psido_pars->U[m][1]     =psido_pars->U[m][1]/n_QU;
            psido_pars->QUtest[m][0]=psido_pars->QUtest[m][0]/n_QU;
            psido_pars->QUtest[m][1]=psido_pars->QUtest[m][1]/n_QU;
        }
        // read in supplied coefficient if flag==0
        if (flag==0){
          float * coeff_prev=(float *)malloc(lenc*sizeof(float));
          err=rsfread(coeff_prev,gsc_arr,nc_arr,cname.c_str(),1.0,stream,panelindex);
	  if (err) {
	    RVLException e;
	    e<<"Error: Error: PsiDiffOp::apply, rsfread input coefficient of PsiDO\n";
	    e<<"  failed to read panel "<<panelindex<<" from "<<xname<<"\n";
	    throw e;
	  }
#pragma ivdep
          for(int m=0;m<lenc;m++){
            x[m]=coeff_prev[m];
          }
          free(coeff_prev);
        }
        // initialize coefficient to 0 otherwise
        // and run lbfgs inversion to compute coefficients
        else {
#pragma ivdep
          for(int m=0;m<lenc;m++)
            x[m]=0;

          //========================================================//
          // Initialize the parameters for the L-BFGS optimization. //
          //========================================================//
          lbfgs_parameter_init(&param);
          //param.orthantwise_c = 1;
          //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE
          param.max_iterations=valparse<int>(*pars,"PsiDOmaxiters",3);
          param.epsilon=valparse<ireal>(*pars,"PsiDOepsilon",1e-8);
          param.xtol=valparse<ireal>(*pars,"PsiDOxtol",1e-8);
          param.m=valparse<int>(*pars,"PsiDOm",10);
          // param.gtol=0.1;
          param.linesearch = 0;
          // param.delta=100;
       
          //========================================================//
          // Start the L-BFGS optimization; this will invoke the callback functions
          //  evaluate() and progress() when necessary.
          //========================================================//
          ret = lbfgs(lenc, x, &fx, evaluate, progress, psido_pars, &param);
 
          // Report the result.
          fprintf(stream, " L-BFGS optimization terminated with status code = %d\n", ret);
          fprintf(stream, " fx = %f, x[0] = %f, x[center] = %f\n", fx, x[0], x[(nc_arr[1]+1)*nc_arr[0]/2+(nc_arr[0]+1)/2]);
        }
        // array to store inverted image
        fftwf_complex  *inv = (fftwf_complex *)fftwf_malloc(len*sizeof(fftwf_complex));
#pragma ivdep
        for(int m=0;m<lenc;m++){
          psido_pars->coeff_c[m]=x[m];
        }
        if (flag!=0){
          // write PsiDO scale into file
          err=rsfwrite(psido_pars->coeff_c,gsc_arr,nc_arr,cname.c_str(),1.0,stream,panelindex);
          if (err) {
            RVLException e;
            e<<"Error: Error: PsiDiffOp::apply, rsfwrite output coefficient of PsiDO\n";
            e<<"  failed to read panel "<<panelindex<<" from "<<cname<<"\n";
            throw e;
          }
        }
        make_inv_real(psido_pars->QUtest,inv,psido_pars->coeff_c,nc_arr[1],nc_arr[0],psido_pars->k,o_arr[1],n_arr[1],d_arr[1],o_arr[0],n_arr[0],d_arr[0],psido_pars->m_ord);

        float * inv_real=(float *)malloc(len*sizeof(float));
#pragma ivdep
        for(int m=0;m<len;m++){
          inv_real[m]=inv[m][0]*n_QU;
        }

        //cerr << "inv_real[825400] = " << inv_real[825400] << endl;
        err=rsfwrite(inv_real,gs_arr,n_arr,yname.c_str(),1.0,stream,panelindex);
        if (err) {
          RVLException e;
          e<<"Error: Error: PsiDiffOp::apply, rsfwrite output scaled image\n";
          e<<"  failed to read panel "<<panelindex<<" from "<<xname<<"\n";
          throw e;
        } 

        lbfgs_free(x);
        free(Ureal);
        free(QUreal);
        free(inv);
        free(inv_real);
        free(psido_pars->coeff_c);
        free(psido_pars->Grad);
        free(psido_pars->Obj);
        fftwf_free(psido_pars->U);
        fftwf_free(psido_pars->QUtest);

        userfree_(psido_pars);


      } // end of for panelindex loop


#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from TSOpt::PsiDiffOp::apply\n";
      throw e;
    }


  }

  void PsiDiffOp::applyAdj(const Vector<ireal> & x,
                           Vector<ireal> & dy) const {
    try{
      RVLException e;
      e<<"Error: PsiDiffOp::applyAdj\n";
      e<<"  adjoint method has not been implemented\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from TSOpt::PsiDiffOp::applyAdj\n";
      throw e;
    }
  }



}
