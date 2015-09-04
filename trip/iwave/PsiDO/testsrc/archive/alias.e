#include "functions.H"
#include "cubicfo.H"

using namespace RVL;
using namespace RVLGrid;

int main() {

  try{
  Grid<float> gin("fine.sep");
  Grid<float> gout("outgrid.sep");
  GridData<float> din(gin);
  GridData<float> dout(gout);

  GridLoad<float> dload("fine.sep");
  dload(din);

  CubicFwdInterp fwd;
  CubicAdjInterp adj;
  
  fwd(dout,din);

  GridSave<float> dsave("fwdinterp.sep");
  dsave(dout);

  GridData<float> din1(gin);
  
  adj(din1,dout);

  GridSave<float> dsave1("adjinterp.sep");
  dsave1(din1);

  float dtin=gin.get_d(0);
  float dtout=gout.get_d(0);

  RVLL2innerProd<float> ipin;
  RVLL2innerProd<float> ipout;
  ipout(dout,dout);
  cout<<"inner product in range: "<<dtout*ipout.getResult()<<endl;
  ipin(din,din1);
  cout<<"inner product in domain: "<<dtin*ipin.getResult()<<endl;
  }
  catch (RVLException & e) {
    e.write(cerr);
  }
}




