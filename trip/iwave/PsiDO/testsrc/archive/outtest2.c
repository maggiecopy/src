#include "linop.H"
#include "cubicfo.H"

using namespace RVL;
using namespace RVLGrid;

int main() {

  vector<int> n(2);
  vector<float> d(2);
  vector<float> o(2);  

  n[0]=5;
  n[1]=4;
  d[0]=5;
  d[1]=5;
  o[0]=0.0;
  o[1]=0.0;

  Grid<float> gin(n,d,o);
  GridSpace<float> dom(gin);

  //  n[0]=240;
  //  n[1]=480;
  //  d[0]=4.16667;
  //  d[1]=4.16667;
  n[0]=7;
  n[1]=5;
  d[0]=4.0;
  d[1]=4.0;
  o[0]=0.0;
  o[1]=0.0;

  Grid<float> gout(n,d,o);
  GridSpace<float> rng(gout);

  BiCubicFwdInterp<float> fwd;
  BiCubicAdjInterp<float> adj;
  
  LinearOpFO<float> op(dom,rng,fwd,adj);

  LocalVector<float> v1(dom);
  LocalVector<float> v2(dom);
  LocalVector<float> w1(rng);
  LocalVector<float> w2(rng);

  v1.zero();
  v2.zero();
  w1.zero();
  w2.zero();

  v1.getData()[3] = 1.0;
  w2.getData()[3] = 1.0;

  op.apply(v1,w1);
  op.applyAdj(w2,v2);

  cout << "forward input:"<<endl;
  v1.write(cout);
  cout << "forward output:"<<endl;
  w1.write(cout);
  cout << "adjoint input:"<<endl;
  w2.write(cout);
  cout << "adjoint output:"<<endl;
  v2.write(cout);
  float ip1 = w2.inner(w1);
  float ip2 = v2.inner(v1);
  cout << "<forward output, adjoint input> = "<<ip1<<endl;
  cout << "<adjoint output, forward input> = "<<ip2<<endl;

}
