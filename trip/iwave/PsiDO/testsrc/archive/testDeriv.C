#include "deriv.H"

using namespace RVL;
using namespace RVLGrid;

int main() {

  try {

    srand(getpid());

    vector<int> n(2);
    vector<float> d(2);
    vector<float> o(2);  

    n[0]=21;
    n[1]=41;
    d[0]=5;
    d[1]=5;
    o[0]=0.0;
    o[1]=0.0;

    Grid<float> g(n,d,o);
    int indx = 1;

    DerivOp op(g,indx);
  
    RVLRandomize<float> rnd;

    op.checkAdjointRelation(rnd,cout);

  }
  catch (RVLException & e) {
    e.write(cout);
    exit(1);
  }
}

