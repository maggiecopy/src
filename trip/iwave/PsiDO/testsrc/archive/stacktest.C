#include "stackop.H"

using namespace RVL;
using namespace RVLGrid;

int main() {
 
  try {

    srand(getpid());

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

    SpreadOp op(dom,1);

    RVLRandomize<float> rnd;

    op.checkAdjointRelation(rnd,cout);

  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
