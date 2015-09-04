#include "std_cpp_includes.H"
#include "grid.H"

using namespace RVL;
using namespace RVLGrid;

int main() {
  try {
    GridSpace<float> sp("newtest.sep",3);
    sp.write(cout);
    int n1;
    sp.getAttr("n1",n1);
    cout<<"attribute n1 is "<<n1<<endl;
    int naxes;
    sp.getAttr("naxes",naxes);
    cout<<"attribute naxes is "<<naxes<<endl;
    cout<<"size is "<<sp.getSize()<<endl;
    float d2;
    sp.getAttr("d2",d2);
    cout<<"attribute d2 is "<<d2<<endl;
    sp.getAttr("flower",n1);
    cout<<"attribute flower is "<<n1<<endl;
  }
  catch (RVLException & e) {
    e.write(cerr);
  }
}
