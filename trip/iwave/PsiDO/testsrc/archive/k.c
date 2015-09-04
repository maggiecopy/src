#include "std_cpp_includes.H"
#include "grid.H"

using namespace RVL;
using namespace RVLGrid;

int main() {
  try {
    GridSpace<double> sp("test1.sep",3);
    Vector<double> vec(sp);
    RVLRandomize<double> f;
    vec.eval(f);
    GridSave<double> svr("test.sep");
    vec.eval(svr);
    vec.write(cout);
    Vector<double> vec0(sp);
    vec0.copy(vec);
    Vector<double> vec1(sp);
    GridLoad<double> ldr("test.sep");
    vec1.eval(ldr);
    cout<<"before diff"<<endl;
    vec1.write(cout);
    vec1.linComb(1,vec1,-1,vec0);
    cout<<"after diff"<<endl;
    vec1.write(cout);
    cout<<"mean square difference = "<<vec1.inner(vec1)<<" should = 0"<<endl;
  }
  catch (RVLException &e) {
    e.write(cerr);
  }
}
