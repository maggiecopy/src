
#include "grid.H"

using namespace RVL;

template<class Scalar>
class Laplace2: public BinaryFunctionObject<Scalar> {
private:
  vector<int> n;
  Scalar c0;
  vector<Scalar> c1;
  int bc; // 0 = Dirichlet, 1 = Neumann
  int gsize;
  Laplace2() {}
  Laplace2(const Laplace2<Scalar> &) {}
public:
  Laplace2(Grid<Scalar> & g, int _bc=0)
    : n(g.getNaxes()), c0(0.0), c1(g.getNaxes()), bc(_bc), gsize(g.getSize()) {
    for (int i=0;i<n.size();i++) {
      n[i]=g.get_n()[i];
      if (bc<0||bc>1) {
	RVLException e; e<<"Error: Laplace2 constructor\n";
	e<<"only Dirichlet (bc=0) and Neumann (bc=1) conditions implemented\n";
	throw e;
      }
      Scalar dr;
      if (ProtectedDivision<Scalar>(1.0,g.get_d()[i],dr)) {
	RVLException e; e<<"Error: Laplace2 constructor\n";
	e<<"step "<<i+1<<" of grid causes zerodivide\n";
	g.write(e);
	throw e;
      }
      c1[i]=dr*dr;
      c0  += -2*c1[i];
    }
  }
  ~Laplace2() {}
  virtual void operator()(LocalDataContainer<Scalar> & d1,
			  LocalDataContainer<Scalar> & d2) {
    if (gsize != d1.getSize() || gsize != d2.getSize()) { 
      RVLException e; e<<"Error: Laplace2::operator()\n";
      e<<"input LocalDataContainers must have size of defining grid\n";
      e<<"grid size = "<<gsize<<"\n";
      e<<"output data container:\n";
      d1.write(e);
      e<<"input data container:\n";
      d2.write(e);
      throw e;
    }
    Scalar * p1 = d1.getData();
    Scalar * p2 = d2.getData();
    if (n.size()==1) {

      int n0=n[0];
      Scalar co1 = c1[0];
      for (int i=1;i<n0-1;i++) {
	p1[i]=co1*(p2[i-1]+p2[i+1])+c0*p2[i];
      }
      if (bc == 0) {
	p1[0]=0.0;
	p1[n0-1]=0.0;
      }
      else if (bc == 1) {
	p1[0]=2*co1*p2[1]+c0*p2[0];
	p1[n0-1]=2*co1*p2[n0-2]+c0*p2[n0-1];
      }
    }
    else if (n.size()==2) {
      int n0=n[0];
      int n1=n[1];
      Scalar co1 = c1[0];
      Scalar co2 = c1[1];
      for (int i=1;i<n0-1;i++) {
	for (int j=1;i<n1-1;j++) {
	  p1[j*n0+i]=co1*(p2[j*n0+i-1]+p2[j*n0+i+1])+
	    co2*(p2[(j-1)*n0+i]+p2[(j+1)*n0+i])+c0*p2[j*n0+i];
	}
      }
      if (bc == 0) {
	for (int i=0;i<n0;i++) {
	  p1[i]=0.0;
	  p1[(n1-1)*n0+i]=0.0;
	}
	for (int j=0;j<n1;j++) {
	  p1[j*n0]=0.0;
	  p1[(j+1)*n0-1]=0.0;
	}
      }
      else if (bc == 1) {
	// sides
	for (int i=1;i<n0-1;i++) {
	  p1[i]=co1*(p2[i-1]+p2[i+1])+
	    2*co2*p2[n0+i] +c0*p2[i];
	  p1[(n1-1)*n0+i]=co1*(p2[(n1-1)*n0+i-1]+p2[(n1-1)*n0+i+1])+
	    2*co2*p2[(n1-2)*n0+i] + c0*p2[(n1-1)*n0+i];
	}
	for (int j=1;j<n1-1;j++) {
	  p1[j*n0]=2*co1*p2[j*n0+1]+
	    co2*(p2[(j+1)*n0]+p2[(j-1)*n0]) + c0*p2[j*n0];
	  p1[(j+1)*n0-1]=2*co1*p2[(j+1)*n0-2]+
	    co2*(p2[j*n0-1]+p2[(j+2)*n0-1]) + c0*p2[(j+1)*n0-1];
	}
	// corners
	p1[0]=2*co1*p2[1]+2*co2*p2[n0]+c0*p2[0];
	p1[n1-1]=2*co1*p2[n0-2]+2*co2*p2[2*n0-1]+c0*p2[n0-1];
	p1[(n1-1)*n0]=2*co1*p2[(n1-1)*n0+1]+
	  2*co2*p2[(n1-2)*n0]+c0*p2[(n1-1)*n0];
	p1[n0*n1-1]=2*co1*p2[n0*n1-2]+
	  2*co2*p2[(n1-1)*n0-1]+c0*p2[n0*n1-1];
      }
    }
    else {
      RVLException e; e<<"Error: Laplace2::operator()\n";
      e<<"dimension "<<(int)(n.size())<<"not implemented\n";
      throw e;
    }
  }

  string name() { return "Laplace2"; }

};

int main() {
  try {
    GridSpace<float> sp("spike1d.sep");
    Vector<float> vin(sp);
    GridLoad<float> ldr("spike1d.sep");
    vin.eval(ldr);
    cout<<"INPUT:"<<endl;
    vin.write(cout);
    Laplace2<float> lap(sp.getGrid());
    Vector<float> vout(sp);
    vout.eval(lap,vin);
    cout<<"OUTPUT:"<<endl;
    vout.write(cout);
  }
  catch (RVLException & e) {
    e.write(cerr);
  }
}
