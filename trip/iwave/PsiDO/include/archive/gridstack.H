#ifndef __RVL_GRID_FLOAT_STACK
#define __RVL_GRID_FLOAT_STACK

#include "griddata.H"

namespace RVLGrid {

  using RVL::BinaryLocalFunctionObject;
  using RVL::RVLException;

  /** Extends input GridData d1 over d2 having larger-dimensional
      grid:
      
      d2[i_0,i_1,...,i_n1,i_{n1+1},...i_n2]
        = d1[i_0,i_1,...,i_n1]

      all i_0,...,i_n2. Used to treat lower-dimensional data as
      effectively higher dimensional.

      Preconditions: both d1 and d2 are actually GridData objects, and
      d2 has at least as many axes (n2) as d1 (n1). Checked by requiring
      (1) d2.size() >= d1.size(), (2) nd, d=1,...,n1 same for d1 and d2.

      Postconditions: data related as above.

      Note: aspects of grid geometry other than number of gridpoints
      per axis not checked for compatibility.
      
      Spread and Stack are adjoints */

  class Spread: public BinaryLocalFunctionObject<float> {

  public:

    Spread() {}
    Spread(const Spread &) {}
    virtual ~Spread() {}

    void operator()(LocalDataContainer<float> & d2,
		    LocalDataContainer<float> const & d1);

    string getName() const { string str="RVLGrid::Spread"; return str; }

  };

  /** Volume-averages input GridData d1, dimension n1, over last n1-n2
      axes and stores average in output grid d2 having dimension n2:
      
      d1[i_0,i_1,...,i_n1]
        = dx_n1*dx_{n1+1}*...*dx_n2 sum_{i_n1,...i_n2} 
          d2[i_0,i_1,...,i_n1,i_{n1+1},...i_n2]
       

      all i_0,...,i_n2. Adjoint to spread.

      Preconditions: both d1 and d2 are actually GridData objects, and
      d1 has at least as many axes (n1) as d2 (n2). Checked by requiring
      (1) d1.size() >= d2.size(), (2) nd, d=1,...,n1 same for d1 and d2.

      Postconditions: data related as above.

      Note: aspects of grid geometry other than number of gridpoints
      per axis not checked for compatibility.
  */
  
  class Stack:public BinaryLocalFunctionObject<float> {

  public:

    Stack() {}
    Stack(const Stack &) {}
    virtual ~Stack() {}

    void operator()(LocalDataContainer<float> & d2,
		    LocalDataContainer<float> const & d1);

    string getName() const { string str="RVLGrid::Stack"; return str; }

  };

}

#endif

	
	
