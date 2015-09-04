#include "table.H"
#include "stackop.H"

using namespace RVL;
using namespace RVLGrid;

int main(int argc, char ** argv) {
 
  try {

    if (argc != 2) {
      cerr<<"Spread.x: create nD GridData object from\n";
      cerr<<"mD GridData object, m <= n\n";
      cerr<<"usage: Spread.x <parfile name>\n";
      cerr<<"required parameters, identified in key=value form in parfile:\n";
      cerr<<"  InFile = <name of input data file> (SEP header file)\n";
      cerr<<"  TargetGrid = <name of file specifying target grid> (SEP header file)\n";
      cerr<<"  OutFile = <name of output file>\n";
      cerr<<"InFile and TargetGrid must specify existing files. Infile\n";
      cerr<<"must specify data, whereas TargetGrid may be header only\n";
      exit(1);
    }

    Table par(argv[1]);

    string inname;
    if (par.getValue("InFile",inname)) {
      RVLException e;
      e<<"Error: Spread.x\n";
      e<<"InFile not specified in parameter table\n";
      throw e;
    }
    Grid<float> gin(inname);
    GridLoad<float> gld(inname.c_str());
    GridSpace<float> ginsp(gin);
    Vector<float> ginv(ginsp);
    ginv.eval(gld);

    string tgtname;
    if (par.getValue("TargetGrid",tgtname)) {
      RVLException e;
      e<<"Error: Spread.x\n";
      e<<"TargetGrid not specified in parameter table\n";
      throw e;
    }
    Grid<float> gtgt(tgtname);
    
    int defect = gtgt.getNaxes() - gin.getNaxes();
    if (defect < 0) {
      RVLException e;
      e<<"Error: Spread.x\n";
      e<<"TargetGrid has fewer axes than input grid\n";
      throw e;
    }      

    GridSpace<float> goutsp(gtgt);
    Vector<float> goutv(goutsp);

    SpreadOp op(goutsp,defect);
    
    op.apply(ginv,goutv);

    string outname;
    if (par.getValue("OutFile",outname)) {
      RVLException e;
      e<<"Error: Spread.x\n";
      e<<"OutFile not specified in parameter table\n";
      throw e;
    }    
    GridSave<float> gsv(outname);
    goutv.eval(gsv);
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
