#ifdef DOUBLE_PRECISION
    typedef double data_t;
#else
    typedef float data_t;
#endif

#include <string.h>
#include <unistd.h>
#include "operator.hpp"
#include "bsplines.hpp"
#include "param.hpp"
#include "IO.hpp"
#include "seplib.h"

typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;


void printdoc(){
    std::string doc = "\nDescription:\n"
    "   Apply cubic B-splines smoothing.\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'.\n"
    "\nParameters:\n"
    "   nx,nz - int :\n\t\tnumber of control points in each direction. If provided, will override the parameters below.\n"
    "   controlx,controlz - [float] :\n\t\tarrays of control points manually entered. Must contain the first and last physical points. To be used in conjunction with mx,mz.\n"
    "   mx,mz - [int] :\n\t\tarrays of control points multiplicity manually entered.\n"
    "\nNote:\n"
    "   The number of control points must be >= 3 in all cases.\n"
    "\nExamples:\n"
    "   BSPLINES.x < infile.H nx=11 nz=23 > oufile.H.\n"
    "   BSPLINES.x < infile.H controlx=0,1,5.5,10 controlz=0,5,20 mx=2,1,1,2 mz=2,1,2 > oufile.H.\n"
    "   BSPLINES.x < infile.H nx=11 controlz=0,5,20 mz=2,1,2 > oufile.H.\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);

    std::string input_file="in", output_file="out";
    int nx=1, nz=1;
    std::vector<data_t> controlx={0}, controlz={0}; // locations of the control points that define the B-splines
    std::vector<int> mx={0}, mz={0}; // multiplicity of the control points that define the B-splines

	readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
    readParam<int>(argc, argv, "nx", nx);
    readParam<int>(argc, argv, "nz", nz);
    readParam<int>(argc, argv, "mx", mx);
    readParam<int>(argc, argv, "mz", mz);
    readParam<data_t>(argc, argv, "controlx", controlx);
    readParam<data_t>(argc, argv, "controlz", controlz);
    
    std::shared_ptr<vec> input = sepRead<data_t>(input_file);
    std::shared_ptr<vec> output = std::make_shared<vec>(*input->getHyper());
    const data_t* pinput = input->getCVals();
    data_t* poutput = output->getVals();

    // Build the knots, multiplicity, and the B-splines smoothing operator
// ----------------------------------------------------------------------------------------//
    if (nx > 1){
        axis<data_t> X = input->getHyper()->getAxis(2);
        controlx.resize(nx); mx.resize(nx);
        controlx[0] = X.o;
        data_t dx = X.d*(X.n-1)/(nx-1);
        for (int i=1; i<nx-1; i++) {
            controlx[i] = controlx[0] + i*dx;
            mx[i] = 1;
        }
        controlx[nx-1] = X.o + (X.n-1)*X.d;
        mx[0] = 2;
        mx[nx-1] = 2;
    }
    else successCheck(mx.size() == controlx.size(),__FILE__,__LINE__,"Multiplicity and control vectors must be of the same size\n");
    if (nz > 1){
        axis<data_t> Z = input->getHyper()->getAxis(1);
        controlz.resize(nz); mz.resize(nz);
        controlz[0] = Z.o;
        data_t dz = Z.d*(Z.n-1)/(nz-1);
        for (int i=1; i<nz-1; i++) {
            controlz[i] = controlz[0] + i*dz;
            mz[i] = 1;
        }
        controlz[nz-1] = Z.o + (Z.n-1)*Z.d;
        mz[0] = 2;
        mz[nz-1] = 2;
    }
    else successCheck(mz.size() == controlz.size(),__FILE__,__LINE__,"Multiplicity and control vectors must be of the same size\n");
    
    std::vector<data_t> kx;
    std::vector<data_t> kz;
    setKnot(kx,controlx,mx);
    setKnot(kz,controlz,mz);

    std::vector<axis<data_t> > axes = input->getHyper()->getAxes();
    axes[0].n=controlz.size(); axes[1].n=controlx.size();
    std::shared_ptr<vecReg<data_t> > bsmodel (new vecReg<data_t>(hypercube<data_t>(axes)));
    fillin(bsmodel,input,controlx,controlz);
  
    duplicate D(*bsmodel->getHyper(),mx,mz);
    bsplines3 B(*D.getRange(),*input->getHyper(),kx,kz);
    chainLOper * BD = new chainLOper(&B,&D);
    BD->forward(false,bsmodel,output);

    delete BD;
// ----------------------------------------------------------------------------------------//

    sepWrite<data_t>(output,output_file);

    return 0;
}



