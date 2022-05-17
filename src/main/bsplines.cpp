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
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'. '<' and '>' are valid for SEPlib format only.\n"
    "   bsmodel - string - ['none'] : optional B-splines model to map back to the input space. Must be compatible with the input and parameters below.\n"
    "\nParameters:\n"
    "   nx,nz - int - [3] :\n\t\tnumber of control points in each direction. If provided, will override the parameters below.\n"
    "   controlx,controlz - [float] :\n\t\tarrays of control points manually entered. Must contain the first and last physical points. To be used in conjunction with mx,mz.\n"
    "   mx,mz - [int] :\n\t\tarrays of control points multiplicity manually entered.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries.\n"
    "\nNote:\n"
    "   The number of control points must be >= 3 in all cases, unless nx=0 in which case the smoothing will be applied in the z-direction only.\n"
    "\nExamples:\n"
    "   BSPLINES.x < infile.H nx=11 nz=23 > oufile.H.\n"
    "   BSPLINES.x < infile.H controlx=0,1,5.5,10 controlz=0,5,20 mx=2,1,1,2 mz=2,1,2 > oufile.H\n"
    "   BSPLINES.x < infile.H nx=11 controlz=0,5,20 mz=2,1,2 > oufile.H\n"
    "   BSPLINES.x < infile.H bsmodel=bsm.H nx=11 nz=23 > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);
    omp_set_num_threads(1);

    std::string input_file="in", bsmodel_file="none", output_file="out", datapath="none";
    int nx=3, nz=3;
    bool format=0;
    std::vector<data_t> controlx={0}, controlz={0}; // locations of the control points that define the B-splines
    std::vector<int> mx={0}, mz={0}; // multiplicity of the control points that define the B-splines

	readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "bsmodel", bsmodel_file);
	readParam<std::string>(argc, argv, "output", output_file);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<int>(argc, argv, "nx", nx);
    readParam<int>(argc, argv, "nz", nz);
    readParam<int>(argc, argv, "mx", mx);
    readParam<int>(argc, argv, "mz", mz);
    readParam<bool>(argc, argv, "format", format);
    readParam<data_t>(argc, argv, "controlx", controlx);
    readParam<data_t>(argc, argv, "controlz", controlz);
    
    std::shared_ptr<vec> input = read<data_t>(input_file,format);
    std::shared_ptr<vec> output = std::make_shared<vec>(*input->getHyper());
    const data_t* pinput = input->getCVals();
    data_t* poutput = output->getVals();

    // Build the knots, multiplicity, and the B-splines smoothing operator
// ----------------------------------------------------------------------------------------//
    if (nx > 2){
        ax X = input->getHyper()->getAxis(2);
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
    if (nz > 2){
        ax Z = input->getHyper()->getAxis(1);
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

    std::vector<ax> axes = input->getHyper()->getAxes();
    axes[0].n=controlz.size(); 
    if (nx!=0) axes[1].n=controlx.size();
    
    std::shared_ptr<vec> bsmodel;
    if (bsmodel_file=="none"){
        bsmodel = std::make_shared<vec> (hyper(axes));
        if (nx!=0) fillin(bsmodel,input,controlx,controlz);
        else fillin1d(bsmodel,input,controlz);
    }
    else{
        bsmodel = read<data_t>(bsmodel_file,format);
        successCheck(bsmodel->getHyper()->isCompatible(hyper(axes)),__FILE__,__LINE__,"The B-splines model is not compatible with the input and B-splines parameters\n");
    }
  
    if (nx!=0){
        duplicate D(*bsmodel->getHyper(),mx,mz);
        bsplines3 B(*D.getRange(),*input->getHyper(),kx,kz);
        chainLOper * BD = new chainLOper(&B,&D);
        BD->forward(false,bsmodel,output);
        delete BD;
    }
    else{
        duplicate1d D(*bsmodel->getHyper(),mz);
        bsplines31d B(*D.getRange(),*input->getHyper(),kz);
        chainLOper * BD = new chainLOper(&B,&D);
        BD->forward(false,bsmodel,output);
        delete BD;
    }
// ----------------------------------------------------------------------------------------//

    if (output_file!="none") write<data_t>(output,output_file,format,datapath);

    return 0;
}



