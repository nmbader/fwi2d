#include <unistd.h>

#include "seplib.h"
#include "param.hpp"
#include "operator.hpp"
#include "IO.hpp"

typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;

void printdoc(){
    std::string doc = "\nDescription:\n"
    "   Integrate along the fast axis (often time), in time or frequency domain.\n"
    "   For time domain, a trapezoidal quadrature is used for the integration.\n"
    "   For frequency domain, a damping is added to avoid dividing by zero.\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'.  '<' and '>' are valid for SEPlib format only.\n"
    "\nParameters:\n"
    "   domain - string - ['time']:\n\t\t'time' or 'frequency'.\n"
    "   eps - positive float - [1e-06]:\n\t\tdamping value.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries.\n"
    "\nExample:\n"
    "   INTEGRAL.x < infile.H domain=frequency eps=1e-04 > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);
    omp_set_num_threads(1);

    std::string input_file="in", output_file="out", domain="time", datapath="none";
    bool format=0;
    data_t eps=1e-6;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "domain", domain);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<data_t>(argc, argv, "eps", eps);
    readParam<bool>(argc, argv, "format", format);

    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");
    
    std::shared_ptr<vec> input = read<data_t>(input_file, format);

    if (domain=="time")
    {
        integral op(*input->getHyper());
        op.forward(false,input,input);
    }
    else if (domain=="frequency")
    {
        fxTransform fx(*input->getHyper());
        std::shared_ptr<cvec> fxvec = std::make_shared<cvec> (*fx.getRange());
        fxvec->zero();
        fx.forward(false,input,fxvec);

        std::complex<data_t> * pfx = fxvec->getVals();
        data_t omega=0;
        ax F = fxvec->getHyper()->getAxis(1);
        int nx = input->getN123() / F.n;
        for (int i=0; i<nx; i++){
            for (int j=0; j<F.n; j++){
                omega = 2*M_PI*j*F.d;
                pfx[i*F.n+j] /= std::complex<data_t>(0,omega+eps);
            }
        }

        fx.inverse(false,input,fxvec);
    }
    else{
        successCheck(false, __FILE__,__LINE__,"domain must be time or frequency\n");
    }
          
    if (output_file!="none") write<data_t>(input, output_file, format, datapath);

    return 0;
}