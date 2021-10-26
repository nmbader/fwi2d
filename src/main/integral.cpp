#ifdef DOUBLE_PRECISION
    typedef double data_t;
#else
    typedef float data_t;
#endif

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
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'.\n"
    "\nParameters:\n"
    "   domain - string - ['time']:\n\t\t'time' or 'frequency'.\n"
    "   eps - positive float - [1e-06]:\n\t\tdamping value.\n"
    "\nExample:\n"
    "   INTEGRAL.x < infile.H domain=frequency eps=1e-04 > oufile.H.\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);

    std::string input_file="in", output_file="out", domain="time";
    data_t eps=1e-6;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "domain", domain);
    readParam<data_t>(argc, argv, "eps", eps);

    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");
    
    std::shared_ptr<vec> input = sepRead<data_t>(input_file);

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
          
    if (output_file!="none") sepWrite<data_t>(input, output_file);

    return 0;
}