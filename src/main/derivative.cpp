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
    "   Apply first derivative along the fast axis (often time), in time or frequency domain.\n"
    "   For time domain, a centered difference is applied in the interior and Euler forward/backward at the boundaries.\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'. '<' and '>' are valid for SEPlib format only.\n"
    "\nParameters:\n"
    "   domain - string - ['time']:\n\t\t'time' or 'frequency'.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries when format=1 is used.\n"
    "\nExample:\n"
    "   DERIVATIVE.x < infile.H domain=frequency > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);

    std::string input_file="in", output_file="out", domain="time", datapath="none";
    bool format=0;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "domain", domain);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<bool>(argc, argv, "format", format);

    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");
    
    std::shared_ptr<vec> input = read<data_t>(input_file,format);

    if (domain=="time")
    {
        std::shared_ptr<vec> temp = input->clone();
        ax T = input->getHyper()->getAxis(1);
        int nx = input->getN123() / T.n; 
        Dt(false, false, temp->getCVals(), input->getVals(), nx, T.n, T.d, 0, nx);
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
                pfx[i*F.n+j] *= std::complex<data_t>(0,omega);
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