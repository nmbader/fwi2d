#ifdef DOUBLE_PRECISION
    typedef double data_t;
#else
    typedef float data_t;
#endif

#include <unistd.h>

#include "param.hpp"
#include "operator.hpp"
#include "IO.hpp"
#include "seplib.h"

typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;

void printdoc(){
    std::string doc = "\nDescription:\n"
    "   Compute deconvolution filters between consecutive traces: filter_n=trace_n/trace_n-1 (division in frequency domain)\n"
    "   First filter is always a Dirac.\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'.\n"
    "\nParameters:\n"
    "   fmin - 0<=float<=1 - [0.0]:\n\t\tlowest frequency for spectral division.\n"
    "   fmax - 0<=float<=1 - [1.0]:\n\t\thighest frequency for spectral division.\n"
    "   eps - 0<=float<=1 - [0.1]:\n\t\tadd white noise expressed in percentage of maximum amplitude spectrum of the input.\n"
    "   smth_half_length - int>0 - [0]:\n\t\thalf length in nb of samples of the triangular smoothing window applied to the amplitude spectrum after division.\n"
    "   filter - string - ['none'] : optional filters (previously computed) to map back to the input space. If provided, the Jacobian of the operator will be applied to the filters wrt input.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries when format=1 is used.\n"
    "\nExample:\n"
    "   TR2TRDECON.x < infile.H fmin=0.1 fmax=0.5 eps=0.05 smth_half_length=3 > oufile.H\n"
    "   TR2TRDECON.x < infile.H filter=filters.H fmin=0.1 fmax=0.5 eps=0.05 smth_half_length=3 > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);

    std::string input_file="in", output_file="out", filter_file="none", datapath="none";
    data_t fmin=0, fmax=1, eps=0.1;
    int hl=0;
    bool format=0;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
	readParam<std::string>(argc, argv, "filter", filter_file);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<data_t>(argc, argv, "fmin", fmin);
    readParam<data_t>(argc, argv, "fmax", fmax);
    readParam<data_t>(argc, argv, "eps", eps);
    readParam<int>(argc, argv, "smth_half_length", hl);
    readParam<bool>(argc, argv, "format", format);

    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");
    
    std::shared_ptr<vec> input = read<data_t>(input_file, format);
    std::shared_ptr<vec> filter;
    if (filter_file !="none") filter = read<data_t>(filter_file, format);

    // compute the maxium of the amplitude spectrum of the input
    fxTransform fx(*input->getHyper());
    std::shared_ptr<cvec> fxvec = std::make_shared<cvec> (*fx.getRange());
    fxvec->zero();
    fx.forward(false,input,fxvec);
    std::shared_ptr<vec> spec=fxvec->modulus();
    data_t max_spec=spec->max();

    // perform the trace to trace deconvolution or its jacobian
    tr2trDecon decon(*input->getHyper(),fmin,fmax,eps*max_spec,hl);
    std::shared_ptr<vec> output = std::make_shared<vec>(*decon.getRange());
    output->zero();

    if (filter_file=="none") decon.forward(false,input,output);
    else decon.jacobianT(false,output,input,filter);

    output->setHyper(*input->getHyper());
    if (output_file!="none") write<data_t>(output, output_file, format, datapath);

    return 0;
}