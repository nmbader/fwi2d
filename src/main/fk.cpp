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
    "   Compute 2D spectra using Fourier transform along the two fastest axes (often time and first spatial dimension).\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'.  '<' and '>' are valid for SEPlib format only.\n"
    "\nParameters:\n"
    "   type - string - ['amplitude']:\n\t\toptions: 'amplitude', 'power', 'phase', 'real', 'imag'.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries.\n"
    "\nExample:\n"
    "   FK.x < infile.H type=power > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);
    omp_set_num_threads(1);

    std::string input_file="in", output_file="out", type="amplitude", datapath="none";
    bool format=0;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "type", type);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<bool>(argc, argv, "format", format);

    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");
    
    std::shared_ptr<vec> input = read<data_t>(input_file, format);

    fkTransform fk(*input->getHyper());
    std::shared_ptr<cvec> fkvec = std::make_shared<cvec> (*fk.getRange());
    fkvec->zero();
    fk.forward(false,input,fkvec);

    std::shared_ptr<vec> output;
    if(type=="power") output = fkvec->modulus2();
    else if (type=="phase") output = fkvec->arg();
    else if (type=="real") output = fkvec->real();
    else if (type=="imag") output = fkvec->imag();
    else output = fkvec->modulus();

    output->setHyper(*fk.getRange());
     
    if (output_file!="none") write<data_t>(output, output_file, format, datapath);

    return 0;
}