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
    "   Hilbert transform (90 deg phase shift) along the fast axis (often time).\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'.  '<' and '>' are valid for SEPlib format only.\n"
    "\nParameters:\n"
    "   type - string - ['hilbert']:\n\t\toptions: 'hilbert', 'envelop', 'iphase'.\n"
    "   wrapped - bool - ['1']:\n\t\twrap or unwrap the instantaneous phase.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries when format=1 is used.\n"
    "\nExample:\n"
    "   HILBERT.x < infile.H > oufile.H.\n"
    "   HILBERT.x < infile.H type=iphase wrapped=0 > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);

    std::string input_file="in", output_file="out", type="hilbert", datapath="none";
    bool wrapped=true, format=0;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "type", type);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<bool>(argc, argv, "wrapped", wrapped);
    readParam<bool>(argc, argv, "format", format);

    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");
    
    std::shared_ptr<vec> input = read<data_t>(input_file, format);

// transfer input to a complex vector
    std::shared_ptr<cvec> cinput = std::make_shared<cvec> (*input->getHyper());
    data_t * pin = input->getVals();
    std::complex<data_t> * pcin = cinput->getVals();
    for (int i=0; i<input->getN123(); i++) pcin[i] = pin[i];

// forward transform
    fxTransform fx(*cinput->getHyper());
    std::shared_ptr<cvec> fxvec = std::make_shared<cvec >(*fx.getCRange());
    fx.cforward(false, cinput, fxvec);
    
// zero the negative frequencies
    int nf = fxvec->getHyper()->getAxis(1).n;
    int nx = fxvec->getN123() / nf;
    std::complex<data_t> * p = fxvec->getVals();
    for (int i=0; i<nx; i++){
        for (int j=0; j<(nf-1)/2; j++) p[i*nf+j] *= 0.0;
    }

// inverse transform
    fx.cinverse(false,cinput,fxvec);
    std::shared_ptr<vec> hil = cinput->imag();
    hil->scale(2);
    data_t * phil = hil->getVals();

    if (type=="envelop")
    {
        for (int i=0; i<input->getN123(); i++) pin[i] = sqrt(pin[i]*pin[i] + phil[i]*phil[i]);
    }
    else if (type == "iphase") // unwrapped instantaneous phase
    {
        if (wrapped)
        {
            for (int i=0; i<input->getN123(); i++) 
            {
                pin[i] = std::atan2(phil[i],pin[i]);
            }
        }
        else
        {
            std::complex<data_t> val1(pin[0],phil[0]);
            std::complex<data_t> val2(pin[0],phil[0]);
            pin[0] = std::arg(val1);
            //pin[0] = 0;
            for (int i=1; i<input->getN123(); i++) 
            {
                val1=val2;
                val2.real(pin[i]);
                val2.imag(phil[i]);
                pin[i] = pin[i-1] + std::arg(val2*std::conj(val1));
            }
        }
    }   
    else // Hilbert transform (90 deg phase shift)
    {
        for (int i=0; i<input->getN123(); i++) pin[i] = phil[i];
    }      

    if (output_file!="none") write<data_t>(input, output_file, format, datapath);

    return 0;
}