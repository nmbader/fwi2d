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
    "   Apply a provided filter along the fast axis (often time).\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'.  '<' and '>' are valid for SEPlib format only.\n"
    "\nParameters:\n"
    "   filter - string - ['none']:\n\t\t1D filter file to be applied to every input trace. If 'none', then output=input\n"
    "   phase - string - ['default']:\n\t\toptions: 'default', 'zero', 'minimum'.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries.\n"
    "\nExample:\n"
    "   FX_FILTER.x < infile.H filter=filfile.H domain=time > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);
    omp_set_num_threads(1);

    std::string input_file="in", output_file="out", filter_file="none", phase="default", domain="time", datapath="none";
    bool format=0;
    data_t eps=1e-07;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
	readParam<std::string>(argc, argv, "filter", filter_file);
    readParam<std::string>(argc, argv, "phase", phase);
    readParam<data_t>(argc, argv, "epsilon", eps); // used in Kolgomoroff factorization for minimum phase
    readParam<std::string>(argc, argv, "domain", domain);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<bool>(argc, argv, "format", format);

    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");
    
    std::shared_ptr<vec> input = read<data_t>(input_file, format);
    std::shared_ptr<vec> output;
    std::shared_ptr<vec> filter;
    if (filter_file != "none") 
    {
        filter = read<data_t>(filter_file, format);
        if (phase=="zero") filter = zero_phase(filter);
        else if (phase=="minimum") filter = minimum_phase(filter,eps);
        ax Tf = filter->getHyper()->getAxis(1);
        ax T = input->getHyper()->getAxis(1);
        if (domain == "time"){
            Tf.d=T.d;
            filter->setHyper(hyper(Tf)); 
            conv1dnd op(*input->getHyper(), filter, phase!="minimum");
            output = std::make_shared<vec>(*op.getRange());
            output->zero();
            op.forward(false,input,output);
        }
        else{
            // pad or cut the filter if necessary
            int nx = input->getN123()/T.n;
            int ntf = filter->getN123();
            std::shared_ptr<vec> fil = std::make_shared<vec>(hyper(T));
            fil->zero();
            memcpy(fil->getVals(), filter->getVals(), sizeof(data_t)*std::min(T.n,ntf));

            // forward Fourier transform
            fxTransform fx(*input->getHyper());
            fxTransform fxf(*fil->getHyper());
            std::shared_ptr<cvec> fxvec = std::make_shared<cvec> (*fx.getRange());
            std::shared_ptr<cvec> fxfvec = std::make_shared<cvec> (*fxf.getRange());
            fxvec->zero();
            fxfvec->zero();
            fx.forward(false,input,fxvec);
            fxf.forward(false,fil,fxfvec);
            int nf = fxf.getRange()->getAxis(1).n;

            // multiply by the filter amplitude spectrum (zero phase convolution)
            for (int ix=0; ix<nx; ix++){
                for (int f=0; f<nf; f++) fxvec->getVals()[ix*nf+f] *= std::abs(fxfvec->getVals()[f]);
            }

            // inverse Fourier transform
            fx.inverse(false,input,fxvec);
            output=input;
        }
    }
    else output = input;
     
    if (output_file!="none") write<data_t>(output, output_file, format, datapath);

    return 0;
}