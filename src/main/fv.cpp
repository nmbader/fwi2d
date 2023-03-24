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
    "   Compute f-v spectra using Fourier transform along the fast axis (often time) followed by phase shifts and stacking along the second fast axis.\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'.  '<' and '>' are valid for SEPlib format only.\n"
    "\nParameters:\n"
    "   type - string - ['amplitude']:\n\t\toptions: 'amplitude', 'power', 'phase', 'real', 'imag'.\n"
    "   vmin - float - [1000]:\n\t\tMinimum phase velocity to scan.\n"
    "   vmax - float - [5000]:\n\t\tMaximum phase velocity to scan (must be > vmin).\n"
    "   nv - int - [101]:\n\t\tNumber of velocity values to scan (must be > 1).\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries.\n"
    "\nExample:\n"
    "   FV.x < infile.H type=power > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);
    omp_set_num_threads(1);

    std::string input_file="in", output_file="out", type="amplitude", datapath="none";
    bool format=0;
    int nv=101;
    data_t vmin=1000, vmax=5000;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "type", type);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<data_t>(argc, argv, "vmin", vmin);
    readParam<data_t>(argc, argv, "vmax", vmax);
    readParam<int>(argc, argv, "nv", nv);
    readParam<bool>(argc, argv, "format", format);

    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");
    successCheck((vmax>vmin) && (nv>1),__FILE__,__LINE__,"Must have vmin>0, vmax>vmin, and nv>1\n");
    
    std::shared_ptr<vec> input = read<data_t>(input_file, format);

    // Fourier transform in time
    fxTransform fx(*input->getHyper());
    std::shared_ptr<cvec> fxvec = std::make_shared<cvec> (*fx.getRange());
    fxvec->zero();
    fx.forward(false,input,fxvec);

    // Prepare the FV vector 
    std::vector<ax> axes = fxvec->getHyper()->getAxes();
    data_t dv = (vmax-vmin)/(nv-1);
    ax F = axes[0];
    ax X = axes[1];
    axes[1].n=nv;
    axes[1].o=vmin;
    axes[1].d=dv;
    std::shared_ptr<cvec> fvvec = std::make_shared<cvec>(hyper(axes));
    fvvec->zero();

    // Populate the FV vector
    int ny = fvvec->getN123() / (F.n*nv);
    const std::complex<data_t> (* pin) [X.n][F.n] = (const std::complex<data_t> (*) [X.n][F.n]) fxvec->getCVals();
    std::complex<data_t> (* pout) [nv][F.n] = (std::complex<data_t> (*) [nv][F.n]) fvvec->getVals();

    for (int iy=0; iy<ny; iy++)
    {
        for (int iv=0; iv<nv; iv++)
        {
            data_t v = vmin + dv*iv;
            for (int ix=0; ix<X.n; ix++)
            {
                data_t x = X.o + ix*X.d;
                for (int fi=0; fi<F.n; fi++)
                {
                    data_t arg = 2*M_PI*(F.o+fi*F.d)*x/v;
                    pout[iy][iv][fi] += pin[iy][ix][fi] * std::complex<data_t>(cos(arg), sin(arg));
                }
            }
        }
    }

    std::shared_ptr<vec> output;
    if(type=="power") output = fvvec->modulus2();
    else if (type=="phase") output = fvvec->arg();
    else if (type=="real") output = fvvec->real();
    else if (type=="imag") output = fvvec->imag();
    else output = fvvec->modulus();
     
    if (output_file!="none") write<data_t>(output, output_file, format, datapath);

    return 0;
}