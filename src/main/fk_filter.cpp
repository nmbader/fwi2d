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
    "   Apply F-K filter to a dataset (tapered mute in the F-K domain).\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'.  '<' and '>' are valid for SEPlib format only.\n"
    "   Input first axis is interpreted as time, second axis as space.\n"
    "\nParameters:\n"
    "   kmin - -0.5<=float<=0.5 - [-0.5]:\n\t\tlowest wavenumber to fully preserve (-0.5 means no boxcar filtering).\n"
    "   kmax - kmin<=float<=0.5 - [0.5]:\n\t\thighest wavenumber to fully preserve (0.5 means no boxcar filtering).\n"
    "   vmin - non-negative float - [0]:\n\t\thighest absolute phase velocity (m/s) to fully remove (0 means no cone filtering).\n"
    "   vmax - vmin<=float - [0]:\n\t\tlowest absolute phase velocity (m/s) to fully preserve (0 means no cone filtering).\n"
    "   taper - non-negative float - [0.1]:\n\t\tlength of the cosine taper, applied in the wavenumber direction.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries.\n"
    "\nExample:\n"
    "   FK_FILTER.x < infile.H kmin=-0.3 kmax=0.3 vmin=1000 vmax=10000 taper=0.05 > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);
    omp_set_num_threads(1);

    std::string input_file="in", output_file="out", datapath="none";
    bool format=0;
    data_t kLow=-0.5, kHigh=0.5, taper=0.1, vmin=0, vmax=0;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<bool>(argc, argv, "format", format);
    readParam<data_t>(argc, argv, "kmin", kLow);
    readParam<data_t>(argc, argv, "kmax", kHigh);
    readParam<data_t>(argc, argv, "vmin", vmin);
    readParam<data_t>(argc, argv, "vmax", vmax);
    readParam<data_t>(argc, argv, "taper", taper);

    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");
    successCheck((abs(kLow)<=0.5) && (abs(kHigh)<=0.5) && (kLow<=kHigh),__FILE__,__LINE__,"The cutoff wavenumbers must obey -0.5 <= kLow <= kHigh <= 0.5\n");

    vmax = abs(vmax);
    vmin = abs(vmin);
    successCheck(vmin<=vmax,__FILE__,__LINE__,"The cutoff velocities must obey 0 <= vmin <= vmax\n");

    std::shared_ptr<vec> input = read<data_t>(input_file, format);

    int n123 = input->getN123();
    std::vector<ax > axes = input->getHyper()->getAxes();

    ax T = axes[0];
    ax X = axes[1];
    int ny = round(n123/(T.n*X.n));
    ax Y(ny,0,1);

    fkTransform fk(*input->getHyper());
    ax F = fk.getRange()->getAxis(1);
    X = fk.getRange()->getAxis(2);

    std::shared_ptr<cvec> fkvec = std::make_shared<cvec >(*fk.getRange());
    fk.forward(false, input, fkvec);

    std::complex<data_t> (*pvec) [X.n][F.n] = (std::complex<data_t> (*) [X.n][F.n]) fkvec->getVals();

    for (int i=0; i<ny; i++){

        int nx1 = (0.5+kLow-taper)*X.n;
        nx1 = std::max(nx1, 0);
        int nx2 = (0.5+kLow)*X.n;
        int nx3 = std::ceil((0.5+kHigh)*X.n);
        int nx4 = (0.5+kHigh+taper)*X.n;
        nx4 = std::min(nx4, X.n);

        if ((kLow != -0.5) || (kHigh != 0.5)){

            for (int ix=0; ix<nx1; ix++){
                for (int iz=0; iz<F.n; iz++){
                    pvec[i][ix][iz] = 0;
                }
            }
            for (int ix=nx1; ix<nx2; ix++){
                for (int iz=0; iz<F.n; iz++){
                    pvec[i][ix][iz] *= (1 + std::cos(M_PI*(kLow + 0.5 - (1.0*ix)/X.n)/taper))/2.0;
                }
            }
            for (int ix=nx3; ix<nx4; ix++){
                for (int iz=0; iz<F.n; iz++){
                    pvec[i][ix][iz] *= (1 + std::cos(M_PI*((1.0*ix)/X.n - 0.5 - kHigh)/taper))/2.0;
                }
            }
            for (int ix=nx4; ix<X.n; ix++){
                for (int iz=0; iz<F.n; iz++){
                    pvec[i][ix][iz] = 0;
                }
            }
        }

        if (vmin != 0){
            data_t f, k, k1, k2;
            int imin, imax;
            // negative wavenumbers
            for (int ix=0; ix<X.n/2; ix++){
                k = X.o+ix*X.d;
                imin = floor((-k*vmin)/(2*M_PI*F.d));
                imin = std::min(F.n-1, imin);
                imax = floor((-k*vmax)/(2*M_PI*F.d));
                imax = std::min(F.n-1, imax);
                for (int iz=0; iz<=imin; iz++) pvec[i][ix][iz] = 0;
                for (int iz=imin+1; iz<=imax; iz++){
                    f = iz*F.d;
                    k1 = -2*M_PI*f / vmax;
                    k2 = -2*M_PI*f / vmin;
                    pvec[i][ix][iz] *= 0.5*(1+cos(M_PI*(k1-k)/(k1-k2)));
                }
            }
            // positive wavenumbers
            for (int ix=X.n/2; ix<X.n; ix++){
                k = X.o+ix*X.d;
                imin = floor((k*vmin)/(2*M_PI*F.d));
                imin = std::min(F.n-1, imin);
                imax = floor((k*vmax)/(2*M_PI*F.d));
                imax = std::min(F.n-1, imax);
                for (int iz=0; iz<=imin; iz++) pvec[i][ix][iz] = 0;
                for (int iz=imin+1; iz<=imax; iz++){
                    f = iz*F.d;
                    k1 = 2*M_PI*f / vmax;
                    k2 = 2*M_PI*f / vmin;
                    pvec[i][ix][iz] *= 0.5*(1+cos(M_PI*(k-k1)/(k2-k1)));
                }
            }
        }
    }

    fk.inverse(false, input, fkvec);
     
    if (output_file!="none") write<data_t>(input, output_file, format, datapath);

    return 0;
}