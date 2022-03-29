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

static data_t sinc(data_t x){
    if (x==0.0) return 1.0;
    else return sin(x)/x;
}

static data_t bessel0( data_t x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=0.  */
/*------------------------------------------------------------*/
{
   data_t ax,ans;
   data_t y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } else {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))));
   }
   return ans;
}

void printdoc(){
    std::string doc = "\nDescription:\n"
    "   Resample along the fast axis (often time), in time or frequency domain.\n"
    "   The origin remains fixed and the maximum time of the output is <= then that of the input.\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'. '<' and '>' are valid for SEPlib format only.\n"
    "\nParameters:\n"
    "   domain - string - ['time']:\n\t\t'time' or 'frequency'.\n"
    "   type - string - ['sinc']:\n\t\t'sinc' or 'kaiser' or 'linear'. 'linear' and 'kaiser' are valid for time domain only. 'kaiser' uses a sinc filter with a Kaiser window.\n"
    "   factor - int - [2]:\n\t\tfactor by which to interpolate.\n"
    "   sinc_half_length - int - [21]:\n\t\thalf length in nb of samples of the filter used for sinc and kaiser interpolation.\n"
    "   si - non-negative float - [0]:\n\t\tdesired output sampling. If > 0 it overwrites 'factor' and is valid for time domain only.\n"
    "   alpha - 0<=float<=1 or >0 - [0.5]:\n\t\tused in the cosine window for the sinc filter (1 means no taper). It is also used for kaiser window; recommended value 10.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries.\n"
    "\nExample:\n"
    "   RESAMPLE.x < infile.H domain=time type=linear si=0.1 > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);
    omp_set_num_threads(1);

    std::string input_file="in", output_file="out", domain="time", type="sinc", datapath="none";
    bool format=0;
    int a=2, hl=21;
    data_t si=0, alpha=0.5;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "domain", domain);
    readParam<std::string>(argc, argv, "type", type);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<bool>(argc, argv, "format", format);
    readParam<int>(argc, argv, "factor", a);
    readParam<int>(argc, argv, "sinc_half_length", hl);
    readParam<data_t>(argc, argv, "si", si);
    readParam<data_t>(argc, argv, "alpha", alpha);

    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");
    if ((domain=="frequency") && ((si>0) || (type=="linear"))) {
        fprintf(stderr,"\nParameter 'si>0' and 'type=linear' are to be used for time domain resampling only.\n");
        return 0;
    }

    std::shared_ptr<vec> input = read<data_t>(input_file, format);

    int n123 = input->getN123();
    std::vector<ax > axes = input->getHyper()->getAxes();

    ax T = axes[0];
    int nx = round(n123/T.n);
    if (si>0) {axes[0].n= floor(((axes[0].n-1)*axes[0].d)/si)+1; axes[0].d=si;}
    else{
        axes[0].n = (axes[0].n-1)*a+1;
        axes[0].d = axes[0].d/a;
    }

    ax T2 = axes[0];
    std::shared_ptr<vec> output = std::make_shared<vec>(hyper(axes));
    output->zero();

    if (domain=="time"){
        if (type == "linear"){
            linear_resampler R(*input->getHyper(), *output->getHyper());
            R.forward(false, input, output);
        }
        else if(type=="sinc"){
            sinc_resampler R(*input->getHyper(), *output->getHyper(), hl, alpha);
            R.forward(false, input, output);
        }
        else{
            data_t t;
            int itm0;
            const data_t * pin = input->getCVals();
            data_t * pout = output->getVals();
            for (int i=0; i<nx; i++){
                for (int itd=0; itd<T2.n; itd++){
                    t = itd * T2.d;
                    itm0 = floor(t / T.d);
                    for (int itm = std::max(0, itm0 -hl+1) ; itm <= std::min(T.n-1, itm0 + hl); itm++){
                        pout[i*T2.n+itd] += pin[i*T.n+itm]*sinc(M_PI * (t-itm*T.d) / T.d)*bessel0(alpha*std::sqrt(1-(itm-t/T.d)*(itm-t/T.d)/(hl*hl)))/bessel0(alpha);
                    }
                }
            }
        }
    }

    else{

        fxTransform fx(*input->getHyper());
        ax F = fx.getRange()->getAxis(1);

// 1D Fourier transform over the time axis
        std::shared_ptr<cvec> fxvec = std::make_shared<cvec >(*fx.getRange());
        fx.forward(false,input,fxvec);

// Pad with zeros
        ax F2 = F; F2.n = axes[0].n/2 + 1;
        axes[0] = F2;
        std::shared_ptr<cvec> fxvec2 = std::make_shared<cvec>(hyper(axes));
        std::complex<data_t> c0 (0,0);
        std::complex<data_t> * pf = fxvec->getVals();
        std::complex<data_t> * pf2 = fxvec2->getVals();
        for (int i=0; i<nx; i++){
            for (int j=0; j<F.n; j++){
                pf2[i*F2.n+j]=(data_t)a*pf[i*F.n+j];
            }
            for (int j=F.n; j<F2.n; j++){
                pf2[i*F2.n+j]=c0;
            }
        }

// Inverse Fourier transform
        fxTransform fx2(*output->getHyper());
        fx2.inverse(false,output,fxvec2);
    }
          
    if (output_file!="none") write<data_t>(output, output_file, format, datapath);

    return 0;
}