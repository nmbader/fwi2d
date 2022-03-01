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

inline data_t sinc(data_t x){
    if (x==0.0) return 1.0;
    else return sin(x)/x;
}

std::shared_ptr<vec> rickerWavelet(data_t wc , int N){
    std::shared_ptr<vec> v = std::make_shared<vec> (hyper(2*N+1));
    data_t * p = v->getVals();
    data_t num = 0;
	data_t sig = sqrt(2)/(M_PI*wc);
	for (int iz=N; iz<2*N+1; iz++){ 
		num = (1-(iz-N)*(iz-N)/(sig*sig))*exp(-0.5*(iz-N)*(iz-N)/(sig*sig));
		p[iz] = num;
		p[2*N-iz] = num;
	}
	return v;
}
std::shared_ptr<vec> gaussianWavelet(data_t sig, int N){
    std::shared_ptr<vec> v = std::make_shared<vec> (hyper(2*N+1));
    data_t * p = v->getVals();
    data_t num = 0;
    for (int iz=N; iz<2*N+1; iz++){ 
		num = exp(-0.5*(iz-N)*(iz-N)/(sig*sig));
		p[iz] = num;
		p[2*N-iz] = num;
	}
	return v;
}
std::shared_ptr<vec> sincWavelet(data_t wc, int N, data_t alpha){
    std::shared_ptr<vec> v = std::make_shared<vec> (hyper(2*N+1));
    data_t * p = v->getVals();
    data_t num;
	for (int iz=0; iz<2*N+1; iz++){
		num = wc*sinc(wc*M_PI*(iz-N))*(alpha +(1-alpha)*cos(M_PI*(iz-N)/N));
		p[iz] = num;
	}
	return v;
}
void iirButterworth(data_t * p, int nt, data_t wc, bool lp, int ho){
    if (lp == false) wc = 1-wc;

    // Initialization
    int M = 2*ho;
    data_t Omega_c = tan(M_PI*wc/2);
    data_t a0, a1, a2, b1, b2, c;
    data_t xi, xi1, xi2;

    // Apply the same filter ho times
    for (int k=0; k<ho; k++){

        // compute filter coefficients for Low-Pass
        c=1+2*cos(M_PI*(2*k+1)/(2*M))*Omega_c + Omega_c*Omega_c;
        a0=Omega_c*Omega_c/c;
        a1 = 2 * a0;
        a2 = a0;
        b1=2*(Omega_c*Omega_c-1)/c;
        b2=(1-2*cos(M_PI*(2*k+1)/(2*M))*Omega_c +Omega_c*Omega_c)/c;

        // deduce filter coefficients for High-Pass
        if (lp == false){
            a0 = a0;
            a1 = - a1;
            a2 = a2;
            b1 = - b1;
            b2 = b2;
        }

        // apply filter recursively
        xi = p[1];
        xi1 = p[0];
        xi2 = 0;
        p[0] = a0*xi1;
        p[1] = a0*xi + a1*xi1 - b1*p[0];

        for (int iz=2; iz<nt; iz++){
            xi2 = xi1;
            xi1 = xi;			
            xi = p[iz];
            p[iz] = a0*xi + a1*xi1 + a2*xi2 
                    - b1*p[iz-1] -b2*p[iz-2];
        }
    }
}
void iirBpButterworth(data_t * p, int nt, data_t wcl, data_t wch, int ho){
    
    // Initialization
    int M = 2*ho;
    data_t Omega_ch = tan(M_PI*wch/2);
    data_t Omega_cl = tan(M_PI*(1-wcl)/2);
    data_t a0h, a1h, a2h, b1h, b2h, ch;
    data_t a0l, a1l, a2l, b1l, b2l, cl;
    data_t a0, a1, a2, a3, a4, b1, b2, b3, b4;
    data_t xi, xi1, xi2, xi3, xi4;

    // Apply the same filter ho times
    for (int k=0; k<ho; k++){

        // compute filter coefficients for Low-Pass
        ch=1+2*cos(M_PI*(2*k+1)/(2*M))*Omega_ch + Omega_ch*Omega_ch;
        a0h=Omega_ch*Omega_ch/ch;
        a1h = 2 * a0h;
        a2h = a0h;
        b1h=2*(Omega_ch*Omega_ch-1)/ch;
        b2h=(1-2*cos(M_PI*(2*k+1)/(2*M))*Omega_ch +Omega_ch*Omega_ch)/ch;

        // compute filter coefficients for High-Pass
        cl=1+2*cos(M_PI*(2*k+1)/(2*M))*Omega_cl + Omega_cl*Omega_cl;
        a0l=Omega_cl*Omega_cl/cl;
        a1l = - 2 * a0l;
        a2l = a0l;
        b1l=-2*(Omega_cl*Omega_cl-1)/cl;
        b2l=(1-2*cos(M_PI*(2*k+1)/(2*M))*Omega_cl +Omega_cl*Omega_cl)/cl;

        // deduce filter coefficients for Band-Pass
        a0 = a0h*a0l;
        a1 = a0h*a1l + a0l*a1h;
        a2 = a0h*a2l + a1h*a1l + a2h*a0l;
        a3 = a1h*a2l + a2h*a1l;
        a4 = a2h*a2l;
        b1 = b1l + b1h;
        b2 = b2l + b1h*b1l + b2h;
        b3 = b1h*b2l + b2h*b1l;
        b4 = b2h*b2l;
        
        // apply filter recursively
        xi = p[3];
        xi1 = p[2];
        xi2 = p[1];
        xi3 = p[0];
        xi4 = 0;
        p[0] = a0*xi3;
        p[1] = a0*xi2 + a1*xi3 - b1*p[0];
        p[2] = a0*xi1 + a1*xi2 + a2*xi3 - b1*p[1] - b2*p[0];
        p[3] = a0*xi + a1*xi1 + a2*xi2 + a3*xi3 - b1*p[2] - b2*p[1] - b3*p[0];

        for (int iz=4; iz<nt; iz++){
            xi4 = xi3;
            xi3 = xi2;
            xi2 = xi1;
            xi1 = xi;			
            xi = p[iz];
            p[iz] = a0*xi + a1*xi1 + a2*xi2 + a3*xi3 + a4*xi4
                    - b1*p[iz-1] -b2*p[iz-2] -b3*p[iz-3] -b4*p[iz-4];
        }
    }
}
std::shared_ptr<vec> zero_phase(const std::shared_ptr<vec> dat){

    ax T = dat->getHyper()->getAxis(1);
    ax X(1,0,1);
    
    // Make a new vector with odd length
    std::shared_ptr<vec> vec0;
	if (2*(T.n/2)==T.n) {
		T.n+=1;
        vec0 = std::make_shared<vec>(hyper(T));

        for (int iz=0; iz<T.n-1; iz++){
        vec0->getVals()[iz] = dat->getVals()[iz];
        }
        vec0->getVals()[T.n-1] = 0;
    }
	
    else {
        vec0 = dat;
        vec0->setHyper((hyper(T,X)));
    }

    // Fourier transform the vector
    fxTransform fx(hyper(T,X));
    ax F;
    F.d=fx.getRange()->getAxis(1).d;
    F.o=0;
    F.n=T.n/2 + 1;
    std::shared_ptr<cvec> fxvec = std::make_shared<cvec >(hyper(F,X));
    fx.forward(false,vec0,fxvec);

    int N = F.n;
    data_t dw = 2*M_PI/T.n;
    data_t r;
    std::complex<data_t> z (0,0);

    // Modify the FT by setting a linear phase while keeping the same amplitude spectrum
	for (int i=0; i<N; i++){

        r = abs(fxvec->getVals()[i]);

        // add linear phase shift (* exp (-j.omega.dt.(Nt-1)/2)))
        z.real(r*cos(-i*dw*(T.n - 1)/2));
        z.imag(r*sin(-i*dw*(T.n - 1)/2));

        fxvec->getVals()[i] = z;
	}

    // inverse FT and modify the original vector
    fx.inverse(false,vec0,fxvec);
    T.o = (int)(-T.n/2) * T.d; 
    vec0->setHyper(hyper(T));
    
    return vec0;
}
std::shared_ptr<vec> minimum_phase(const std::shared_ptr<vec> dat, const data_t eps){
    
    ax T = dat->getHyper()->getAxis(1);
    ax X(1,0,1);
    
    // Make a new vector with odd length
    std::shared_ptr<vec> vec0;
	if (2*(T.n/2)==T.n) {
		T.n+=1;
        vec0 = std::make_shared<vec>(hyper(T));

        for (int iz=0; iz<T.n-1; iz++){
        vec0->getVals()[iz] = dat->getVals()[iz];
        }
        vec0->getVals()[T.n-1] = 0;
    }
	
    else {
        vec0 = dat;
        vec0->setHyper((hyper(T,X)));
    }

    // Fourier transform the vector
    fxTransform fx(hyper(T,X));
    ax F;
    F.d=fx.getRange()->getAxis(1).d;
    F.o=0;
    F.n=T.n/2 + 1;
    std::shared_ptr<cvec> fxvec = std::make_shared<cvec >(hyper(F,X));
    fx.forward(false,vec0,fxvec);

    int N = F.n;
    data_t r;
    std::complex<data_t> z (0,0);

    // Modify the FT by making it real and equal to the logarithm of the power spectrum
	for (int i=0; i<N; i++){

        r = sqrt(norm(fxvec->getVals()[i]));

        // set FT = log(r+eps)
        z.real(log(r+eps));
        z.imag(0);

        fxvec->getVals()[i] = z;
	}

    // inverse FT and keep the causal part
    fx.inverse(false,vec0,fxvec);
    for (int iz=1; iz<T.n/2+1; iz++){
        vec0->getVals()[iz] *= 2;
    }
    for (int iz=T.n/2+1; iz<T.n; iz++){
        vec0->getVals()[iz] = 0;
    }

    // FT again and take the exponent
    fx.forward(false,vec0,fxvec);
    for (int iz=0; iz<N; iz++){
        fxvec->getVals()[iz] = exp(fxvec->getVals()[iz]);
    }

    // inverse FT 
    fx.inverse(false,vec0,fxvec);
    T.o = 0; 
    vec0->setHyper(hyper(T));
    
    return vec0;
}

void printdoc(){
    std::string doc = "\nDescription:\n"
    "   Generate a 1D wavelet.\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'.  '<' and '>' are valid for SEPlib format only.\n"
    "   Input is required only if the wavelet is built from a provided amplitude spectrum.\n"
    "\nParameters:\n"
    "   type - string - ['ricker']:\n\t\toptions: 'ricker', 'gaussian', 'sinc', 'butterworth', 'exponential', 'bandpass','spectrum','power'.\n"
    "   wc - 0<=float<=1 - [0.1]:\n\t\tcentral frequency for Ricker wavelet, given in fraction of the Nyquist.\n"
    "   sigma - float - [1]:\n\t\tstandard deviation (or damping) in axis units for Gaussian wavelet (or exponential wavelet)\n"
    "   low_cutoff - 0<=float<=1 - [0]:\n\t\tlow cutoff frequency for sinc and Butterworth wavelets, given in fraction of the Nyquist.\n"
    "   high_cutoff - low_cutoff<=float<=1 - [0]:\n\t\thigh cutoff frequency for sinc and Butterworth wavelets, given in fraction of the Nyquist.\n"
    "   half_order - int - [2]:\n\t\thalf of the order used in the default IIR Butterworth.\n"
    "   alpha - 0<=float<=1 - [0.5]:\n\t\tused in the cosine window for the sinc wavelet (1 means no taper).\n"
    "   w1,w2,w3,w4 - 0<=w1<w2<w3<w4<=1 - [0,0.1,0.6,0.9]:\n\t\tfrequencies defining band-pass wavelet, given in fraction of the Nyquist. A cosine taper is applied to the spectrum between w1 and w2, w3 and w4.\n"
    "   pow - float - [1]:\n\t\tpower used for the power function t^pow t>0.\n"
    "   phase - string - ['default']:\n\t\toptions: 'default', 'zero', 'minimum'.\n"
    "   nt - int - [100]:\n\t\tnumber of samples in the wavelet. If it is even, it will be augmented by 1.\n"
    "   dt - positive float - [0.004]:\n\t\tsampling interval.\n"
    "   shift - int - [0]:\n\t\tshift in nb of samples, applied to the origin of the output. Inactive if 'phase' is not provided.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries when format=1 is used.\n"
    "\nExample:\n"
    "   GENERATE_WAVELET.x < infile.H type=spectrum nt=100 dt=0.004 shift=11 phase=minimum > oufile.H\n"
    "   GENERATE_WAVELET.x type=bandpass nt=100 dt=0.004 shift=11 w1=0 w2=0.2 w3=0.6 w4=0.8 > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);

    std::string input_file="in", output_file="out", type="ricker", phase="default", datapath="none";
    data_t wc=0.1, dt=0.004, lwc=0, hwc=1, alpha=0.5, sigma=1, eps=1e-07, w1=0, w2=0.1, w3=0.6, w4=0.9, pow=1;
    int nt=100, shift=0, ho=2;
    bool format=0;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "type", type);
    readParam<std::string>(argc, argv, "phase", phase);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<data_t>(argc, argv, "wc", wc); // central freq for ricker
    readParam<data_t>(argc, argv, "sigma", sigma); // std dev in nb of samples for Gaussian (or exponential)
    readParam<data_t>(argc, argv, "low_cutoff", lwc); // cutoff for sinc and butterworth
    readParam<data_t>(argc, argv, "high_cutoff", hwc);
    readParam<data_t>(argc, argv, "w1", w1); // frequencies for bandpass
    readParam<data_t>(argc, argv, "w2", w2);
    readParam<data_t>(argc, argv, "w3", w3);
    readParam<data_t>(argc, argv, "w4", w4);
    readParam<data_t>(argc, argv, "pow", pow);
    readParam<data_t>(argc, argv, "alpha", alpha); // used for cosine window for sinc
    readParam<data_t>(argc, argv, "epsilon", eps); // used in Kolgomoroff factorization for minimum phase
    readParam<data_t>(argc, argv, "dt", dt);
    readParam<int>(argc, argv, "half_order", ho); // half order used for butterworth
    readParam<int>(argc, argv, "nt", nt);
    readParam<int>(argc, argv, "shift", shift);
    readParam<bool>(argc, argv, "format", format);

    if (type=="spectrum") {
        successCheck(input_file!="none",__FILE__,__LINE__,"Input file is expected for type=spectrum\n");
    }

    std::shared_ptr<vec> dat;

    int half_length = floor(nt/2);
    int length = 2*half_length+1;

    if (type=="ricker")
	    dat = rickerWavelet(wc, half_length);
    else if (type=="gaussian")
        dat=gaussianWavelet(sigma/dt, half_length);
    else if (type=="sinc"){
        std::shared_ptr<vec> lp = sincWavelet(hwc, half_length, alpha);
        std::shared_ptr<vec> hp = sincWavelet(lwc, half_length, alpha);
        hp->scale(-1);
        hp->getVals()[half_length] += 1;
        if (lwc==0)
            dat=lp->clone();
        else if (hwc==1)
            dat=hp->clone();
        else{
            dat=lp->clone();
            dat->zero();
            for (int i=half_length; i<length; i++){
                for (int j=0; j<=i; j++){
                    dat->getVals()[i-half_length] += lp->getVals()[i-j]*hp->getVals()[j];
                }
                dat->getVals()[length-1+half_length-i] = dat->getVals()[i-half_length];
            }
        }
    }
    else if (type=="butterworth"){
        dat=std::make_shared<vec> (hyper(length));
        dat->zero();
        dat->getVals()[0] = 1;
        if (lwc==0)
            iirButterworth(dat->getVals(),length,hwc,true,ho);
        else if (hwc==1)
            iirButterworth(dat->getVals(),length,lwc,false,ho);
        else
            iirBpButterworth(dat->getVals(),length,lwc,hwc,ho);
    }
    else if (type=="exponential"){
        dat=std::make_shared<vec> (hyper(length));
        for (int i=0; i<dat->getHyper()->getN123(); i++)
            dat->getVals()[i] = exp(-i*dt/sigma);
    }
    else if (type=="power"){
        dat=std::make_shared<vec> (hyper(length));
        dat->getVals()[0]=0;
        for (int i=1; i<dat->getHyper()->getN123(); i++)
            dat->getVals()[i] = std::pow(i*dt,pow);
    }
    else if (type=="bandpass"){
        if ((w1<0) || (w1>=w2) || (w2>w3) || (w3>=w4) || (w4>1))
            throw std::logic_error("The bandpass frequencies must obey 0 <= w1 < w2 < w3 < w4 <= 1.\n");
        ax T(length,0,dt);
        ax X(1,0,1);
        dat=std::make_shared<vec> (hyper(T));
        dat->zero();
        fxTransform fx(hyper(T,X));
        ax F = fx.getRange()->getAxis(1);
        std::shared_ptr<cvec> fxvec = std::make_shared<cvec >(hyper(F,X));
        fxvec->zero();
        std::complex<data_t> * pf = fxvec->getVals();
        data_t dw = 2*M_PI/T.n;
        data_t r;

        int iw1, iw2, iw3, iw4;
        data_t W1, W2, W3, W4;
        iw1=floor(w1*(F.n-1));
        iw2=floor(w2*(F.n-1));
        iw3=floor(w3*(F.n-1));
        iw4=floor(w4*(F.n-1));
        W1=w1*M_PI; W2=w2*M_PI; W3=w3*M_PI; W4=w4*M_PI;

        for (int i=iw1; i<=iw2; i++) {
            r = 0.5*(1-cos(M_PI*(i*dw-W1)/(W2-W1)));
            pf[i].real(r*cos(-i*dw*(T.n - 1)/2));
            pf[i].imag(r*sin(-i*dw*(T.n - 1)/2));
        }
        for (int i=iw2+1; i<=iw3; i++) {
            pf[i].real(cos(-i*dw*(T.n - 1)/2));
            pf[i].imag(sin(-i*dw*(T.n - 1)/2));
        }
        for (int i=iw3+1; i<=iw4; i++) {
            r = 0.5*(1+cos(M_PI*(i*dw-W3)/(W4-W3)));
            pf[i].real(r*cos(-i*dw*(T.n - 1)/2));
            pf[i].imag(r*sin(-i*dw*(T.n - 1)/2));
        }

        fx.inverse(false,dat,fxvec);
        dat->scale(1.0/dat->max());
    }
    else if (type=="spectrum"){
        std::shared_ptr<vec> input = read<data_t>(input_file, format);
        ax T(length,0,dt);
        ax X(1,0,1);
        dat=std::make_shared<vec> (hyper(T));  
        dat->zero();
        fxTransform fx(hyper(T,X));
        ax F = fx.getRange()->getAxis(1);
        std::shared_ptr<cvec> fxvec = std::make_shared<cvec >(hyper(F,X));

        const data_t * pin = input->getCVals();
        std::complex<data_t> * pf = fxvec->getVals();
        ax A = input->getHyper()->getAxis(1);
        data_t dw = 2*M_PI/T.n;
        data_t r, f;
        int imin, imax;

        for (int i=0; i<F.n; i++){
            f = F.d*i;
            imin=std::max(0, (int)floor((f-A.o)/A.d));
            imax=std::max(0, (int)round((f-A.o)/A.d));
            imin=std::min(imin, A.n-1);
            imax=std::min(imax, A.n-1);

            // interpolate the amplitude spectrum from the provided one
            r = (1-(f-imin*A.d)/A.d)*pin[imin] + (f-imin*A.d)/A.d*pin[imax];

            // add linear phase shift (* exp (-j.omega.dt.(Nt-1)/2)))
            pf[i].real(r*cos(-i*dw*(T.n - 1)/2));
            pf[i].imag(r*sin(-i*dw*(T.n - 1)/2));
        }

        fx.inverse(false,dat,fxvec);

    }
    else {
        throw std::logic_error("The requested type is not implemented.\n");
    }

// adjust the phase
    ax T;
    if (phase=="zero"){
        dat = zero_phase(dat);
        T.o=(-floor(nt/2)+shift)*dt; // time origin
    }
    else if (phase=="minimum"){
        dat = minimum_phase(dat,eps);
        T.o=shift*dt; // time origin
    }
    else{
        // std::clog << "The phase will be set to default.\n";
    }
        
//  set the hypercube
	T.d=dt; // (in sec)
	T.n=length; // number of samples in time
    hyper hyp(T);
    dat->setHyper(hyp);
     
    if (output_file!="none") write<data_t>(dat, output_file, format, datapath);

    return 0;
}