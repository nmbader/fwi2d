#include <unistd.h>

#include "seplib.h"
#include "param.hpp"
#include "operator.hpp"
#include "optimization.hpp"
#include "lsolver.hpp"
#include "IO.hpp"

typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;

// Window operator (cosine squared) to be applied to the matching filter
class window : public loper {
protected:
    int _t1,_t2,_t3,_t4; // Times defining the window function
public:
    window(){}
    ~window(){}
    window(const hyper &domain, const int t1, const int t2, const int t3, const int t4){
        successCheck((t1>=0) && (t1<=t2) && (t2<=t3) && (t3<=t4) && (t4 <= domain.getAxis(1).n),__FILE__,__LINE__,"Window samples must satisfy 0<=t1<=t2<=t3<=t4<=filter_length\n");
        _domain = domain;
        _range = domain;
        _t1=t1;
        _t2=t2;
        _t3=t3;
        _t4=t4;
    }
    window * clone() const {
        window * op = new window(_domain, _t1,_t2,_t3,_t4);
        return op;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){

        ax T = _domain.getAxis(1);
        int ntr = _domain.getN123()/T.n;

        #pragma omp parallel for
        for (size_t j=0; j<ntr; j++){
            data_t val=0, t=0;
            for (int i=0; i<_t1; i++) pdat[j*T.n+i] = add*pdat[j*T.n+i];
            for (int i=_t1; i<_t2; i++) {
                val = cos(0.5*M_PI*(_t2-i)/(_t2-_t1));
                pdat[j*T.n+i] = add*pdat[j*T.n+i] + pmod[j*T.n+i]*val*val;
            }
            for (int i=_t2; i<_t3; i++) {
                pdat[j*T.n+i] = add*pdat[j*T.n+i] + pmod[j*T.n+i];
            }
            for (int i=_t3; i<_t4; i++) {
                val = cos(0.5*M_PI*(i-_t3)/(_t4-_t3));
                pdat[j*T.n+i] = add*pdat[j*T.n+i] + pmod[j*T.n+i]*val*val;
            }
            for (int i=_t4; i<T.n; i++) pdat[j*T.n+i] = add*pdat[j*T.n+i];
        }
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        apply_forward(add,pdat,pmod);
    }
};

// Operator to apply Gaussian taper in the time domain (to be applied to the matching filter)
class temporal_taper : public loper {
protected:
    data_t _stddev; // Standard deviation of the taper in terms of fraction of filter length
    data_t _bias; // Bias to shift the taper in terms of fraction of filter length
public:
    temporal_taper(){}
    ~temporal_taper(){}
    temporal_taper(const hyper &domain, const data_t stddev, const data_t bias){
        successCheck( stddev>0 ,__FILE__,__LINE__,"Standard deviation must be positive\n");
        _domain = domain;
        _range = domain;
        _stddev=stddev;
        _bias = bias;
    }
    temporal_taper * clone() const {
        temporal_taper * op = new temporal_taper(_domain, _stddev, _bias);
        return op;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){

        ax T = _domain.getAxis(1);
        int ntr = _domain.getN123()/T.n;

        data_t factor = _stddev*(T.n-1);
        factor *= factor;
        data_t shift = _bias*(T.n-1);
        #pragma omp parallel for
        for (int itr=0; itr<ntr; itr++){
            for (int t=0; t<T.n; t++){
                pdat[itr*T.n+t] = add*pdat[itr*T.n+t] + pmod[itr*T.n+t] * exp(-0.5*(t-shift)*(t-shift)/factor);
            }
        }
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        apply_forward(add,pdat,pmod);
    } 
};

// Operator to apply Gaussian taper in the frequency domain (to be applied to the matching filter)
// This operator can be written as: d = F^-1.W.F.m
// F is FT transform (real to real vectors)
// W is a Gaussian taper
// F^-1 is inverse FT transform (real to real vectors)
class spectral_taper : public loper {
protected:
    data_t _stddev; // Standard deviation of the taper in terms of fraction of Nyquist frequency
    data_t _bias; // Bias to shift the taper in terms of fraction of Nyquist frequency
public:
    spectral_taper(){}
    ~spectral_taper(){}
    spectral_taper(const hyper &domain, const data_t stddev, const data_t bias){
        successCheck( stddev>0 ,__FILE__,__LINE__,"Standard deviation must be positive\n");
        _domain = domain;
        _range = domain;
        _stddev=stddev;
        _bias = bias;
    }
    spectral_taper * clone() const {
        spectral_taper * op = new spectral_taper(_domain, _stddev, _bias);
        return op;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){

        // F.m
        fxTransformR fx(_domain);
        std::shared_ptr<vecReg<data_t> > vec1 = std::make_shared<vecReg<data_t> >(*fx.getRange());
        vec1->zero();
        data_t * pvec1=vec1->getVals();
        fx.apply_forward(false,pmod,pvec1);

        // W.F.m
        int ntr = _domain.getN123()/_domain.getAxis(1).n;
        int nf = fx.getRange()->getAxis(1).n / 2;
        data_t df = fx.getRange()->getAxis(1).d;
        data_t (*pv1) [2*nf] = (data_t (*)[2*nf]) pvec1;
        data_t factor = _stddev*(nf-1);
        factor *= factor;
        data_t shift = _bias*(nf-1);
        #pragma omp parallel for
        for (int itr=0; itr<ntr; itr++){
            for (int f=0; f<nf; f++){
                pv1[itr][2*f] *= exp(-0.5*(f-shift)*(f-shift)/factor);
                pv1[itr][2*f+1] *= exp(-0.5*(f-shift)*(f-shift)/factor);
            }
        }

        // F^-1.W.F.m
        ifxTransformR ifx(_range);
        ifx.apply_forward(add,pvec1,pdat);
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        
        // F'^-1.d
        ifxTransformR ifx(_range);
        std::shared_ptr<vecReg<data_t> > vec1 = std::make_shared<vecReg<data_t> >(*ifx.getDomain());
        vec1->zero();
        data_t * pvec1=vec1->getVals();
        ifx.apply_adjoint(false,pvec1,pdat);

        // W'.F'^-1.d
        int ntr = _domain.getN123()/_domain.getAxis(1).n;
        int nf = ifx.getDomain()->getAxis(1).n / 2;
        data_t df = ifx.getDomain()->getAxis(1).d;
        data_t (*pv1) [2*nf] = (data_t (*)[2*nf]) pvec1;
        data_t factor = _stddev*(nf-1);
        factor *= factor;
        data_t shift = _bias*(nf-1);
        #pragma omp parallel for
        for (int itr=0; itr<ntr; itr++){
            for (int f=0; f<nf; f++){
                pv1[itr][2*f] *= exp(-0.5*(f-shift)*(f-shift)/factor);
                pv1[itr][2*f+1] *= exp(-0.5*(f-shift)*(f-shift)/factor);
            }
        }

        // F'.W'.F'^-1.d
        fxTransformR fx(_domain);
        fx.apply_adjoint(add,pmod,pvec1);
    }
};

void printdoc(){
    std::string doc = "\nDescription:\n"
    "   Match input to target data by least-squares filter.\n"
    "   Dimension 1 is for time. Dimension 2 is for traces within the same gather. Dimensions 3 and beyond are for separate gathers.\n"
    "   The number of output filters is equal to the number of elements in Dimensions 3 and beyond.\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H', target data as 'target=file.H', and output matched data as 'output=file.H' or '> file.H'. '<' and '>' are valid for SEPlib format only.\n"
    "\nParameters:\n"
    "   filter_half_length - int - [7]:\n\t\thalf length of the filters (number of samples). The full length is always 2*hl+1.\n"
    "   centered - bool - [1]:\n\t\t0 for causal convolution, 1 for centered (acausal) convolution.\n"
    "   filter - string - ['none']:\n\t\tname of file to output filters.\n"
    "   t1,t2,t3,t4 - int - [0,0,2*hl+1,2*hl+1]:\n\t\ttime samples defining the cosine squared taper of the filter (default no taper).\n"
    "   stddev_t - float - [0]:\n\t\tStandard deviation in fraction of filter length for the temporal Gaussian taper (default no taper).\n"
    "   bias_t - float - [0]:\n\t\tBias (shift) in fraction of filter length for the temporal Gaussian taper.\n"
    "   stddev_f - float - [0]:\n\t\tStandard deviation in fraction of Nyquist frequency for the spectral Gaussian taper (default no taper).\n"
    "   bias_f - float - [0]:\n\t\tBias (shift) in fraction of Nyquist frequency for the spectral Gaussian taper.\n"
    "   solver - string - ['cgls']:\n\t\tlinear solver to use, 'sdls' or 'cgls'.\n"
    "   niter - int - [5]:\n\t\tnumber of iterations.\n"
    "   threshold - float - [0]:\n\t\tconvergence threshold to stop the solver.\n"
    "   verbose - bool - [1]:\n\t\tprint info during iterations.\n"
    "   obj_func - string - ['none']:\n\t\tname of file to output the objective functions.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries.\n"
    "\nExample:\n"
    "   MATCHING.x < infile.H target=target.H filter_half_length=11 niter=20 filter=filters.H obj_func=func.H > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);

    std::string input_file="in", output_file="out", target_file="none", filter_file="none", solver="cgls", obj_func_file="none", datapath="none";
    int fl=7, niter=5, t1=0, t2=0, t3=0, t4=0;
    data_t stddev_t=0, bias_t=0, stddev_f=0, bias_f=0, threshold=0;
    bool format=0, centered=true, verbose=true;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
	readParam<std::string>(argc, argv, "target", target_file);
	readParam<std::string>(argc, argv, "filter", filter_file);
	readParam<std::string>(argc, argv, "solver", solver);
    readParam<std::string>(argc, argv, "obj_func", obj_func_file);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<int>(argc, argv, "filter_half_length", fl);
    readParam<int>(argc, argv, "niter", niter);
    t3=2*fl+1;
    t4=t3;
    readParam<int>(argc, argv, "t1", t1);
    readParam<int>(argc, argv, "t2", t2);
    readParam<int>(argc, argv, "t3", t3);
    readParam<int>(argc, argv, "t4", t4);
    readParam<data_t>(argc, argv, "stddev_t", stddev_t);
    readParam<data_t>(argc, argv, "bias_t", bias_t);
    readParam<data_t>(argc, argv, "stddev_f", stddev_f);
    readParam<data_t>(argc, argv, "bias_f", bias_f);
    readParam<data_t>(argc, argv, "threshold", threshold);
    readParam<bool>(argc, argv, "centered", centered);
    readParam<bool>(argc, argv, "verbose", verbose);
    readParam<bool>(argc, argv, "format", format);

    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");
    successCheck(target_file!="none",__FILE__,__LINE__,"Target file is not provided\n");
    
    std::shared_ptr<vec> input = read<data_t>(input_file,format);
    std::shared_ptr<vec> target = read<data_t>(target_file,format);

    successCheck(input->getN123()==target->getN123(),__FILE__,__LINE__,"Input and target must have the same number of samples\n");
    successCheck(input->getHyper()->getAxis(1).n==target->getHyper()->getAxis(1).n,__FILE__,__LINE__,"Input and target must have the same number of samples along the first axis\n");

    int ntr=1, ns=1, ndim=input->getHyper()->getNdim();
    ax T = input->getHyper()->getAxis(1);
    int nt=T.n;
    ax Tf = T;
    Tf.n = 2*fl+1;
    Tf.o=0;
    if (centered) Tf.o = -fl*Tf.d;
    if (ndim>1) {
        ntr=input->getHyper()->getAxis(2).n;
        ns=input->getN123()/(ntr*nt);
    }

    std::shared_ptr<vec> f = std::make_shared<vec> (hyper(Tf,ax(ns,0,1)));
    f->zero();
    hyper hyp(T,ax(ntr,0,1));
    std::shared_ptr<vec> all_func = std::make_shared<vec> (hyper(niter+1,ns));
    all_func->set(-1);

    // filters preconditioning operator P = W.TT.ST (temporal and spectral Gaussian tapers and temporal cosine squared taper)
    loper * P = nullptr;
    {
        window W(hyper(Tf),t1,t2,t3,t4);
        if (stddev_f>0){
            fprintf(stderr,"Spectral taper activated\n");
            spectral_taper ST(hyper(Tf), stddev_f, bias_f);
            if (stddev_t>0){
                fprintf(stderr,"Temporal taper activated\n");
                temporal_taper TT(hyper(Tf), stddev_t, bias_t);
                chainLOper TT_ST(&TT, &ST);
                P = new chainLOper(&W, &TT_ST);
            }
            else{
                P = new chainLOper(&W, &ST);
            }
        }
        else if (stddev_t>0){
            fprintf(stderr,"Temporal taper activated\n");
            temporal_taper TT(hyper(Tf), stddev_t, bias_t);
            P = new chainLOper(&W, &TT);
        }
        else P = W.clone();
    }

    for (int s=0; s<ns; s++){
        fprintf(stderr,"\n==========================\n Start processing gather %d\n==========================\n",s);

        // copy a single gather
        std::shared_ptr<vec> ds = std::make_shared<vec> (hyp);
        std::shared_ptr<vec> ts = std::make_shared<vec> (hyp);
        memcpy(ds->getVals(), input->getVals()+s*nt*ntr,sizeof(data_t)*nt*ntr);
        memcpy(ts->getVals(), target->getVals()+s*nt*ntr,sizeof(data_t)*nt*ntr);

        // initialize a single filter to a Dirac
        std::shared_ptr<vec> fs = std::make_shared<vec> (hyper(Tf));
        fs->zero();
        if (centered) fs->getVals()[fl]=1;
        else fs->getVals()[0]=1;

        // time convolution operator with preconditioning
        convnd1d conv(*fs->getHyper(), ds, centered);
        chainLOper op(&conv, P);

        // set the least-squares problem
        llsq prob(&op, fs, ts);

        // Set the linear solver
        lsolver * sol;
        if (solver == "sdls") sol = new sdls(niter,threshold);
        else sol = new cgls(niter,threshold);

        // Invert the data
        sol->run(&prob,verbose);

        // copy the objective function
        memcpy(all_func->getVals()+(niter+1)*s, sol->_func.data(), sizeof(data_t)*sol->_func.size());

        // compute the final filter
        std::shared_ptr<vec> temp = fs->clone();
        P->forward(false, fs, temp);
        fs = temp;

        // compute the matched data
        conv.forward(false, fs, ds);

        // save the filter and matched data
        memcpy(input->getVals()+s*nt*ntr, ds->getVals(), sizeof(data_t)*nt*ntr);
        memcpy(f->getVals()+s*Tf.n, fs->getVals(), sizeof(data_t)*Tf.n);

        delete sol;
    }

    delete P;
          
    if (output_file!="none") write<data_t>(input, output_file, format, datapath);
    if (filter_file!="none") write<data_t>(f, filter_file, format, datapath);
    if (obj_func_file!="none") write<data_t>(all_func, obj_func_file, format, datapath);

    return 0;
}