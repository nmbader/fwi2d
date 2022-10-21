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
        int nx = _domain.getN123()/T.n;
        int it1, it2, it3, it4;
        it1=floor((_t1-T.o)/T.d);
        it2=floor((_t2-T.o)/T.d);
        it3=ceil((_t3-T.o)/T.d);
        it4=ceil((_t4-T.o)/T.d);

        #pragma omp parallel for
        for (size_t j=0; j<nx; j++){
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
    "   t1,t2,t3,t4 - int - []:\n\t\ttime samples defining the taper of the filter if any.\n"
    "   solver - string - ['cgls']:\n\t\tlinear solver to use, 'sdls' or 'cgls'.\n"
    "   niter - int - [10]:\n\t\tnumber of iterations.\n"
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
    int fl=7, niter=10, t1=0, t2=0, t3=0, t4=0;
    data_t threshold=0;
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

        // time convolution operator with windowing
        convnd1d conv(*fs->getHyper(), ds, centered);
        window W(*fs->getHyper(),t1,t2,t3,t4);
        chainLOper op(&conv, &W);

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

        // compute the matched data
        op.forward(false, fs, ds);

        // save the filter and matched data
        memcpy(input->getVals()+s*nt*ntr, ds->getVals(), sizeof(data_t)*nt*ntr);
        memcpy(f->getVals()+s*Tf.n, fs->getVals(), sizeof(data_t)*Tf.n);

        delete sol;
    }
          
    if (output_file!="none") write<data_t>(input, output_file, format, datapath);
    if (filter_file!="none") write<data_t>(f, filter_file, format, datapath);
    if (obj_func_file!="none") write<data_t>(all_func, obj_func_file, format, datapath);

    return 0;
}