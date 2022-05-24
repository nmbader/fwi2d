#include <string.h>
#include "we_op.hpp"
#include "bsplines.hpp"
#include "nlsolver.hpp"
#include "IO.hpp"
#include "seplib.h"

typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;

// operator that maps from several models (one for each source) to the physical space with damping (cosine squared) around the sources location 
// m = W.ms where ms is the extended model across sources, 'W' is the row operator of normalized cosine squared weighting operator for each source 's'
class sWeighting : public loper {
protected:
    data_t _dpower; // damping strength for cosine^2
    data_t _dwidthx; // damping width (radius)
    data_t _dwidthz;
    std::shared_ptr<vec> _nW; // normalization for the weighting (sum of all weighting operators)
    std::vector<std::vector<data_t> > _sxz; // sources location

public:
    sWeighting(){}
    ~sWeighting(){}
    sWeighting(const hypercube<data_t> &domain,const std::vector<std::vector<data_t> > &sxz, const data_t dwidthx, const data_t dwidthz, const data_t dpower){
        std::vector<ax > axes = domain.getAxes();
        successCheck(axes.size()>=3,__FILE__,__LINE__,"The domain must contain at least 3 axes\n");
        int ns=domain.getN123()/(axes[0].n*axes[1].n*axes[2].n);
        ax S(ns,0,1);
        successCheck(ns==sxz.size(),__FILE__,__LINE__,"The domain size must be consistent with the number of sources\n");
        
        _domain = hyper(axes[0],axes[1],axes[2],S);
        _range = hyper(axes[0],axes[1],axes[2]);
        _dwidthx = dwidthx;
        _dwidthz = dwidthz;
        _dpower = dpower;
        _sxz=sxz;

        // build the normalization weights
        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        _nW = std::make_shared<vec> (vec(hyper(Z,X)));
        _nW->zero();
        data_t (*pw) [Z.n] = (data_t (*)[Z.n]) _nW->getVals();

        for (int s=0; s<ns; s++){

            // normalization weight
            int ixmin = ceil((sxz[s][0]-dwidthx-X.o)/X.d);
            int ixmax = ceil((sxz[s][0]+dwidthx-X.o)/X.d);
            int izmin = ceil((sxz[s][1]-dwidthz-Z.o)/Z.d);
            int izmax = ceil((sxz[s][1]+dwidthz-Z.o)/Z.d);

            ixmin=std::max(0,ixmin);
            ixmax=std::min(X.n,ixmax);
            izmin=std::max(0,izmin);
            izmax=std::min(Z.n,izmax);

            for (int ix=ixmin; ix<ixmax; ix++){
                data_t x = X.o + ix*X.d;
                data_t rx = (x-_sxz[s][0])/_dwidthx;
                rx*=rx;
                for (int iz=izmin; iz<izmax; iz++){
                    data_t z = Z.o + iz*Z.d;
                    data_t rz = (z-_sxz[s][1])/_dwidthz;
                    rz*=rz;
                    data_t d = sqrt(rx+rz);
                    data_t val = 1;
                    if (d<1) val = cos(_dpower*0.5*M_PI*(1-d));
                    pw[ix][iz] += val*val;
                }
            }

            for (int ix=0; ix<ixmin; ix++){
                for (int iz=0; iz<Z.n; iz++) pw[ix][iz] += 1;
            }
            for (int ix=ixmax; ix<X.n; ix++){
                for (int iz=0; iz<Z.n; iz++) pw[ix][iz] += 1;
            }
            for (int ix=ixmin; ix<ixmax; ix++){
                for (int iz=0; iz<izmin; iz++) pw[ix][iz] += 1;
                for (int iz=izmax; iz<Z.n; iz++) pw[ix][iz] += 1;
            }
        }
    }

    sWeighting * clone() const {
        sWeighting * op = new sWeighting(_domain, _sxz, _dwidthx, _dwidthz, _dpower);
        return op;
    }    
    int getNs(){return _sxz.size();}
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){

        if (!add) memset(pdat,0,_range.getN123()*sizeof(data_t));

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        int ns = _sxz.size();
        const data_t (*pm) [C.n][X.n][Z.n] = (const data_t (*)[C.n][X.n][Z.n]) pmod;
        data_t (*pd) [X.n][Z.n] = (data_t (*)[X.n][Z.n]) pdat;
        const data_t (*pw) [Z.n] = (const data_t (*)[Z.n]) _nW->getVals();

        for (int s=0; s<ns; s++){
            int ixmin = ceil((_sxz[s][0]-_dwidthx-X.o)/X.d);
            int ixmax = ceil((_sxz[s][0]+_dwidthx-X.o)/X.d);
            int izmin = ceil((_sxz[s][1]-_dwidthz-Z.o)/Z.d);
            int izmax = ceil((_sxz[s][1]+_dwidthz-Z.o)/Z.d);
            ixmin=std::max(0,ixmin);
            ixmax=std::min(X.n,ixmax);
            izmin=std::max(0,izmin);
            izmax=std::min(Z.n,izmax);
            for (int c=0; c<C.n; c++){
                #pragma omp parallel for
                for (int ix=ixmin; ix<ixmax; ix++){
                    data_t x = X.o + ix*X.d;
                    data_t rx = (x-_sxz[s][0])/_dwidthx;
                    rx*=rx;
                    for (int iz=izmin; iz<izmax; iz++){
                        data_t z = Z.o + iz*Z.d;
                        data_t rz = (z-_sxz[s][1])/_dwidthz;
                        rz*=rz;
                        data_t d = sqrt(rx+rz);
                        data_t val = 1;
                        if (d<1) val = cos(_dpower*0.5*M_PI*(1-d));
                        pd[c][ix][iz] += pm[s][c][ix][iz]*val*val/pw[ix][iz];
                    }
                }

                #pragma omp parallel for
                for (int ix=0; ix<ixmin; ix++){
                    for (int iz=0; iz<Z.n; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]/pw[ix][iz];
                }
                #pragma omp parallel for
                for (int ix=ixmax; ix<X.n; ix++){
                    for (int iz=0; iz<Z.n; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]/pw[ix][iz];
                }
                #pragma omp parallel for
                for (int ix=ixmin; ix<ixmax; ix++){
                    for (int iz=0; iz<izmin; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]/pw[ix][iz];
                    for (int iz=izmax; iz<Z.n; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]/pw[ix][iz];
                }
            }            
        }
    }

    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        int ns = _sxz.size();
        data_t (*pm) [C.n][X.n][Z.n] = (data_t (*)[C.n][X.n][Z.n]) pmod;
        const data_t (*pd) [X.n][Z.n] = (const data_t (*)[X.n][Z.n]) pdat;
        const data_t (*pw) [Z.n] = (const data_t (*)[Z.n]) _nW->getVals();

        for (int s=0; s<ns; s++){
            int ixmin = ceil((_sxz[s][0]-_dwidthx-X.o)/X.d);
            int ixmax = ceil((_sxz[s][0]+_dwidthx-X.o)/X.d);
            int izmin = ceil((_sxz[s][1]-_dwidthz-Z.o)/Z.d);
            int izmax = ceil((_sxz[s][1]+_dwidthz-Z.o)/Z.d);
            ixmin=std::max(0,ixmin);
            ixmax=std::min(X.n,ixmax);
            izmin=std::max(0,izmin);
            izmax=std::min(Z.n,izmax);
            for (int c=0; c<C.n; c++){
                for (int ix=ixmin; ix<ixmax; ix++){
                    data_t x = X.o + ix*X.d;
                    data_t rx = (x-_sxz[s][0])/_dwidthx;
                    rx*=rx;
                    for (int iz=izmin; iz<izmax; iz++){
                        data_t z = Z.o + iz*Z.d;
                        data_t rz = (z-_sxz[s][1])/_dwidthz;
                        rz*=rz;
                        data_t d = sqrt(rx+rz);
                        data_t val = 1;
                        if (d<1) val = cos(_dpower*0.5*M_PI*(1-d));
                        pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val*val/pw[ix][iz];
                    }
                }

                for (int ix=0; ix<ixmin; ix++){
                    for (int iz=0; iz<Z.n; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]/pw[ix][iz];
                }
                for (int ix=ixmax; ix<X.n; ix++){
                    for (int iz=0; iz<Z.n; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]/pw[ix][iz];
                }
                for (int ix=ixmin; ix<ixmax; ix++){
                    for (int iz=0; iz<izmin; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]/pw[ix][iz];
                    for (int iz=izmax; iz<Z.n; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]/pw[ix][iz];
                }
            }
        }
    }
};

// operator that maps from several models (one for each source) to the physical space with damping (cosine squared) around the sources location, and then spray the difference (initial - collapsed) across sources 
// Op = (I-S.W).ms where ms is the extended model across sources, 'W' is the row operator of normalized cosine squared weighting operator for each source 's', 'S' is a spraying operator
class model_extension : public loper {
protected:
    sWeighting * _sW;

public:
    model_extension(){}
    ~model_extension(){delete _sW;}
    model_extension(sWeighting * sW){
        _sW = sW->clone();
        _domain = *sW->getDomain();
        _range = _domain;
    }
    model_extension * clone() const {
        model_extension * op = new model_extension(_sW);
        return op;
    }
    sWeighting * getsW(){return _sW;}
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        int ns = _sW->getNs();
        const data_t (*pm) [C.n][X.n][Z.n] = (const data_t (*)[C.n][X.n][Z.n]) pmod;
        data_t (*pd) [C.n][X.n][Z.n] = (data_t (*)[C.n][X.n][Z.n]) pdat;

        data_t * vec = new data_t[C.n*X.n*Z.n];
        memset(vec, 0, C.n*X.n*Z.n*sizeof(data_t));
        data_t (*pv) [X.n][Z.n] = (data_t (*)[X.n][Z.n]) vec;

        // W.m
        _sW->apply_forward(false,pmod,vec);

        // (I-S.W).m
        for (int s=0; s<ns; s++){
            for (int c=0; c<C.n; c++){
                #pragma omp parallel for
                for (int ix=0; ix<X.n;ix++){
                    for (int iz=0; iz<Z.n;iz++) pd[s][c][ix][iz] = add*pd[s][c][ix][iz] + pm[s][c][ix][iz] - pv[c][ix][iz];
                }
            }
        }
        delete [] vec;
    }

    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        int ns = _sW->getNs();
        data_t (*pm) [C.n][X.n][Z.n] = (data_t (*)[C.n][X.n][Z.n]) pmod;
        const data_t (*pd) [C.n][X.n][Z.n] = (const data_t (*)[C.n][X.n][Z.n]) pdat;

        data_t * vec = new data_t[C.n*X.n*Z.n];
        memset(vec, 0, C.n*X.n*Z.n*sizeof(data_t));
        data_t (*pv) [X.n][Z.n] = (data_t (*)[X.n][Z.n]) vec;

        // -S'.d and I.d
        for (int s=0; s<ns; s++){
            for (int c=0; c<C.n; c++){
                #pragma omp parallel for
                for (int ix=0; ix<X.n;ix++){
                    for (int iz=0; iz<Z.n;iz++) {
                        pv[c][ix][iz] -= pd[s][c][ix][iz];
                        pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[s][c][ix][iz];
                    }
                }
            }
        }

        // (I-W'.S').d
        _sW->apply_adjoint(true,pmod,vec);
        delete [] vec;
    }
};

// soft clip of the elastic model and the Vs/Vp ratio for the extended model
class emodelSoftClipExt : public nloper {
    emodelSoftClip * _S;
public:
    emodelSoftClipExt(){}
    ~emodelSoftClipExt(){delete _S;}
    emodelSoftClipExt(const hypercube<data_t> &domain, emodelSoftClip * S){
        successCheck((domain.getNdim()>=3) && (domain.getAxis(3).n>=3),__FILE__,__LINE__,"The domain must be at least 3D with 3rd dimension containing at least 3 fields\n");
        _domain = domain;
        _range = domain;
        hyper hyp(domain.getAxis(1),domain.getAxis(2),domain.getAxis(3));
        _S = S->clone();
    }
    emodelSoftClipExt * clone() const {
        emodelSoftClipExt * op = new emodelSoftClipExt(_domain, _S);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){
        int n=_domain.getAxis(1).n*_domain.getAxis(2).n*_domain.getAxis(3).n;
        int ns=_domain.getN123()/n;
        for (int s=0; s<ns; s++) _S->apply_forward(add,pmod+s*n,pdat+s*n);
    }
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat){
        int n=_domain.getAxis(1).n*_domain.getAxis(2).n*_domain.getAxis(3).n;
        int ns=_domain.getN123()/n;
        for (int s=0; s<ns; s++) _S->apply_jacobianT(add,pmod+s*n,pmod0+s*n,pdat+s*n);
    }
};

// Executable to run 2D FWI

int main(int argc, char **argv){

int rank=0, size=0;
#ifdef ENABLE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    fprintf (stderr,"\n====================\nSize of MPI communicator = %d ; current rank = %d\n====================\n",size,rank);
#endif
    
    initpar(argc,argv);

// Read parameters for wave propagation and inversion
    param par;
    readParameters(argc, argv, par);
    int verbose=par.verbose;
    if (rank>0) par.verbose=0;
    par.device+=rank;

// Set the maximum number of threads
    if (par.nthreads>0) omp_set_num_threads(par.nthreads);

// Read inputs/outputs files
    std::string source_file="none", model_file="none", data_file="none", output_file="none", ioutput_file="none", obj_func_file="none";
    readParam<std::string>(argc, argv, "source", source_file);
    readParam<std::string>(argc, argv, "model", model_file);
    readParam<std::string>(argc, argv, "data", data_file);
    readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "ioutput", ioutput_file);
    readParam<std::string>(argc, argv, "obj_func", obj_func_file);

// Read some parameters
    data_t dwidthx=0, dwidthz=0, dpower=0;
    readParam<data_t>(argc, argv, "damping_widthx", dwidthx);
    readParam<data_t>(argc, argv, "damping_widthz", dwidthz);
    readParam<data_t>(argc, argv, "damping_power", dpower);

    successCheck(source_file!="none",__FILE__,__LINE__,"Source wavelet is not provided\n");
    successCheck(model_file!="none",__FILE__,__LINE__,"Earth model is not provided\n");
    successCheck(data_file!="none",__FILE__,__LINE__,"Data to be inverted is not provided\n");

    std::shared_ptr<vec> src = read<data_t>(source_file, par.format);
    std::shared_ptr<vec> data = read<data_t>(data_file, par.format);
    std::shared_ptr<vec> model = read<data_t>(model_file, par.format);
    hyper hyp0 = *model->getHyper();

    std::shared_ptr<vec> gmask = nullptr;
    std::shared_ptr<vec> w = nullptr;
    std::shared_ptr<vec> invDiagH = nullptr;
    std::shared_ptr<vec> prior = nullptr;
    if (par.mask_file!="none") gmask = read<data_t>(par.mask_file, par.format);
    if (par.weights_file!="none") w = read<data_t>(par.weights_file, par.format);
    if (par.inverse_diagonal_hessian_file!="none") invDiagH = read<data_t>(par.inverse_diagonal_hessian_file, par.format);

// Analyze the inputs and parameters and modify if necessary
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0);
    analyzeBsplines(*model->getHyper(),par);
    analyzeNLInversion(par);
    par.sextension=true;

// ----------------------------------------------------------------------------------------//
// Extend model along sources, extend mask and inverse Hessian if necessary
// ----------------------------------------------------------------------------------------//
    std::vector<ax > axes = model->getHyper()->getAxes();
    int n=model->getN123();
{
    ax S(par.ns,0,1);
    axes.push_back(S);
    std::shared_ptr<vec> tmp = std::make_shared<vec>(vec(hyper(axes)));
    for (int s=0; s<par.ns; s++) memcpy(tmp->getVals()+s*n,model->getVals(),model->getN123()*sizeof(data_t));
    model=tmp->clone();
}
    if (par.mask_file != "none"){
        if (gmask->getN123() == n){
            std::shared_ptr<vec> tmp = std::make_shared<vec>(vec(hyper(axes)));
            for (int s=0; s<par.ns; s++) memcpy(tmp->getVals()+s*n,gmask->getVals(),gmask->getN123()*sizeof(data_t));
            gmask=tmp;
        }
        else successCheck(gmask->getN123()==n*par.ns,__FILE__,__LINE__,"The extended gradient mask has an incorrect size\n");
    }
    if (par.inverse_diagonal_hessian_file != "none"){
        if (invDiagH->getN123() == n){
            std::shared_ptr<vec> tmp = std::make_shared<vec>(vec(hyper(axes)));
            for (int s=0; s<par.ns; s++) memcpy(tmp->getVals()+s*n,invDiagH->getVals(),invDiagH->getN123()*sizeof(data_t));
            invDiagH=tmp;
        }
        else successCheck(gmask->getN123()==n*par.ns,__FILE__,__LINE__,"The extended inverse Hessian has an incorrect size\n");
    }

    sWeighting * sW = new sWeighting(*model->getHyper(), par.sxz, dwidthx, dwidthz, dpower);
    model_extension * E = new model_extension(sW);
    delete sW;
  
// ----------------------------------------------------------------------------------------//
// Build model precon for model extension along sources with B-splines included
// ----------------------------------------------------------------------------------------//
    std::shared_ptr<vec> bsmodel = model;
    std::shared_ptr<vec> bsmask = gmask;
    std::shared_ptr<vec> bsinvDiagH = invDiagH;
    std::shared_ptr<vec> bsprior = prior;
    loper * BD;

if (par.bsplines)
{
    std::vector<data_t> kx;
    std::vector<data_t> kz;
    setKnot(kx,par.bs_controlx,par.bs_mx);
    setKnot(kz,par.bs_controlz,par.bs_mz);

    bsfillin F(*model->getHyper(),par.bs_controlx,par.bs_controlz);
    bsmodel = std::make_shared<vec>(vec(*F.getRange()));
    bsmodel->zero();
    F.apply_forward(false, model->getVals(), bsmodel->getVals());
    if (par.mask_file != "none") {
        bsmask = std::make_shared<vec>(vec(*F.getRange()));
        bsmask->zero();
        F.apply_forward(false, gmask->getVals(), bsmask->getVals());
    }
    if (par.inverse_diagonal_hessian_file != "none") {
        bsinvDiagH = std::make_shared<vec>(vec(*F.getRange()));
        bsinvDiagH->zero();
        F.apply_forward(false, invDiagH->getVals(), bsinvDiagH->getVals());
    }
    if (prior != nullptr) {
        bsprior = std::make_shared<vec>(vec(*F.getRange()));
        bsprior->zero();
        F.apply_forward(false, prior->getVals(), bsprior->getVals());
    }
   
    duplicate D(*F.getRange(),par.bs_mx,par.bs_mz);
    bsplines3 B(*D.getRange(),*model->getHyper(),kx,kz);
    BD = new chainLOper(&B,&D);
}
 
// ----------------------------------------------------------------------------------------//
// Build model precon if soft clipping is activated
// ----------------------------------------------------------------------------------------//
    emodelSoftClipExt * S;
    if (par.soft_clip) {
        emodelSoftClip S0(hyp0, par.vpmin, par.vpmax, par.vsmin, par.vsmax, par.rhomin, par.rhomax, 1/sqrt(2.00001), 9, 9);
        S = new emodelSoftClipExt(*model->getHyper(), &S0);
        if (par.verbose>0) fprintf(stderr,"Soft clipping is added to the inversion. It overrides the hard clipping\n");
    }
// ----------------------------------------------------------------------------------------//


    if (rank>0) par.verbose=verbose;
    nloper * op = nullptr;
    nl_we_op * L;
    if (par.nmodels==2) L=new nl_we_op_a(*model->getHyper(),allsrc,par);
    else if (par.nmodels==3 && !par.acoustic_elastic) L=new nl_we_op_e(*model->getHyper(),allsrc,par);
    else if (par.nmodels==3 && par.acoustic_elastic) L=new nl_we_op_ae(*model->getHyper(),allsrc,par);
    else if (par.nmodels==5) L=new nl_we_op_vti(*model->getHyper(),allsrc,par);

    if (rank>0) par.verbose=0;

    if (par.bsplines)
    {
        if (par.soft_clip) op  = new chainNLOper(S,BD);
        else op = BD->clone();
    }
    else
    {
        if (par.soft_clip) op = S->clone();
    }
    if (par.bsplines) delete BD;
    if (par.soft_clip) delete S;

    nloper * D = nullptr;
    
    if (op==nullptr) D = E->clone();
    else D = new chainNLOper(E, op);

    nlls_fwi_reg * prob = new nlls_fwi_reg(L, D, bsmodel, data, par.lambda, bsprior, op, bsmask, w);

    lsearch * ls;
    if (par.lsearch=="weak_wolfe") ls = new weak_wolfe(par.ls_c1, par.ls_a0, par.ls_a1, par.ls_version);
    else if(par.lsearch=="strong_wolfe") ls = new strong_wolfe(par.ls_c1, par.ls_c2, par.ls_a0, par.ls_a1, par.ls_max_step, par.ls_version);
    else ls = new regular_wolfe(par.ls_c1, par.ls_c2, par.ls_a0, par.ls_a1, par.ls_max_step, par.ls_version);

    nlsolver * solver;
    if (par.nlsolver=="nlsd") solver = new nlsd(par.niter, par.max_trial, par.threshold, ls); 
    else if (par.nlsolver=="nlcg") solver = new nlcg(par.niter, par.max_trial, par.threshold, ls); 
    else if (par.nlsolver=="bfgs") solver = new bfgs(par.niter, par.max_trial, par.threshold, ls); 
    else solver = new lbfgs(par.niter, par.max_trial, par.threshold, ls, bsinvDiagH, par.lbfgs_m); 
    
    solver->run(prob, par.verbose>0, ioutput_file, par.isave, par.format, par.datapath);
    
    if (D != nullptr) delete D;
    if (op != nullptr)
    {
        op->forward(false, bsmodel, model);
        delete op;
    }
    
    std::shared_ptr<vec> fmodel = std::make_shared<vec>(vec(*E->getsW()->getRange()));
    fmodel->zero();
    E->getsW()->forward(false,model,fmodel);

    if (rank==0 && output_file!="none") {
        write<data_t>(fmodel, output_file, par.format, par.datapath);
        write<data_t>(model, output_file+".ext", par.format, par.datapath);
    }
    if (rank==0 && obj_func_file!="none") {
        std::shared_ptr<vec > func = std::make_shared<vec > (hyper(solver->_func.size()));
        memcpy(func->getVals(), solver->_func.data(), solver->_func.size()*sizeof(data_t));
        write(func,obj_func_file, par.format, par.datapath);
        if (D != nullptr){
            std::shared_ptr<vec > dfunc = std::make_shared<vec > (hyper(prob->_dfunc.size()));
            std::shared_ptr<vec > mfunc = std::make_shared<vec > (hyper(prob->_mfunc.size()));
            memcpy(dfunc->getVals(), prob->_dfunc.data(), prob->_dfunc.size()*sizeof(data_t));
            memcpy(mfunc->getVals(), prob->_mfunc.data(), prob->_mfunc.size()*sizeof(data_t));
            write(dfunc,obj_func_file+".d", par.format, par.datapath);
            write(mfunc,obj_func_file+".m", par.format, par.datapath);
        }
    }

    delete E;
    delete L;
    delete prob;

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif

return 0;
}