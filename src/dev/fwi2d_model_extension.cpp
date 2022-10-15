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

// operator that takes and extended model (one for each source) and damp (cosine squared) around the sources location while borrowing from other sources location
// m = W.ms where ms is the extended model across sources, 'W' is a non-diagonal operator of normalized cosine squared weighting for each source 's'
class sWeighting : public loper {
protected:
    data_t _dpower; // damping strength for cosine^2
    data_t _dwidthx; // damping width (radius)
    data_t _dwidthz;
    std::vector<std::vector<data_t> > _sxz; // sources location

public:
    sWeighting(){}
    ~sWeighting(){}
    sWeighting(const hypercube<data_t> &domain,const std::vector<std::vector<data_t> > &sxz, const data_t dwidthx, const data_t dwidthz, const data_t dpower){
        std::vector<ax > axes = domain.getAxes();
        successCheck(axes.size()>=3,__FILE__,__LINE__,"The domain must contain at least 3 axes\n");
        int ns=domain.getN123()/(axes[0].n*axes[1].n*axes[2].n);
        successCheck(ns>=2,__FILE__,__LINE__,"The number of sources in the extension must be at least 2\n");
        ax S(ns,0,1);
        successCheck(ns==sxz.size(),__FILE__,__LINE__,"The domain size must be consistent with the number of sources\n");
        
        _domain = hyper(axes[0],axes[1],axes[2],S);
        _range = _domain;
        _dwidthx = dwidthx;
        _dwidthz = dwidthz;
        _dpower = dpower;
        _sxz=sxz;
    }

    sWeighting * clone() const {
        sWeighting * op = new sWeighting(_domain, _sxz, _dwidthx, _dwidthz, _dpower);
        return op;
    }    
    int getNs(){return _sxz.size();}
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){

        data_t * temp;
        if (add) {
            temp = new data_t [_domain.getN123()];
            memcpy(temp, pmod, _domain.getN123()*sizeof(data_t));
        }
        else {
            memcpy(pdat, pmod, _domain.getN123()*sizeof(data_t));
            temp = pdat;
        }

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        int ns = _sxz.size();
        const data_t (*pm) [C.n][X.n][Z.n] = (const data_t (*)[C.n][X.n][Z.n]) pmod;
        data_t (*pd) [C.n][X.n][Z.n] = (data_t (*)[C.n][X.n][Z.n]) temp;

        // loop over source-extended model 
        for (int s=0; s<ns; s++){

            // find the box surrounding the source location
            int ixmin = ceil((_sxz[s][0]-_dwidthx-X.o)/X.d);
            int ixmax = ceil((_sxz[s][0]+_dwidthx-X.o)/X.d);
            int izmin = ceil((_sxz[s][1]-_dwidthz-Z.o)/Z.d);
            int izmax = ceil((_sxz[s][1]+_dwidthz-Z.o)/Z.d);
            ixmin=std::max(0,ixmin);
            ixmax=std::min(X.n,ixmax);
            izmin=std::max(0,izmin);
            izmax=std::min(Z.n,izmax);

            // damp the model for source s using an ellispsis cosine squared and borrow model from the other sources
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
                        data_t val = 1.0;
                        if (d<1) val = cos(_dpower*0.5*M_PI*(1-d));
                        val *= val;

                        // damp with ellipsis cosine squared
                        pd[s][c][ix][iz] *= val;

                        // borrow a weighted averaged model from other sources
                        data_t sum = 0;
                        for (int so=0; so<ns; so++) sum += pm[so][c][ix][iz];
                        sum -= pm[s][c][ix][iz];
                        pd[s][c][ix][iz] += sum*(1-val)/(ns-1);
                    }
                }
            }
        }

        if (add) {
            #pragma omp parallel for
            for (int i=0; i<_domain.getN123(); i++) pdat[i] += temp[i];
            delete [] temp;
        }
    }

    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

        data_t * temp;
        if (add) {
            temp = new data_t [_domain.getN123()];
            memcpy(temp, pdat, _domain.getN123()*sizeof(data_t));
        }
        else {
            memcpy(pmod, pdat, _domain.getN123()*sizeof(data_t));
            temp = pmod;
        }

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        int ns = _sxz.size();
        data_t (*pm) [C.n][X.n][Z.n] = (data_t (*)[C.n][X.n][Z.n]) temp;
        const data_t (*pd) [C.n][X.n][Z.n] = (const data_t (*)[C.n][X.n][Z.n]) pdat;

        // set the output to zero in the source regions
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
                    for (int iz=izmin; iz<izmax; iz++){
                        pm[s][c][ix][iz] = 0;
                    }
                }
            }
        }

        // loop over source-extended model 
        for (int s=0; s<ns; s++){

            // find the box surrounding the source location
            int ixmin = ceil((_sxz[s][0]-_dwidthx-X.o)/X.d);
            int ixmax = ceil((_sxz[s][0]+_dwidthx-X.o)/X.d);
            int izmin = ceil((_sxz[s][1]-_dwidthz-Z.o)/Z.d);
            int izmax = ceil((_sxz[s][1]+_dwidthz-Z.o)/Z.d);
            ixmin=std::max(0,ixmin);
            ixmax=std::min(X.n,ixmax);
            izmin=std::max(0,izmin);
            izmax=std::min(Z.n,izmax);

            // damp the model for source s using an ellispsis cosine squared and spread model to the other sources
            #pragma omp parallel for
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
                        data_t val = 1.0;
                        if (d<1) val = cos(_dpower*0.5*M_PI*(1-d));
                        val *= val;

                        // spread a weighted model to the other sources
                        for (int so=0; so<ns; so++) pm[so][c][ix][iz] += pd[s][c][ix][iz]*(1-val)/(ns-1);
                        pm[s][c][ix][iz] += pd[s][c][ix][iz]*(val - (1-val)/(ns-1));
                    }
                }
            }
        }

        if (add) {
            #pragma omp parallel for
            for (int i=0; i<_domain.getN123(); i++) pmod[i] += temp[i];
            delete [] temp;
        }
    }
};

// operator that outputs the "average" of an extended model along sources (the extension is assumed along the last axis) 
// the average is taken with sqrt(nb of sources) to maintain the adjointness of the adjoint
class sSum : public loper {
public:
    sSum(){}
    ~sSum(){}
    sSum(const hypercube<data_t> &domain){
        successCheck(domain.getNdim()>=2,__FILE__,__LINE__,"The domain must contain at least 2 axes\n");
        _domain = domain;
        std::vector<ax > axes = domain.getAxes();
        axes.pop_back();
        _range = hyper(axes);
    }
    sSum * clone() const {
        sSum * op = new sSum(_domain);
        return op;
    }

    void apply_forward(bool add, const data_t * pmod, data_t * pdat){
        
        int n = _range.getN123();
        int ns = _domain.getN123()/n;

        const data_t (*pm) [n] = (const data_t (*)[n]) pmod;
        data_t sq_ns = sqrt(ns);

        #pragma omp parallel for
        for (int i=0; i<n; i++){
            data_t sum = 0;
            for (int s=0; s<ns; s++) sum += pm[s][i];
            pdat[i] = add*pdat[i] + sum/sq_ns;
        }
    }

    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

        if (!add) memset(pmod, 0, _domain.getN123()*sizeof(data_t));
        
        int n = _range.getN123();
        int ns = _domain.getN123()/n;

        data_t (*pm) [n] = (data_t (*)[n]) pmod;
        data_t sq_ns = sqrt(ns);

        #pragma omp parallel for
        for (int i=0; i<n; i++){
            for (int s=0; s<ns; s++) pm[s][i] += pdat[i]/sq_ns; 
        }
    }
};

// Regularization operator used in source-extended FWI
// Op = (I-S'.S).W.ms where ms is the extended model across sources, 'W' and 'S' are the operators defined above
class model_extension : public loper {
protected:
    sWeighting * _W;
    sSum * _S;
    std::shared_ptr<vec> _v1;
    std::shared_ptr<vec> _v2;

public:
    model_extension(){}
    ~model_extension(){delete _W; delete _S;}
    model_extension(sWeighting * W){
        _W = W->clone();
        _domain = *W->getDomain();
        _range = _domain;
        _S = new sSum(*W->getRange());
        _v1 = std::make_shared<vec>(_domain);
        _v2 = std::make_shared<vec>(*_S->getRange());
        _v1->zero();        
        _v2->zero();        
    }
    model_extension * clone() const {
        model_extension * op = new model_extension(_W);
        return op;
    }
    sWeighting * getW(){return _W;}
    sSum * getS(){return _S;}

    void apply_forward(bool add, const data_t * pmod, data_t * pdat){

        _W->apply_forward(false, pmod, _v1->getVals());
        _S->apply_forward(false, _v1->getCVals(), _v2->getVals());
        _v2->scale(-1);
        _S->apply_adjoint(add, pdat, _v2->getCVals());
        
        data_t * pv1 = _v1->getVals();

        #pragma omp parallel for
        for (int i=0; i<_domain.getN123(); i++) pdat[i] += pv1[i];
    }

    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

        _S->apply_forward(false, pdat, _v2->getVals());
        _S->apply_adjoint(false, _v1->getVals(), _v2->getCVals());

        data_t * pv1 = _v1->getVals();
        data_t * pv2 = _v2->getVals();

        #pragma omp parallel for
        for (int i=0; i<_domain.getN123(); i++) pv1[i] = pdat[i] - pv1[i];

        _W->apply_adjoint(add, pmod, pv1);
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
    std::vector<ax > axes = model->getHyper()->getAxes();
    int n=model->getN123();

// Analyze the inputs and parameters and modify if necessary
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0);
    analyzeBsplines(*model->getHyper(),par);
    analyzeNLInversion(par);
    par.sextension=true;

    std::shared_ptr<vec> gmask = nullptr;
    std::shared_ptr<vec> w = nullptr;
    std::shared_ptr<vec> filter = nullptr;
    std::shared_ptr<vec> invDiagH = nullptr;
    std::shared_ptr<vec> prior = nullptr;
    if (par.mask_file!="none") {gmask = read<data_t>(par.mask_file, par.format); successCheck(gmask->getN123()==n || gmask->getN123()==n*par.ns,__FILE__,__LINE__,"Gradient mask must have the same number of samples as the model\n");}
    if (par.weights_file!="none") {w = read<data_t>(par.weights_file, par.format); successCheck(w->getN123()==data->getN123(),__FILE__,__LINE__,"Data weights must have the same number of samples as the data\n");}
    if (par.filter_file!="none") {filter = read<data_t>(par.filter_file, par.format); successCheck(filter->getHyper()->getAxis(1).d==data->getHyper()->getAxis(1).d,__FILE__,__LINE__,"Filter and data must have the same sampling rate\n");}
    if (par.inverse_diagonal_hessian_file!="none") {invDiagH = read<data_t>(par.inverse_diagonal_hessian_file, par.format); successCheck(invDiagH->getN123()==n || invDiagH->getN123()==n*par.ns,__FILE__,__LINE__,"Inverse diagonal Hessian must have the same number of samples as the model\n");}
    if (par.prior_file!="none") {prior = read<data_t>(par.prior_file, par.format); successCheck(prior->getN123()==n || prior->getN123()==n*par.ns,__FILE__,__LINE__,"Prior model must have the same number of samples as the model\n");}


// ----------------------------------------------------------------------------------------//
// Extend model along sources, extend prior, mask and inverse Hessian if necessary
// ----------------------------------------------------------------------------------------//
    int nxz=n/par.nmodels;
    data_t vs0 = model->sum(nxz, 2*nxz);
    data_t rho0 = model->sum(2*nxz, 3*nxz);
    vs0/=nxz;
    rho0/=nxz;
{
    ax S(par.ns,0,1);
    axes.push_back(S);
    std::shared_ptr<vec> tmp = std::make_shared<vec>(vec(hyper(axes)));
    for (int s=0; s<par.ns; s++) memcpy(tmp->getVals()+s*n,model->getVals(),model->getN123()*sizeof(data_t));
    model=tmp;
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
        else successCheck(invDiagH->getN123()==n*par.ns,__FILE__,__LINE__,"The extended inverse Hessian has an incorrect size\n");
    }
    if (par.prior_file != "none"){
        if (prior->getN123() == n){
            std::shared_ptr<vec> tmp = std::make_shared<vec>(vec(hyper(axes)));
            for (int s=0; s<par.ns; s++) memcpy(tmp->getVals()+s*n,prior->getVals(),prior->getN123()*sizeof(data_t));
            prior=tmp;
        }
        else successCheck(prior->getN123()==n*par.ns,__FILE__,__LINE__,"The prior model has an incorrect size\n");
    }

    sWeighting * sW = new sWeighting(*model->getHyper(), par.sxz, dwidthx, dwidthz, dpower);
    model_extension * E = new model_extension(sW);
    delete sW;
// ----------------------------------------------------------------------------------------//

// Build model parameterization precon
// ----------------------------------------------------------------------------------------//
    std::shared_ptr<vec> model_temp = model;
    std::shared_ptr<vec> prior_temp = prior;

    nloper * P = nullptr;
    if (par.model_parameterization==0) P = new lam_mu_rho(*model->getHyper());
    else if (par.model_parameterization==2) P = new ip_is_rho(*model->getHyper());
    else if (par.model_parameterization==3) {
        P = new vs_vpvs_rho(*model->getHyper(), vs0, rho0);
        if (par.verbose>0) fprintf(stderr,"Average Vs and rho used in model parameterization are %.3f and %.3f\n",vs0/nxz, rho0/nxz);
    }
    if (P != nullptr){
        model_temp = std::make_shared<vec>(*model->getHyper());
        model_temp->zero();
        P->inverse(false,model_temp,model);
        if (par.prior_file != "none") {
            prior_temp = std::make_shared<vec>(*model->getHyper());
            prior_temp->zero();
            P->inverse(false,prior_temp,prior);
        }
    }
    model = model_temp;
    prior = prior_temp;
// ----------------------------------------------------------------------------------------//
  
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
    if (par.soft_clip) 
    {
        if (par.model_parameterization==0) {
            data_t mu_min = par.rhomin*par.vsmin*par.vsmin;
            data_t mu_max = par.rhomax*par.vsmax*par.vsmax;
            data_t lam_min = par.rhomin*(par.vpmin*par.vpmin - 2*par.vsmax*par.vsmax);
            data_t lam_max = par.rhomax*(par.vpmax*par.vpmax - 2*par.vsmin*par.vsmin);
            lam_min = std::max((data_t)0.0 , lam_min);
            emodelSoftClip S0(hyp0, lam_min, lam_max, mu_min, mu_max, par.rhomin, par.rhomax, 1, 9, 9);
            S = new emodelSoftClipExt(*model->getHyper(), &S0);
        }
        else if (par.model_parameterization==2){
            data_t ip_min = par.rhomin*par.vpmin;
            data_t ip_max = par.rhomax*par.vpmax;
            data_t is_min = par.rhomin*par.vsmin;
            data_t is_max = par.rhomax*par.vsmax;
            emodelSoftClip S0(hyp0, ip_min, ip_max, is_min, is_max, par.rhomin, par.rhomax, 1/sqrt(2.00001), 9, 9);
            S = new emodelSoftClipExt(*model->getHyper(), &S0);
        }
        else if (par.model_parameterization==3){
            emodelSoftClip S0(hyp0, log(par.vsmin/vs0), log(par.vsmax/vs0), -10.0, log(par.vpmax/par.vsmin - sqrt(2)), log(par.rhomin/rho0), log(par.rhomax/rho0), 1, 9, 9);
            S = new emodelSoftClipExt(*model->getHyper(), &S0);
        }
        else {
            emodelSoftClip S0(hyp0, par.vpmin, par.vpmax, par.vsmin, par.vsmax, par.rhomin, par.rhomax, 1/sqrt(2.00001), 9, 9);
            S = new emodelSoftClipExt(*model->getHyper(), &S0);
        }
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
        if (P != nullptr)
        {
            if (par.soft_clip) {
                chainNLOper SBD(S,BD);
                op  = new chainNLOper(P,&SBD);
            } 
            else op = new chainNLOper(P,BD);
        }
        else
        {
            if (par.soft_clip) op  = new chainNLOper(S,BD);
            else op = BD->clone();
        }
    }
    else
    {
        if (P != nullptr){
            if (par.soft_clip) op = new chainNLOper(P,S);
            else op = P->clone();
        }
        else{
            if (par.soft_clip) op = S->clone();
        }
    }
    if (par.bsplines) delete BD;
    if (P != nullptr) delete P;
    if (par.soft_clip) delete S;

    nloper * D = nullptr;
    loper * R;
    if (par.regularization==0) R = new identity (*model->getHyper());
    else if (par.regularization==1) R = new gradient2d(*model->getHyper(),par.reg_xweight,par.reg_zweight);
    else if (par.regularization==2) R = new laplacian2d(*model->getHyper(),par.reg_xweight,par.reg_zweight);
    else R = nullptr;

    if (R!=nullptr){
        if (op==nullptr) D = new chainNLOper(R,E);
        else {
            chainNLOper Eop(E,op);
            D = new chainNLOper(R, &Eop);
        }
        delete R;
    }
    else {
        if (op==nullptr) D = E->clone();
        else D = new chainNLOper(E, op);
    }

    nlls_fwi_reg * prob = new nlls_fwi_reg(L, D, bsmodel, data, par.lambda, bsprior, op, bsmask, w, filter);

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
    model = bsmodel;
    
    if (D != nullptr) delete D;
    if (op != nullptr)
    {
        std::shared_ptr<vec> tmp = std::make_shared<vec>(*op->getRange());
        tmp->zero();
        op->forward(false, model, tmp);
        model = tmp;
        delete op;
    }
    
// collapse the extended model into a single physical model mf = S.W.ms / sqrt(ns)
    std::shared_ptr<vec> fmodel = std::make_shared<vec>(*E->getS()->getRange());
    fmodel->zero();
    {
        std::shared_ptr<vec> temp = std::make_shared<vec>(*E->getW()->getRange());
        temp->zero();
        E->getW()->forward(false, model, temp);
        E->getS()->forward(false, temp, fmodel);
        fmodel->scale(1.0/sqrt(E->getW()->getNs()));
    }

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