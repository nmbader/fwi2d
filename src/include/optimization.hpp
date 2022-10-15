#pragma once

#define ZERO 1e-16
#define M_INF -1e+16

#include <time.h>
#include "operator.hpp"
#include "we_op.hpp"

#ifdef CUDA
    #include "cudaMisc.h"
#endif

// General class for optimization problems: min(f(m))
class optimization {
protected:
    std::shared_ptr<vecReg<data_t> > _d; // data vector
    std::shared_ptr<vecReg<data_t> > _m; // model vector
    std::shared_ptr<vecReg<data_t> > _r; // residual vector
    std::shared_ptr<vecReg<data_t> > _g; // functional gradient
public:
    optimization(){}
    virtual ~optimization(){}
    // default functional = 1/2 ||r||^2
    virtual data_t getFunc() {
        return 0.5*_r->norm2();
    }
    void initGrad() {
        _g = std::make_shared<vecReg<data_t> > (*_m->getHyper());
        _g->zero();
    }
    virtual void initRes() {
        _r = std::make_shared<vecReg<data_t> > (*_d->getHyper());
        _r->zero();
    }
    std::shared_ptr<vecReg<data_t> > getMod(){return _m;}
    std::shared_ptr<vecReg<data_t> > getGrad(){return _g;}
    std::shared_ptr<vecReg<data_t> > getRes(){return _r;}
    virtual void res() {}
    virtual void grad() {}
    virtual void hessian(const std::shared_ptr<vecReg<data_t> > H){};
    virtual data_t getZero(){return ZERO;}
};

// #################################################### Linear part ##################################################################
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// #################################################### ############### ############################################################## 

// Linear least-squares problem: f(m)=1/2.|Lm-d|^2
class llsq : public optimization{
protected:
    std::shared_ptr<vecReg<data_t> > _dr; // residual direction
    loper * _L; // linear operator

public:
    llsq(){}
    virtual ~llsq(){}
    llsq(loper * L, std::shared_ptr<vecReg<data_t> > m, std::shared_ptr<vecReg<data_t> > d){
       successCheck(L->checkDomainRange(m,d),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        _L = L;
        _m = m;
        _d = d;        
    }
    virtual void initDRes() {
        _dr = std::make_shared<vecReg<data_t> > (*_d->getHyper());
        _dr->zero();
    }
    std::shared_ptr<vecReg<data_t> > getDat(){return _d;}
    std::shared_ptr<vecReg<data_t> > getDRes(){return _dr;}

    // Lm - d
    virtual void res(){
        _L->forward(false,_m,_r);
        _r->scaleAdd(_d,1,-1);
    }
    // L'.(Lm-d)
    virtual void grad() {
        _L->adjoint(false,_g,_r);
    }
    virtual void dres(const std::shared_ptr<vecReg<data_t> > p){
        _L->forward(false,p,_dr);
    }
};

// Linear least-squares problem with regularization: f(m)=1/2.|Lm-d|^2 + 1/2.lambda^2.|D(m-m_prior)|^2
class llsq_reg : public llsq{
protected:
    loper * _D; // regularization linear operator
    data_t _lambda; // regularization damping parameter
    std::shared_ptr<vecReg<data_t> > _Dmp; // prior model pre-multiplied by D

public:
    llsq_reg(){}
    virtual ~llsq_reg(){}
    llsq_reg(loper * L, loper * D, std::shared_ptr<vecReg<data_t> > m, std::shared_ptr<vecReg<data_t> > d, const data_t lambda, const std::shared_ptr<vecReg<data_t> > mprior=nullptr){
       successCheck(L->checkDomainRange(m,d),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
       successCheck(D->checkDomainRange(m,m),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        successCheck(m->getN123()==mprior->getN123(),__FILE__,__LINE__,"Model and prior model do not have the same number of samples\n");
        _L = L;
        _D = D;
        _m = m;
        _d = d;
        _lambda = lambda;

        _Dmp = std::make_shared<vecReg<data_t> > (*_D->getRange());
        _Dmp->zero();
        if (mprior==nullptr) _D->forward(false, _m, _Dmp); // mprior assumed = m0 if not provided
        else _D->forward(false, mprior, _Dmp);
    }
    void initRes() {
        _r = std::make_shared<vecReg<data_t> >(hypercube<data_t>(_d->getN123()+_m->getN123()));
        _r->zero();
    }
    void initDRes() {
        _dr = std::make_shared<vecReg<data_t> >(hypercube<data_t>(_d->getN123()+_m->getN123()));
        _dr->zero();
    }
    // Lm - d ; lambda*D(m - m_prior)
    void res(){
        int nd = _d->getN123();
        int nm = _m->getN123();
        data_t * pr = _r->getVals();
        const data_t * pd = _d->getCVals();
        const data_t * pm = _m->getCVals();
        const data_t * pdmp = _Dmp->getCVals();
        
        _L->apply_forward(false,pm,pr);
        _D->apply_forward(false,pm,pr+nd);

        #pragma omp parallel for
        for (int i=0; i<nd; i++) pr[i] -= pd[i];
        #pragma omp parallel for
        for (int i=0; i<nm; i++) pr[nd+i] = _lambda*(pr[nd+i] - pdmp[i]);
    }
    // L'(Lm-d) + lambda.D'.D(m-m_prior)
    void grad() {
        _D->apply_adjoint(false,_g->getVals(),_r->getVals()+_d->getN123());
        _g->scale(_lambda);
        _L->apply_adjoint(true,_g->getVals(),_r->getVals());
    }
    void dres(const std::shared_ptr<vecReg<data_t> > p){
        int nd = _d->getN123();
        int nm = _m->getN123();
        data_t * pdr = _dr->getVals();
        _L->apply_forward(false,p->getVals(),_dr->getVals());
        _D->apply_forward(false,p->getVals(),_dr->getVals()+nd);
        _dr->scale(_lambda,nd,nd+nm);
    }

};

// Linear least-squares problem for SPD system: f(m)=1/2.m'.L.m -d'.m   L is SPD
// The minimizer is solution to Lm = d (m = L^-1.d)
class lspd : public llsq{
public:
    lspd(){}
    virtual ~lspd(){}
    lspd(loper * L, std::shared_ptr<vecReg<data_t> > m, std::shared_ptr<vecReg<data_t> > d){
       successCheck(L->checkDomainRange(m,d),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
       successCheck(m->getN123()==d->getN123(),__FILE__,__LINE__,"Model and data vectors must be of same size for SPD system\n");
        _L = L;
        _m = m;
        _d = d;
    }
    data_t getFunc(){
        data_t * pr = _r->getVals();
        data_t * pm = _m->getVals();
        data_t * pd = _d->getVals();
        data_t f = 0;
        for (int i=0; i<_m->getN123(); i++) f += 0.5*pm[i]*(pr[i]-pd[i]);
        return f;
    }
    // Lm - d
    void grad() {
        _g = _r;
    }
    data_t getZero(){return M_INF;}
};

// #################################################### Non-linear part ################################################################## //
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// //
// #################################################### ############### ################################################################## //

// Rosenbrock function: f(x,y)=(a-x)^2 + b(y-x^2)^2     m = (x,y)
// minimizer at (a,a^2) and fmin=0
// grad(f)=( 2(x-a) + 4bx(x^2-y) , 2b(y-x^2) )
class rosenbrock : public optimization{
protected:
    data_t _a, _b;
public:
    rosenbrock(){}
    virtual ~rosenbrock(){}
    rosenbrock(data_t a, data_t b, std::shared_ptr<vecReg<data_t> > m){
        successCheck(m->getN123()==2,__FILE__,__LINE__,"Rosenbrock model must contain two elements\n");
        _a = a;
        _b = b;
        _m = m;
        _d = std::make_shared<vecReg<data_t> >(hypercube<data_t>(1));
        _d->set(0);
    }
    // f(m)
    data_t getFunc(){
        data_t x = _m->getVals()[0];
        data_t y = _m->getVals()[1];
        return (_a-x)*(_a-x) + _b*(y-x*x)*(y-x*x);
    }
    // grad(f)
    void grad() {
        data_t x = _m->getVals()[0];
        data_t y = _m->getVals()[1];
        _g->getVals()[0] = 2*(x-_a) + 4*_b*x*(x*x-y);
        _g->getVals()[1] = 2*_b*(y-x*x); 
    }
    // Hessian(f)
    void hessian(std::shared_ptr<vecReg<data_t> > H){
        data_t x = _m->getVals()[0];
        data_t y = _m->getVals()[1];
        H->getVals()[0] = 2+4*_b*(x*x-y) + 8*_b*x*x;
        H->getVals()[1] = -4*_b*x;
        H->getVals()[2] = H->getVals()[1];
        H->getVals()[3] = 2*_b;
    }
};

// Paraboloid: f(x,y)=1/2(a.x)^2 + 1/2(b.y)^2     m = (x,y)
// minimizer at (0,0) and fmin=0
// grad(f)=(a^2.x,b^2.y)
class paraboloid : public optimization{
protected:
    data_t _a, _b;
public:
    paraboloid(){}
    virtual ~paraboloid(){}
    paraboloid(data_t a, data_t b, std::shared_ptr<vecReg<data_t> > m){
        successCheck(m->getN123()==2,__FILE__,__LINE__,"Paraboloid model must contain two elements\n");
        _a = a;
        _b = b;
        _m = m;
        _d = std::make_shared<vecReg<data_t> >(hypercube<data_t>(1));
        _d->set(0);
    }
    // f(m)
    data_t getFunc(){
        data_t x = _m->getVals()[0];
        data_t y = _m->getVals()[1];
        return 0.5*((_a*_a*x*x) + (_b*_b*y*y));
    }
    // grad(f)
    void grad() {
        data_t x = _m->getVals()[0];
        data_t y = _m->getVals()[1];
        _g->getVals()[0] = _a*_a*x;
        _g->getVals()[1] = _b*_b*y; 
    }
    // Hessian(f)
    void hessian(std::shared_ptr<vecReg<data_t> > H){
        H->getVals()[0] = _a*_a;
        H->getVals()[1] = 0;
        H->getVals()[2] = 0;
        H->getVals()[3] = _b*_b;
    }
};

// Non-linear least-squares problem: f(m)=1/2.|L(m)-d|^2
class nlls : public optimization{
protected:
    nloper * _L; // non-linear operator

public:
    nlls(){}
    virtual ~nlls(){}
    nlls(nloper * L, std::shared_ptr<vecReg<data_t> > m, std::shared_ptr<vecReg<data_t> > d){
       successCheck(L->checkDomainRange(m,d),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        _L = L;
        _m = m;
        _d = d;
    }
    // L(m) - d
    virtual void res(){
        _L->forward(false,_m,_r);
        _r->scaleAdd(_d,1,-1);
    }
    // (dL(m)/dm)'.(Lm-d)
    virtual void grad() {
        _L->jacobianT(false,_g,_m,_r);
    }
};

// Non-linear least-squares problem with time quadrature, data weighting and gradient mask for FWI
// f(m)=1/2.(L(m)-d)'.Ht.(L(m)-d) = 1/2.r'.Ht.r
// optionally other functional may be considered:
// 1) with weighting: r <- W.r
// 2) with trace normalization: r <- (u/|u| - d/|d|) where u = L(m)
// 3) envelop or envelop squared: r <- E(u) - E(d)
class nlls_fwi : public optimization{
protected:
    nloper * _L; // non-linear WE operator
    std::shared_ptr<vecReg<data_t> > _gmask; // gradient mask vector
    std::shared_ptr<vecReg<data_t> > _w; // residual weighting vector
    bool _normalize; // apply or not trace normalization
    bool _integrate; // apply or not time integration
    std::shared_ptr<vecReg<data_t> > _norms; // traces norm needed for the Jacobian of the normalization
    std::shared_ptr<vecReg<data_t> > _syn1; // synthetic data needed for the Jacobian of the normalization
    int _envelop; // compute or not the envelop (or envelop squared) of the data
    std::shared_ptr<vecReg<data_t> > _syn2; // synthetic data needed for the Jacobian of the envelop

public:
    nlls_fwi(){}
    virtual ~nlls_fwi(){}
    nlls_fwi(nloper * L, std::shared_ptr<vecReg<data_t> > m, std::shared_ptr<vecReg<data_t> > d, std::shared_ptr<vecReg<data_t> > gmask = nullptr, std::shared_ptr<vecReg<data_t> > w = nullptr, bool normalize=false, bool integrate=false, int envelop=0){
       successCheck(L->checkDomainRange(m,d),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        if (gmask != nullptr) {
            successCheck(m->getN123()==gmask->getN123(),__FILE__,__LINE__,"The gradient mask and model vectors must have the same size\n");
            successCheck(gmask->min()>=0,__FILE__,__LINE__,"Gradient mask must be non-negative\n");
        }
        if (w != nullptr) {
            successCheck(d->getN123()==w->getN123(),__FILE__,__LINE__,"The weighting and data vectors must have the same size\n");
            successCheck(w->min()>=0,__FILE__,__LINE__,"Data weights  must be non-negative\n");
        }
        _L = L;
        _m = m;
        _d = d;
        _gmask = gmask;
        _w = w;
        _normalize=normalize;
        _integrate=integrate;
        if (envelop==1 || envelop==2) _envelop = envelop;
        else _envelop=0;

        if (_integrate){
            integral S(*_d->getHyper());
            S.forward(false, _d, _d);
        }
        if (_w != nullptr) _d->mult(_w);
        if (_normalize) {
            int nt=_d->getHyper()->getAxis(1).n;
            int ntr=_d->getN123()/nt;
            data_t * norms = new data_t[ntr];
            ttnormalize(_d->getVals(), norms, nt, ntr);
            delete [] norms;
            _norms = std::make_shared<vecReg<data_t> >(hypercube<data_t>(ntr));
            _syn1 = std::make_shared<vecReg<data_t> >(hypercube<data_t>(_d->getN123()));
        }
        if (_envelop==1) envelop1(_d);
        else if (_envelop==2) envelop2(_d);
        if (_envelop!=0) _syn2 = std::make_shared<vecReg<data_t> >(hypercube<data_t>(_d->getN123()));
    }

    // r = (L(m) - d) = u - d or W.(u-d) or (u/|u|-d/|d|) or (E(u)-E(d)) 
    void res(){
        _L->forward(false,_m,_r);

        if (_integrate){
            integral S(*_r->getHyper());
            S.forward(false, _r, _r);
        }

        if (_w != nullptr) _r->mult(_w);

        if (_normalize){
            int ntr = _norms->getN123();
            int nt = _r->getN123()/ntr;
            memcpy(_syn1->getVals(), _r->getVals(), _r->getN123()*sizeof(data_t));
            ttnormalize(_r->getVals(), _norms->getVals(), nt, ntr);
        }

        if (_envelop != 0){
            memcpy(_syn2->getVals(), _r->getVals(), _r->getN123()*sizeof(data_t));
            if (_envelop==1) envelop1(_r);
            else envelop2(_r);
        }

        _r->scaleAdd(_d,1,-1);
    }
    data_t getFunc(){
        int n123 = _r->getN123();
        int nt = _r->getHyper()->getAxis(1).n;
        int nx = n123 / nt;
        data_t dt = _r->getHyper()->getAxis(1).d;
        data_t * pr = _r->getVals();
        data_t f = _r->norm2();
        for (int i=0; i<nx; i++) f = f - 0.5*(pr[i*nt]*pr[i*nt] + pr[i*nt+nt-1]*pr[i*nt+nt-1]);
        return 0.5*dt*f;
    }
    // (dL(m)/dm)'.Ht.r or (dL(m)/dm)'.W'.Ht.r
    // or (dL(m)/dm)'.( I/||u|| - uu'/||u||^3 ).Ht.r  u = syn1
    // or  (dL(m)/dm)'.[Diag(u/E(u)) - H.Diag((H.u)/E(u))].Ht.r for envelop or (dL(m)/dm)'.2[Diag(u) - H.Diag(H.u)].Ht.r for squared envelop
    // H. is the Hilbert transform operator ; E is the envelop ; u = syn2
    void grad() {
        int n123 = _r->getN123();
        int nt = _r->getHyper()->getAxis(1).n;
        int nx = n123 / nt;
        data_t dt = _r->getHyper()->getAxis(1).d;
        std::shared_ptr<vecReg<data_t> > r = _r->clone();
        data_t * pr = r->getVals();
        applyHt(false, false, pr, pr, nx, nt, dt, 0, nx);

        if (_envelop != 0){
            std::shared_ptr<vecReg<data_t> > temp = _syn2->clone();
            hilbert(temp);
            std::shared_ptr<vecReg<data_t> > temp2;
            data_t * ptemp = temp->getVals();
            data_t * psyn = _syn2->getVals();
            if (_envelop==1) {
                temp2 = temp->clone();
                data_t * ptemp2 = temp2->getVals();
                for (int i=0; i<temp->getN123(); i++) {
                    ptemp2[i] = std::max((data_t)ZERO,(data_t)sqrt(ptemp[i]*ptemp[i]+psyn[i]*psyn[i]));
                    ptemp[i] *= pr[i]/ptemp2[i];
                }
            }
            else{
                for (int i=0; i<temp->getN123(); i++) ptemp[i] *= pr[i];
            }
            hilbert(temp);
            ptemp = temp->getVals();
            if (_envelop==1) {
                data_t * ptemp2 = temp2->getVals();
                for (int i=0; i<temp->getN123(); i++) pr[i] = psyn[i]/ptemp2[i] * pr[i] - ptemp[i];
            }
            else {
                for (int i=0; i<temp->getN123(); i++) pr[i] = 2*(psyn[i]*pr[i]-ptemp[i]);
            }
        }

        if (_normalize){
            int ntr = _norms->getN123();
            int nt = _r->getN123()/ntr;
            data_t * psyn = _syn1->getVals();
            data_t * pnorm = _norms->getVals();
            for (int ix=0; ix<ntr; ix++){
                data_t val = 0;
                int i;
                for (int it=0; it<nt; it++){
                    i = ix*nt+it;
                    val += pr[i] * psyn[i]; // u'.r
                }
                for (int it=0; it<nt; it++){
                    i = ix*nt+it;
                    pr[i] = pr[i]/pnorm[ix] - val/(pnorm[ix]*pnorm[ix]*pnorm[ix]) * psyn[i];
                }
            }
        }

        if (_w != nullptr) r->mult(_w);

        if (_integrate){
            integral S(*r->getHyper());
            S.adjoint(false, r, r);
        }

        _L->jacobianT(false,_g,_m,r);

        if (_gmask != nullptr) _g->mult(_gmask);
    }
};


// Same class as nlls_fwi but re-written so that the gradient is computed simultaneously with the residual for better scalability over the number of shots
// the functional is normalized: f(m)=1/2.[(L(m)-d)'.Ht.(L(m)-d)]/(d'.Ht.d) 
// The WE operator and model preconditioner are provided separately
class nlls_fwi_eco : public optimization{
protected:
    nl_we_op * _L; // non-linear WE operator
    nloper * _P; // model preconditioner
    std::shared_ptr<vecReg<data_t> > _p; // preconditioned model
    std::shared_ptr<vecReg<data_t> > _pg; // gradient before applying preconditioner Jacobian
    std::shared_ptr<vecReg<data_t> > _gmask; // gradient mask vector
    std::shared_ptr<vecReg<data_t> > _w; // residual weighting vector
    std::shared_ptr<vecReg<data_t> > _filter; // 1D filter
    data_t _dnorm; // data normalization factor
    data_t _f; // objective function
    bool _flag; // flag to re-evaluate or not the objective function
    int _scale_source_times; // for source inversion by reduced Variable Projection

public:
    std::vector <data_t> _dfunc; // vector that stores data objective function values for all trials (used in the regularization class)
    std::vector <data_t> _mfunc; // vector that stores model objective function values for all trials (used in the regularization class)
    std::vector<data_t> _scalers; // scalers applied to each source at a given trial
    nlls_fwi_eco(){}
    virtual ~nlls_fwi_eco(){}
    nlls_fwi_eco(nl_we_op * L, std::shared_ptr<vecReg<data_t> > m, std::shared_ptr<vecReg<data_t> > d, nloper * P = nullptr, std::shared_ptr<vecReg<data_t> > gmask = nullptr, std::shared_ptr<vecReg<data_t> > w = nullptr, std::shared_ptr<vecReg<data_t> > filter = nullptr){
        _L = L;
        _P = P;
        _m = m;
        _d = d;
        _gmask = gmask;
        _w = w;
        _filter = filter;
        _f = 0;
        _flag = true;
        _d->setHyper(*L->getRange());
       
        if (P!= nullptr){
            _p = std::make_shared<vecReg<data_t> >(*P->getRange());
            _pg = std::make_shared<vecReg<data_t> >(*P->getRange());
            _p->zero();
            _pg->zero();
            successCheck(P->checkDomainRange(m,_p),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        }
        else{
            _p = _m;
            _pg = _g;
        }

        successCheck(L->checkDomainRange(_p,_d),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        if (gmask != nullptr) {
            successCheck(m->getN123()==gmask->getN123(),__FILE__,__LINE__,"The gradient mask and model vectors must have the same size\n");
            successCheck(gmask->min()>=0,__FILE__,__LINE__,"Gradient mask must be non-negative\n");
        }
        if (w != nullptr) {
            successCheck(d->getN123()==w->getN123(),__FILE__,__LINE__,"The weighting and data vectors must have the same size\n");
            successCheck(w->min()>=0,__FILE__,__LINE__,"Data weights  must be non-negative\n");
        }
        
        if (_L->_par.integrate){ // apply time integration
            integral S(*_d->getHyper());
            S.forward(false, _d, _d);
        }
        if (_L->_par.double_difference){ // apply double difference between consecutive traces
            xdifference xdiff(*_d->getHyper());
            xdiff.forward(false, _d, _d);
        }
        if (_w != nullptr) _d->mult(_w); // apply weighting

        if (_filter != nullptr){ // apply filtering
            data_t eps=1e-07;
            if (_L->_par.filter_phase=="zero") _filter = zero_phase(filter);
            else if (_L->_par.filter_phase=="minimum") _filter = minimum_phase(filter,eps);
            axis<data_t> Tf = _filter->getHyper()->getAxis(1);
            axis<data_t> T = _d->getHyper()->getAxis(1);
            Tf.d=T.d;
            _filter->setHyper(hypercube<data_t>(Tf)); 
            conv1dnd op(*_d->getHyper(), _filter, _L->_par.filter_phase!="minimum");
            std::shared_ptr<vecReg<data_t> > output = std::make_shared<vecReg<data_t> >(*op.getRange());
            output->zero();
            op.forward(false,_d,output);
            _d=output;
        }

        if (_L->_par.normalize) { // apply trace normalization
            int nt=_d->getHyper()->getAxis(1).n;
            int ntr=_d->getN123()/nt;
            data_t * norms = new data_t[ntr];
            ttnormalize(_d->getVals(), norms, nt, ntr);
            delete [] norms;
        }
        if (_L->_par.envelop==1) envelop1(_d); // compute the envelop (or envelop squared) of the data
        else if (_L->_par.envelop==2) envelop2(_d);

        if (_L->_par.interferometry){ // trace to trace deconvolution, for interferometric FWI
            // compute the maxium of the amplitude spectrum of the data
            fxTransform fx(*_d->getHyper());
            std::shared_ptr<cvecReg<data_t> > fxvec = std::make_shared<cvecReg<data_t> > (*fx.getRange());
            fxvec->zero();
            fx.forward(false,_d,fxvec);
            std::shared_ptr<vecReg<data_t> > spec=fxvec->modulus();
            data_t max_spec=spec->max();
            _L->_par.epsilon *= max_spec;

            tr2trDecon decon(*_d->getHyper(),_L->_par.wmin,_L->_par.wmax, _L->_par.epsilon, _L->_par.smth_half_length);
            std::shared_ptr<vecReg<data_t> > output = std::make_shared<vecReg<data_t> >(*decon.getRange());
            output->zero();
            decon.forward(false,_d,output);
            _d = output;
        }

        if (_L->_par.normalize_obj_func){
            _dnorm = _d->norm2();
            _dnorm *= _d->getHyper()->getAxis(1).d;
            if (_dnorm<ZERO) _dnorm=1;
        }
        else _dnorm=1.0;

        _scale_source_times=0;
        if (_L->_par.scale_source_times>0) {
            for (int s=0; s<_L->_par.ns; s++) _scalers.push_back(1.0);
        }
    }

    void compute_res_and_grad(data_t * r){       

        _g->zero();

        if (_P != nullptr)  {
            _P->forward(false, _m, _p);
            _pg->zero();
        }
        else _pg=_g;

        hypercube<data_t> hyp(_p->getHyper()->getAxis(1),_p->getHyper()->getAxis(2),_p->getHyper()->getAxis(3));
        int ncxz = hyp.getN123();

        int ns = _L->_par.ns;
        int nt = _d->getHyper()->getAxis(1).n;
        data_t dt = _d->getHyper()->getAxis(1).d;

        int size=1, rank=0;
#ifdef ENABLE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
        time_t t = time(NULL);
        if (_L->_par.verbose>1 && rank==0) fprintf(stderr,"\n====================\n%s\n====================\n",ctime(&t));

        for (int s=rank; s<ns; s+=size)
        {
            if (_L->_par.verbose>1) fprintf(stderr,"Start processing shot %d by process %d\n",s, rank);
            
            // cumulative number of receivers
            int nr=0;
            if (s>0) {for (int i=0; i<s; i++) nr += _L->_par.rxz[i].size();}

            // adjust param object to a single shot
            param par = _L->_par;
            par.sxz = {_L->_par.sxz[s]};
            par.rxz = {_L->_par.rxz[s]};
            par.ns=1;
            par.nr=par.rxz[s].size();
            par.skip_mpi=true;
            //par.verbose=0;
            if (par.gl>0) par.rdip = std::vector<data_t>(_L->_par.rdip.begin()+nr, _L->_par.rdip.begin()+nr+par.nr);

            if (s>0) par.verbose=0;

            // extract the time functions for a single shot
            hypercube<data_t> hyper_s = *_L->_allsrc->getHyper();
            std::shared_ptr<vecReg<data_t> > src = std::make_shared<vecReg<data_t> > (hypercube<data_t>(hyper_s.getAxis(1),axis<data_t>(1,0,1),hyper_s.getAxis(3)));
            memcpy(src->getVals(), _L->_allsrc->getVals() + s*nt, sizeof(data_t)*nt);
            if ((par.nmodels>=3 && !par.acoustic_elastic) || (par.acoustic_elastic && !par.acoustic_source) ){
                memcpy(src->getVals()+nt, _L->_allsrc->getVals() + (ns+s)*nt, sizeof(data_t)*nt);
                if (par.mt) memcpy(src->getVals()+2*nt, _L->_allsrc->getVals() + (2*ns+s)*nt, sizeof(data_t)*nt);
            }

            // build the we operator for a single shot
            nl_we_op * L;
            if (par.nmodels==2) L=new nl_we_op_a(hyp,src,par);
            else if (par.nmodels==3 && !par.acoustic_elastic) L=new nl_we_op_e(hyp,src,par);
            else if (par.nmodels==3 && par.acoustic_elastic) L=new nl_we_op_ae(hyp,src,par);
            else if (par.nmodels==5) L=new nl_we_op_vti(hyp,src,par);
            std::shared_ptr<vecReg<data_t> > rs = std::make_shared<vecReg<data_t> >(*L->getRange());
            int ntr = rs->getN123()/nt;
            
            L->apply_forward(false,_p->getVals()+s*par.sextension*ncxz,rs->getVals());
            data_t * pr = rs->getVals();

            if (par.integrate){
                integral S(*rs->getHyper());
                S.forward(false, rs, rs);
            }

            if (par.double_difference){
                xdifference xdiff(*rs->getHyper());
                xdiff.forward(false, rs, rs);
            }

            if (_w != nullptr) {
                data_t * pw = _w->getVals()+nr*nt;
                // first component
                #pragma omp parallel for
                for (int i=0; i<par.nr*nt; i++) pr[i] *= pw[i];

                // second component
                if ((par.nmodels>=3 && par.gl==0) || par.acoustic_elastic){
                    #pragma omp parallel for
                    for (int i=0; i<par.nr*nt; i++) pr[par.nr*nt+i] *= pw[_L->_par.nr*nt+i];
                }

                // third component
                if (par.gl==0 && par.acoustic_elastic){
                    #pragma omp parallel for
                    for (int i=0; i<par.nr*nt; i++) pr[2*par.nr*nt+i] *= pw[2*_L->_par.nr*nt+i];
                }
            }

            if (_filter != nullptr){
                data_t eps=1e-07;
                conv1dnd op(*rs->getHyper(), _filter, _L->_par.filter_phase!="minimum");
                std::shared_ptr<vecReg<data_t> > output = std::make_shared<vecReg<data_t> >(*op.getRange());
                output->zero();
                op.forward(false,rs,output);
                rs=output;
                pr = rs->getVals();
            }

            // rescale the synthetics by u'd/u'u (equivalent to Variable Projection with the variable being a single scaler multiplying the source time function)
            data_t * pd = _d->getVals()+nr*nt;
            if (par.scale_source_times>0){
                if (par.scale_source_times>_scale_source_times){
                    data_t scaler1=0;
                    data_t scaler2=0;
                    // first component
                    #pragma omp parallel for reduction(+: scaler1,scaler2)
                    for (int i=0; i<par.nr*nt; i++) {
                        scaler1 += pr[i]*pd[i];
                        scaler2 += pr[i]*pr[i];
                    }
                    // second component
                    if ((par.nmodels>=3 && par.gl==0) || par.acoustic_elastic){
                        #pragma omp parallel for reduction(+: scaler1,scaler2)
                        for (int i=0; i<par.nr*nt; i++) {
                            scaler1 += pr[par.nr*nt+i]*pd[_L->_par.nr*nt+i];
                            scaler2 += pr[par.nr*nt+i]*pr[par.nr*nt+i];
                        }
                    }
                    // third component
                    if (par.gl==0 && par.acoustic_elastic){
                        #pragma omp parallel for reduction(+: scaler1,scaler2)
                        for (int i=0; i<par.nr*nt; i++) {
                            scaler1 += pr[2*par.nr*nt+i]*pd[2*_L->_par.nr*nt+i];
                            scaler2 += pr[2*par.nr*nt+i]*pr[2*par.nr*nt+i];
                        }
                    }
                    _scalers[s] = scaler1/scaler2;
                    if (std::abs(std::log(_scalers[s]))>par.scale_source_log_clip) _scalers[s]=std::exp(((_scalers[s]>=1) - (_scalers[s]<1))*par.scale_source_log_clip);
                }
                rs->scale(_scalers[s]);
                if (_L->_par.verbose>0) fprintf(stderr,"Shot %d rescaled by a factor of %f\n",s,_scalers[s]);
            }

            std::shared_ptr<vecReg<data_t> > norms;
            std::shared_ptr<vecReg<data_t> > syn1;
            if (par.normalize){
                norms = std::make_shared<vecReg<data_t> >(hypercube<data_t>(ntr));
                syn1 = rs->clone();
                ttnormalize(rs->getVals(), norms->getVals(), nt, ntr);
            }

            std::shared_ptr<vecReg<data_t> > syn2;
            if (par.envelop != 0){
                syn2 = rs->clone();
                if (par.envelop==1) envelop1(rs);
                else envelop2(rs);
            }

            std::shared_ptr<vecReg<data_t> > syn3;
            if (_L->_par.interferometry){
                syn3 = rs->clone();
                tr2trDecon decon(*rs->getHyper(),_L->_par.wmin,_L->_par.wmax, _L->_par.epsilon, _L->_par.smth_half_length);
                std::shared_ptr<vecReg<data_t> > output = std::make_shared<vecReg<data_t> >(*decon.getRange());
                output->zero();
                decon.forward(false,rs,output);
                rs = output;
            }

            pr = rs->getVals();

            // compute the residual u-d
            // first component
            #pragma omp parallel for
            for (int i=0; i<par.nr*nt; i++) pr[i] -= pd[i];

            // second component
            if ((par.nmodels>=3 && par.gl==0) || par.acoustic_elastic){
                #pragma omp parallel for
                for (int i=0; i<par.nr*nt; i++) pr[par.nr*nt+i] -= pd[_L->_par.nr*nt+i];
            }
            // third component
            if (par.gl==0 && par.acoustic_elastic){
                #pragma omp parallel for
                for (int i=0; i<par.nr*nt; i++) pr[2*par.nr*nt+i] -= pd[2*_L->_par.nr*nt+i];
            }

            memcpy(r+nr*nt, rs->getVals(), par.nr*nt*sizeof(data_t));
            if ((par.nmodels>=3 && par.gl==0) || par.acoustic_elastic) memcpy(r+(_L->_par.nr+nr)*nt, rs->getVals()+par.nr*nt, par.nr*nt*sizeof(data_t));
            if (par.gl==0 && par.acoustic_elastic) memcpy(r+(2*_L->_par.nr+nr)*nt, rs->getVals()+2*par.nr*nt, par.nr*nt*sizeof(data_t));

            // compute the gradient per shot
            applyHt(false, false, pr, pr, ntr, nt, dt, 0, ntr);

            if (_L->_par.interferometry){
                tr2trDecon decon(*rs->getHyper(),_L->_par.wmin,_L->_par.wmax, _L->_par.epsilon, _L->_par.smth_half_length);
                std::shared_ptr<vecReg<data_t> > output = std::make_shared<vecReg<data_t> >(*decon.getDomain());
                output->zero();
                decon.jacobianT(false,output,syn3,rs);
                rs=output;
                pr = rs->getVals();
            }

            if (par.envelop != 0){
                std::shared_ptr<vecReg<data_t> > temp = syn2->clone();
                hilbert(temp);
                std::shared_ptr<vecReg<data_t> > temp2;
                data_t * ptemp = temp->getVals();
                data_t * psyn = syn2->getVals();
                if (par.envelop==1) {
                    temp2 = temp->clone();
                    data_t * ptemp2 = temp2->getVals();
                    for (int i=0; i<temp->getN123(); i++) {
                        ptemp2[i] = std::max((data_t)ZERO,(data_t)sqrt(ptemp[i]*ptemp[i]+psyn[i]*psyn[i]));
                        ptemp[i] *= pr[i]/ptemp2[i];
                    }
                }
                else{
                    for (int i=0; i<temp->getN123(); i++) ptemp[i] *= pr[i];
                }
                hilbert(temp);
                ptemp = temp->getVals();
                if (par.envelop==1) {
                    data_t * ptemp2 = temp2->getVals();
                    for (int i=0; i<temp->getN123(); i++) pr[i] = psyn[i]/ptemp2[i] * pr[i] - ptemp[i];
                }
                else {
                    for (int i=0; i<temp->getN123(); i++) pr[i] = 2*(psyn[i]*pr[i]-ptemp[i]);
                }
            }

            if (par.normalize){
                data_t * psyn = syn1->getVals();
                data_t * pnorm = norms->getVals();
                for (int ix=0; ix<ntr; ix++){
                    data_t val = 0;
                    int i;
                    for (int it=0; it<nt; it++){
                        i = ix*nt+it;
                        val += pr[i] * psyn[i]; // u'.r
                    }
                    for (int it=0; it<nt; it++){
                        i = ix*nt+it;
                        pr[i] = pr[i]/pnorm[ix] - val/(pnorm[ix]*pnorm[ix]*pnorm[ix]) * psyn[i];
                    }
                }
            }

            if (_filter != nullptr){
                data_t eps=1e-07;
                conv1dnd op(*rs->getHyper(), _filter, _L->_par.filter_phase!="minimum");
                std::shared_ptr<vecReg<data_t> > output = std::make_shared<vecReg<data_t> >(*op.getDomain());
                output->zero();
                op.adjoint(false,output,rs);
                rs=output;
                pr = rs->getVals();
            }

            if (_w != nullptr) {
                data_t * pw = _w->getVals()+nr*nt;
                // first component
                #pragma omp parallel for
                for (int i=0; i<par.nr*nt; i++) pr[i] *= pw[i];

                // second component
                if ((par.nmodels>=3 && par.gl==0) || par.acoustic_elastic){
                    #pragma omp parallel for
                    for (int i=0; i<par.nr*nt; i++) pr[par.nr*nt+i] *= pw[_L->_par.nr*nt+i];
                }

                // third component
                if (par.gl==0 && par.acoustic_elastic){
                    #pragma omp parallel for
                    for (int i=0; i<par.nr*nt; i++) pr[2*par.nr*nt+i] *= pw[2*_L->_par.nr*nt+i];
                }
            }

            if (par.double_difference){
                xdifference xdiff(*rs->getHyper());
                xdiff.adjoint(false, rs, rs);
            }

            if (par.integrate){
                integral S(*rs->getHyper());
                S.adjoint(false, rs, rs);
            }

            rs->scale(1.0/_dnorm);

            L->apply_jacobianT(true,_pg->getVals()+s*par.sextension*ncxz,_p->getVals()+s*par.sextension*ncxz,rs->getVals());

            delete L;

            if (_L->_par.verbose>1) fprintf(stderr,"Finish processing shot %d by process %d\n",s, rank);

        } // end of loop over shots

#ifdef ENABLE_MPI
        // Sum all gradients
        data_t * gtemp = new data_t[_pg->getN123()];
        memset(gtemp, 0, _pg->getN123()*sizeof(data_t));
        if (sizeof(data_t)==8) MPI_Allreduce(_pg->getVals(), gtemp, _pg->getN123(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        else MPI_Allreduce(_pg->getVals(), gtemp, _pg->getN123(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        memcpy(_pg->getVals(), gtemp, _pg->getN123()*sizeof(data_t));
        delete [] gtemp;
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        
        if (_L->_par.scale_source_times>0) _scale_source_times++;
    }

    virtual void res(){
        _r->zero();
        compute_res_and_grad(_r->getVals()); _flag = true;
#ifdef ENABLE_MPI
        // gather all residual
        data_t * temp = new data_t[_r->getN123()];
        memset(temp, 0, _r->getN123()*sizeof(data_t));
        if (sizeof(data_t)==8) MPI_Allreduce(_r->getVals(), temp, _r->getN123(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        else MPI_Allreduce(_r->getVals(), temp, _r->getN123(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        memcpy(_r->getVals(), temp, _r->getN123()*sizeof(data_t));
        delete [] temp;
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        }

    virtual data_t getFunc(){
        if (_flag)
        {
            int n123 = _r->getN123();
            int nt = _r->getHyper()->getAxis(1).n;
            int nx = n123 / nt;
            data_t dt = _r->getHyper()->getAxis(1).d;
            data_t * pr = _r->getVals();
            data_t f = _r->norm2();
            for (int i=0; i<nx; i++) f = f - 0.5*(pr[i*nt]*pr[i*nt] + pr[i*nt+nt-1]*pr[i*nt+nt-1]);
            _f = 0.5*dt*f/_dnorm;
            _flag = false;
        }
        return _f;
    }
    virtual void grad() {
        // already computed in the compute_res_and_grad() method above
        if (_P != nullptr) _P->jacobianT(false,_g,_m,_pg);
        if (_gmask != nullptr) _g->mult(_gmask);
    }
};

// Same class as nlls_fwi_eco but with model regularizations
// f(m)=1/2.[(L(m)-d)'.Ht.(L(m)-d)]/(d'.Ht.d)  + 1/2.lambda^2.|D(m)-D(m_prior)|^2 / |D(m_prior)|^2
class nlls_fwi_reg : public nlls_fwi_eco{
protected:
    nloper * _D; // regularization operator
    data_t _lambda; // regularization damping parameter
    std::shared_ptr<vecReg<data_t> > _Dmp; // prior model pre-multiplied by D
    std::shared_ptr<vecReg<data_t> > _dg; // gradient component corresponding to the regularization D
    data_t _mnorm; // model normalization factor

public:
    nlls_fwi_reg(){}
    virtual ~nlls_fwi_reg(){}
    nlls_fwi_reg(nl_we_op * L, nloper * D, std::shared_ptr<vecReg<data_t> > m, std::shared_ptr<vecReg<data_t> > d, data_t lambda, std::shared_ptr<vecReg<data_t> > mprior = nullptr, nloper * P = nullptr, std::shared_ptr<vecReg<data_t> > gmask = nullptr, std::shared_ptr<vecReg<data_t> > w = nullptr, std::shared_ptr<vecReg<data_t> > filter = nullptr){
        _L = L;
        _D = D;
        _P = P;
        _m = m;
        _d = d;
        _lambda = lambda;    
        _gmask = gmask;
        _w = w;
        _filter = filter;
        _f = 0;
        _flag = true;
        _d->setHyper(*L->getRange());

        _Dmp = std::make_shared<vecReg<data_t> > (*_D->getRange());
        _Dmp->zero();
        if (mprior!=nullptr) _D->forward(false, mprior, _Dmp); // mprior assumed = 0 if not provided

        _dg = std::make_shared<vecReg<data_t> >(*m->getHyper());
        _dg->zero();
       
        if (P!= nullptr){
            _p = std::make_shared<vecReg<data_t> >(*P->getRange());
            _pg = std::make_shared<vecReg<data_t> >(*P->getRange());
            _p->zero();
            _pg->zero();
            successCheck(P->checkDomainRange(m,_p),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        }
        else{
            _p = _m;
            _pg = _g;
        }
        successCheck(L->checkDomainRange(_p,_d),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        if (gmask != nullptr) {
            successCheck(m->getN123()==gmask->getN123(),__FILE__,__LINE__,"The gradient mask and model vectors must have the same size\n");
            successCheck(gmask->min()>=0,__FILE__,__LINE__,"Gradient mask must be non-negative\n");
        }
        if (w != nullptr) {
            successCheck(d->getN123()==w->getN123(),__FILE__,__LINE__,"The weighting and data vectors must have the same size\n");
            successCheck(w->min()>=0,__FILE__,__LINE__,"Data weights  must be non-negative\n");
        }
        
        if (_L->_par.integrate){
            integral S(*_d->getHyper());
            S.forward(false, _d, _d);
        }
        if (_L->_par.double_difference){ // apply double difference between consecutive traces
            xdifference xdiff(*_d->getHyper());
            xdiff.forward(false, _d, _d);
        }
        if (_w != nullptr) _d->mult(_w);

        if (_filter != nullptr){ // apply filtering
            data_t eps=1e-07;
            if (_L->_par.filter_phase=="zero") _filter = zero_phase(filter);
            else if (_L->_par.filter_phase=="minimum") _filter = minimum_phase(filter,eps);
            axis<data_t> Tf = _filter->getHyper()->getAxis(1);
            axis<data_t> T = _d->getHyper()->getAxis(1);
            Tf.d=T.d;
            _filter->setHyper(hypercube<data_t>(Tf)); 
            conv1dnd op(*_d->getHyper(), _filter, _L->_par.filter_phase!="minimum");
            std::shared_ptr<vecReg<data_t> > output = std::make_shared<vecReg<data_t> >(*op.getRange());
            output->zero();
            op.forward(false,_d,output);
            _d=output;
        }

        if (_L->_par.normalize) {
            int nt=_d->getHyper()->getAxis(1).n;
            int ntr=_d->getN123()/nt;
            data_t * norms = new data_t[ntr];
            ttnormalize(_d->getVals(), norms, nt, ntr);
            delete [] norms;
        }
        if (_L->_par.envelop==1) envelop1(_d);
        else if (_L->_par.envelop==2) envelop2(_d);

        if (_L->_par.interferometry){ // trace to trace deconvolution, for interferometric FWI
            // compute the maxium of the amplitude spectrum of the data
            fxTransform fx(*_d->getHyper());
            std::shared_ptr<cvecReg<data_t> > fxvec = std::make_shared<cvecReg<data_t> > (*fx.getRange());
            fxvec->zero();
            fx.forward(false,_d,fxvec);
            std::shared_ptr<vecReg<data_t> > spec=fxvec->modulus();
            data_t max_spec=spec->max();
            _L->_par.epsilon *= max_spec;

            tr2trDecon decon(*_d->getHyper(),_L->_par.wmin,_L->_par.wmax, _L->_par.epsilon, _L->_par.smth_half_length);
            std::shared_ptr<vecReg<data_t> > output = std::make_shared<vecReg<data_t> >(*decon.getRange());
            output->zero();
            decon.forward(false,_d,output);
            _d = output;
        }

        if (L->_par.normalize_obj_func){
            _dnorm = _d->norm2();
            _dnorm *= _d->getHyper()->getAxis(1).d;
            if (_dnorm<ZERO) _dnorm=1;

            _mnorm = _Dmp->norm2();
            if (_mnorm<ZERO) _mnorm=1;
            //_mnorm = std::max((data_t)1.0, _mnorm);
        }
        else{
            _dnorm=1.0;
            _mnorm=1.0;
        }

        _scale_source_times=0;
        if (_L->_par.scale_source_times>0) {
            for (int s=0; s<_L->_par.ns; s++) _scalers.push_back(1.0);
        }
    }

    void initRes() {
        _r = std::make_shared<vecReg<data_t> >(hypercube<data_t>(_d->getN123()+_D->getRange()->getN123()));
        _r->zero();
        _dfunc.clear();
        _mfunc.clear();
    }

    // res = ( L(m) - d ; lambda*(D(m) - D(m_prior)) )
    // grad = (dL(m)/dm)'.Ht.(L(m)-d)/(d'.Ht.d) + lambda.(dD(m)/dm)'.(D(m)-D(m_prior))/|D(m_prior)|^2
    void res(){
        _r->zero();
        int nd = _d->getN123();
        int nm = _D->getRange()->getN123();
        data_t * pr = _r->getVals();
        const data_t * pd = _d->getCVals();
        const data_t * pdmp = _Dmp->getCVals();
        compute_res_and_grad(_r->getVals());
#ifdef ENABLE_MPI
        // gather all residual
        data_t * temp = new data_t[_d->getN123()];
        memset(temp, 0, _d->getN123()*sizeof(data_t));
        if (sizeof(data_t)==8) MPI_Allreduce(_r->getVals(), temp, _d->getN123(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        else MPI_Allreduce(_r->getVals(), temp, _d->getN123(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        memcpy(_r->getVals(), temp, _d->getN123()*sizeof(data_t));
        delete [] temp;
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        _D->apply_forward(false,_m->getCVals(),pr+nd);
        #pragma omp parallel for
        for (int i=0; i<nm; i++) pr[nd+i] = _lambda*(pr[nd+i] - pdmp[i]);

        _flag = true;
    }

    data_t getFunc(){
        
        int size=1, rank=0;
#ifdef ENABLE_MPI
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (_flag)
        {
            int n123 = _d->getN123();
            int nt = _d->getHyper()->getAxis(1).n;
            int nx = n123 / nt;
            data_t dt = _d->getHyper()->getAxis(1).d;
            data_t * pr = _r->getVals();
            data_t fd = _r->norm2(0,n123);
            data_t fm = _r->norm2(n123,-1);
            for (int i=0; i<nx; i++) fd = fd - 0.5*(pr[i*nt]*pr[i*nt] + pr[i*nt+nt-1]*pr[i*nt+nt-1]);
            fd = 0.5*dt*fd/_dnorm;
            fm = 0.5*fm/_mnorm;
            _f = fd+fm;
            _flag = false;
            
            if (_L->_par.verbose>0 && rank==0) fprintf(stderr,"Data functional = %f; Model functional = %f\n",fd,fm);
            _dfunc.push_back(fd);
            _mfunc.push_back(fm);
        }
        return _f;
    }
    void grad() {
        // the first component is already computed in the compute_res_and_grad() method above
        _D->apply_jacobianT(false,_dg->getVals(),_m->getCVals(),_r->getVals()+_d->getN123());
        
        if (_P != nullptr) _P->jacobianT(false,_g,_m,_pg);
        _g->scaleAdd(_dg,1,_lambda/_mnorm);

        if (_gmask != nullptr) _g->mult(_gmask);
    }
};

#undef ZERO
#undef M_INF