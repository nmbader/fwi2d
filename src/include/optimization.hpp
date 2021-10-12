#pragma once

#define ZERO 1e-16
#define M_INF -1e+16

#include "operator.hpp"

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
    llsq_reg(loper * L, loper * D, std::shared_ptr<vecReg<data_t> > m, std::shared_ptr<vecReg<data_t> > d, const data_t lambda, const std::shared_ptr<vecReg<data_t> > mp){
       successCheck(L->checkDomainRange(m,d),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
       successCheck(D->checkDomainRange(m,m),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        successCheck(m->getN123()==mp->getN123(),__FILE__,__LINE__,"Model and prior model do not have the same number of samples\n");
        _L = L;
        _D = D;
        _m = m;
        _d = d;
        _lambda = lambda;    
        _Dmp = std::make_shared<vecReg<data_t> >(*m->getHyper());
        _Dmp->zero();
        _D->apply_forward(false,mp->getCVals(),_Dmp->getVals());
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

// Non-linear least-squares problem with time quadrature, data weighting and gradient mask for FWI: f(m)=1/2.(L(m)-d)'.W'.Ht.W.(L(m)-d)
class nlls_fwi : public optimization{
protected:
    nloper * _L; // non-linear WE operator
    std::shared_ptr<vecReg<data_t> > _gmask; // gradient mask vector
    std::shared_ptr<vecReg<data_t> > _w; // residual weighting vector

public:
    nlls_fwi(){}
    virtual ~nlls_fwi(){}
    nlls_fwi(nloper * L, std::shared_ptr<vecReg<data_t> > m, std::shared_ptr<vecReg<data_t> > d, std::shared_ptr<vecReg<data_t> > gmask = nullptr, std::shared_ptr<vecReg<data_t> > w = nullptr){
       successCheck(L->checkDomainRange(m,d),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        if (gmask != nullptr) {
            successCheck(m->getN123()==gmask->getN123(),__FILE__,__LINE__,"The gradient mask and model vectors must have the same size\n");
            successCheck(gmask->min()>=0,__FILE__,__LINE__,"Gradient mask must be non-negative\n");
        }
        if (w != nullptr) successCheck(d->getN123()==w->getN123(),__FILE__,__LINE__,"The residual weighting and data vectors must have the same size\n");
        _L = L;
        _m = m;
        _d = d;
        _gmask = gmask;
        _w = w;
    }

    // W.(L(m) - d)
    void res(){
        _L->forward(false,_m,_r);
        _r->scaleAdd(_d,1,-1);
        if (_w != nullptr) _r->mult(_w);
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
    // (dL(m)/dm)'.W'.Ht.W.(Lm-d)
    void grad() {
        int n123 = _r->getN123();
        int nt = _r->getHyper()->getAxis(1).n;
        int nx = n123 / nt;
        data_t dt = _r->getHyper()->getAxis(1).d;
        std::shared_ptr<vecReg<data_t> > r = _r->clone();
        data_t * pr = r->getVals();
        applyHt(false, false, pr, pr, nx, nt, dt, 0, nx);
        if (_w != nullptr) r->mult(_w);
        _L->jacobianT(false,_g,_m,r);
        if (_gmask != nullptr) _g->mult(_gmask);
    }
};

#undef ZERO
#undef M_INF