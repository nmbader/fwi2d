#pragma once

#include "operator.hpp"
#include "injector.hpp"
#include "param.hpp"

// Analyze model, sources,, and receiver geomtery, as well as boundary conditions
void analyzeGeometry(const hypercube<data_t> &model, param &par, bool verbose=true);

// check how many wavelets are provided and how many are needed ; duplicate when necessary
std::shared_ptr<vecReg<data_t> > analyzeWavelet(std::shared_ptr<vecReg<data_t> > src, const param &par, bool verbose=true);

// Analyze model and modify as necessary
void analyzeModel(const hypercube<data_t> &domain, std::shared_ptr<vecReg<data_t> > model, param &par);

// general non-linear wave equation operator: taking a model and computing data
class nl_we_op : virtual public nloper {
public:
    param _par;
    std::shared_ptr<vecReg<data_t> > _allsrc;
    std::shared_ptr<vecReg<data_t> > _full_wfld;
    std::shared_ptr<vecReg<data_t> > _full_wflda; // acoustic wavefield for acoustic-elastic coupled medium

    nl_we_op(){}
    virtual ~nl_we_op(){}
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        
        successCheck( domain == _domain,__FILE__,__LINE__,"New domain must be the same as the original\n");
        successCheck( range.getAxis(1) == _range.getAxis(1),__FILE__,__LINE__,"New range must have the same time axis as the original\n");
        successCheck( range.getN123()/range.getAxis(1).n == _range.getN123()/_range.getAxis(1).n, __FILE__,__LINE__, "New range must have the same number of traces as the original\n");

        _domain=domain;
        _range=range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((*mod->getHyper() != _domain) || ( dat->getHyper()->getAxis(1) != _range.getAxis(1)) || ( dat->getN123()/dat->getHyper()->getAxis(1).n != _range.getN123()/_range.getAxis(1).n) ) return false;
        else return true;
    }
    
    virtual void forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        apply_forward(add, mod->getCVals(), dat->getVals());
    }
    virtual void jacobianT(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > mod0, const std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        successCheck(mod->getHyper()->isCompatible(*mod0->getHyper()),__FILE__,__LINE__,"Input and background models hypercubes are not compatible\n");
        apply_jacobianT(add, mod->getVals(), mod0->getCVals(), dat->getCVals());
    }
    virtual void apply_forward(bool add, const data_t * pmod, data_t * pdat) = 0;
    virtual void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) = 0;
};

// non-linear isotropic elastic wave equation operator: taking an elastic model and computing data
class nl_we_op_e : virtual public nl_we_op {
public:
    nl_we_op_e(){}
    virtual ~nl_we_op_e(){}
    nl_we_op_e(const hypercube<data_t> &domain, const std::shared_ptr<vecReg<data_t> > allsrc, param &par){
        _allsrc = allsrc;
        _domain = domain;
        _par = par;
        axis<data_t> X(par.nr,0,1);
        axis<data_t> C(2,0,1);
        if (par.gl>0) C.n = 1;
        _range = hypercube<data_t>(allsrc->getHyper()->getAxis(1),X,C);
    }
    nl_we_op_e * clone() const {
        param par = _par;
        nl_we_op_e * op = new nl_we_op_e(_domain,_allsrc,par);
        return op;
    }
    
    // convert Vp, Vs, rho to lambda, mu, rho and vice versa
    // mu = rho.vs2
    // lambda = rho.(vp2 - 2.vs2)
    virtual void convert_model(data_t * m, int n, bool forward) const;
    // refer to the SBP notes for gradients expression
    virtual void compute_gradients(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const;
    virtual void propagate(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const;

    virtual void compute_gradients_gpu(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const;
    virtual void propagate_gpu(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const;

    virtual void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    virtual void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
};

// linear isotropic elastic wave equation operator: taking a source function and computing data
class l_we_op_e : virtual public loper, virtual public nl_we_op_e {
    
public:
    std::shared_ptr<vecReg<data_t> > _model;

    l_we_op_e(){}
    virtual ~l_we_op_e(){}
    l_we_op_e(const hypercube<data_t> &domain, std::shared_ptr<vecReg<data_t> > model, param &par){
        _model = model->clone();
        _par = par;
        convert_model(_model->getVals(), _model->getN123()/par.nmodels, true);
        _domain = domain; // domain assumed to have come from the output of analyzeWavelet
        axis<data_t> X(par.nr,0,1);
        axis<data_t> C(2,0,1);
        if (par.gl>0) C.n = 1;
        _range = hypercube<data_t>(domain.getAxis(1),X,C);
        if (par.sub>0) _full_wfld=std::make_shared<vecReg<data_t> > (hypercube<data_t>(_model->getHyper()->getAxis(1),_model->getHyper()->getAxis(2),axis<data_t>(2,0,1), axis<data_t>(1+par.nt/par.sub,0,par.dt*par.sub)));
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        
        successCheck( domain.getAxis(1) == _domain.getAxis(1),__FILE__,__LINE__,"New domain must have the same time axis as the original\n");
        successCheck( domain.getN123()/domain.getAxis(1).n == _domain.getN123()/_domain.getAxis(1).n, __FILE__,__LINE__, "New domain must have the same number of traces as the original\n");
        successCheck( range.getAxis(1) == _range.getAxis(1),__FILE__,__LINE__,"New range must have the same time axis as the original\n");
        successCheck( range.getN123()/range.getAxis(1).n == _range.getN123()/_range.getAxis(1).n, __FILE__,__LINE__, "New range must have the same number of traces as the original\n");

        _domain=domain;
        _range=range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((mod->getHyper()->getAxis(1) != _domain.getAxis(1)) || ( mod->getN123()/mod->getHyper()->getAxis(1).n != _domain.getN123()/_domain.getAxis(1).n) 
            || ( dat->getHyper()->getAxis(1) != _range.getAxis(1)) || ( dat->getN123()/dat->getHyper()->getAxis(1).n != _range.getN123()/_range.getAxis(1).n) ) return false;
        else return true;
    }
    l_we_op_e * clone() const {
        param par = _par;
        std::shared_ptr<vecReg<data_t> > model = _model->clone();
        convert_model(model->getVals(), model->getN123()/par.nmodels, false);
        l_we_op_e * op = new l_we_op_e(_domain,model,par);
        return op;
    }

    void setModel(std::shared_ptr<vecReg<data_t> > model){
        successCheck(*model->getHyper()==*_model->getHyper(),__FILE__,__LINE__,"Models must have the same hypercube\n");
        _model = model->clone();
        analyzeModel(_domain,_model,_par);
        convert_model(_model->getVals(), _model->getN123()/_par.nmodels, true);
    }

    void forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<vecReg<data_t> > dat){
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        loper::forward(add, mod, dat);
    }
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {
        loper::apply_jacobianT(add, pmod, pmod0, pdat);
    }
    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// non-linear VTI elastic wave equation operator: taking an elastic model and computing data
class nl_we_op_vti : virtual public nl_we_op_e {
    
public:
    using nl_we_op_e::nl_we_op_e;
    ~nl_we_op_vti(){}
    nl_we_op_vti * clone() const {
        param par = _par;
        nl_we_op_vti * op = new nl_we_op_vti(_domain,_allsrc,par);
        return op;
    }
    
    // convert Vp, Vs, rho, delta, epsilon to generalized lambda, generalized mu, rho, c13, eps and vice versa
    // c33 = lambda+2.mu
    // c11 = (1+2.eps).c33
    // c55 = mu
    // c13 = sqrt[2.c33.(c33-c55).del + (c33-c55)^2] - c55
    // mu = rho.vs2
    // lambda = rho.(vp2 - 2.vs2)
    void convert_model(data_t * m, int n, bool forward) const;
    // refer to the SBP notes for gradients expression
    void compute_gradients(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const;
    void propagate(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const;

    void compute_gradients_gpu(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const;
    void propagate_gpu(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const;

};

// linear isotropic elastic wave equation operator: taking a source function and computing data
class l_we_op_vti : public l_we_op_e, public nl_we_op_vti {
    
public:
    l_we_op_vti(){}
    virtual ~l_we_op_vti(){}
    l_we_op_vti(const hypercube<data_t> &domain, std::shared_ptr<vecReg<data_t> > model, param &par){
        _model = model->clone();
        _par = par;
        nl_we_op_vti::convert_model(_model->getVals(), _model->getN123()/par.nmodels, true);
        _domain = domain; // domain assumed to have come from the output of analyzeWavelet
        axis<data_t> X(par.nr,0,1);
        axis<data_t> C(2,0,1);
        if (par.gl>0) C.n = 1;
        _range = hypercube<data_t>(domain.getAxis(1),X,C);
        if (par.sub>0) _full_wfld=std::make_shared<vecReg<data_t> > (hypercube<data_t>(_model->getHyper()->getAxis(1),_model->getHyper()->getAxis(2),axis<data_t>(2,0,1), axis<data_t>(1+par.nt/par.sub,0,par.dt*par.sub)));
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        l_we_op_e::setDomainRange(domain, range);
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        return l_we_op_e::checkDomainRange(mod, dat);
    }
    l_we_op_vti * clone() const {
        param par = _par;
        std::shared_ptr<vecReg<data_t> > model = _model->clone();
        nl_we_op_vti::convert_model(model->getVals(), model->getN123()/par.nmodels, false);
        l_we_op_vti * op = new l_we_op_vti(_domain,model,par);
        return op;
    }

    void setModel(std::shared_ptr<vecReg<data_t> > model){
        successCheck(*model->getHyper()==*_model->getHyper(),__FILE__,__LINE__,"Models must have the same hypercube\n");
        _model = model->clone();
        analyzeModel(_domain,_model,_par);
        nl_we_op_vti::convert_model(_model->getVals(), _model->getN123()/_par.nmodels, true);
    }

    void forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<vecReg<data_t> > dat){
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        l_we_op_e::forward(add, mod, dat);
    }
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {
        l_we_op_e::apply_jacobianT(add, pmod, pmod0, pdat);
    }
    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){
        l_we_op_e::apply_forward(add, pmod, pdat);
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        l_we_op_e::apply_adjoint(add, pmod, pdat);
    }
};

// non-linear acoustic wave equation operator: taking an acoustic model and computing data
class nl_we_op_a : virtual public nl_we_op_e {
public:
    std::shared_ptr<vecReg<data_t> > _transmission; // vector containing transmission coefficients (in [0,1]) of the top boundary
    nl_we_op_a(){}
    virtual ~nl_we_op_a(){}
    nl_we_op_a(const hypercube<data_t> &domain, const std::shared_ptr<vecReg<data_t> > allsrc, param &par){
        _allsrc = allsrc;
        _domain = domain;
        _par = par;
        axis<data_t> X(par.nr,0,1);
        axis<data_t> C(1,0,1);
        _range = hypercube<data_t>(allsrc->getHyper()->getAxis(1),X,C);
    }
    nl_we_op_a * clone() const {
        param par = _par;
        nl_we_op_a * op = new nl_we_op_a(_domain,_allsrc,par);
        return op;
    }
    
    // convert Vp, rho to bulk modulus K, 1/rho and vice versa
    void convert_model(data_t * m, int n, bool forward) const;

    // refer to the SBP notes for gradients expression
    void compute_gradients(const data_t * model, const data_t * u_full, const data_t * curr, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const;
    void propagate(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const;

    //void compute_gradients_gpu(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const;
    //void propagate_gpu(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const;
};

// linear acoustic wave equation operator: taking a source function and computing data
class l_we_op_a : public l_we_op_e, public nl_we_op_a {
    
public:
    l_we_op_a(){}
    virtual ~l_we_op_a(){}
    l_we_op_a(const hypercube<data_t> &domain, std::shared_ptr<vecReg<data_t> > model, param &par){
        _model = model->clone();
        _par = par;
        nl_we_op_a::convert_model(_model->getVals(), _model->getN123()/par.nmodels, true);
        _domain = domain; // domain assumed to have come from the output of analyzeWavelet
        axis<data_t> X(par.nr,0,1);
        axis<data_t> C(1,0,1);
        _range = hypercube<data_t>(domain.getAxis(1),X,C);
        if (par.sub>0) _full_wfld=std::make_shared<vecReg<data_t> > (hypercube<data_t>(_model->getHyper()->getAxis(1),_model->getHyper()->getAxis(2), axis<data_t>(1+par.nt/par.sub,0,par.dt*par.sub)));
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        l_we_op_e::setDomainRange(domain, range);
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        return l_we_op_e::checkDomainRange(mod, dat);
    }
    l_we_op_a * clone() const {
        param par = _par;
        std::shared_ptr<vecReg<data_t> > model = _model->clone();
        nl_we_op_a::convert_model(model->getVals(), model->getN123()/par.nmodels, false);
        l_we_op_a * op = new l_we_op_a(_domain,model,par);
        return op;
    }

    void setModel(std::shared_ptr<vecReg<data_t> > model){
        successCheck(*model->getHyper()==*_model->getHyper(),__FILE__,__LINE__,"Models must have the same hypercube\n");
        _model = model->clone();
        analyzeModel(_domain,_model,_par);
        nl_we_op_a::convert_model(_model->getVals(), _model->getN123()/_par.nmodels, true);
    }

    void forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<vecReg<data_t> > dat){
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        l_we_op_e::forward(add, mod, dat);
    }
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {
        l_we_op_e::apply_jacobianT(add, pmod, pmod0, pdat);
    }
    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){
        l_we_op_e::apply_forward(add, pmod, pdat);
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        l_we_op_e::apply_adjoint(add, pmod, pdat);
    }
};

// non-linear acoustic-elastic wave equation operator: taking an elastic model and computing data
// the acoustic layer is homogenous and lies on top of the elastic one
class nl_we_op_ae : virtual public nl_we_op_e {
public:
    nl_we_op_ae(){}
    virtual ~nl_we_op_ae(){}
    nl_we_op_ae(const hypercube<data_t> &domain, const std::shared_ptr<vecReg<data_t> > allsrc, param &par){
        _allsrc = allsrc;
        _domain = domain;
        _par = par;
        axis<data_t> X(par.nr,0,1);
        axis<data_t> C(3,0,1);
        if (par.gl>0) C.n = 2;
        _range = hypercube<data_t>(allsrc->getHyper()->getAxis(1),X,C);
    }
    nl_we_op_ae * clone() const {
        param par = _par;
        nl_we_op_ae * op = new nl_we_op_ae(_domain,_allsrc,par);
        return op;
    }

    void propagate(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const;

    //void propagate_gpu(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const;
};

// linear acoustic-elastic wave equation operator: taking a source function and computing data
// the acoustic layer is homogenous and lies on top of the elastic one
class l_we_op_ae : public l_we_op_e, public nl_we_op_ae {
    
public:
    l_we_op_ae(){}
    virtual ~l_we_op_ae(){}
    l_we_op_ae(const hypercube<data_t> &domain, std::shared_ptr<vecReg<data_t> > model, param &par){
        _model = model->clone();
        _par = par;
        nl_we_op_e::convert_model(_model->getVals(), _model->getN123()/par.nmodels, true);
        _domain = domain; // domain assumed to have come from the output of analyzeWavelet
        axis<data_t> X(par.nr,0,1);
        axis<data_t> C(3,0,1);
        if (par.gl>0) C.n = 2;
        _range = hypercube<data_t>(domain.getAxis(1),X,C);
        if (par.sub>0) _full_wfld=std::make_shared<vecReg<data_t> > (hypercube<data_t>(_model->getHyper()->getAxis(1),_model->getHyper()->getAxis(2),axis<data_t>(2,0,1), axis<data_t>(1+par.nt/par.sub,0,par.dt*par.sub)));
        if (par.sub>0 && _par.acoustic_elastic && _par.acoustic_wavefield) _full_wflda=std::make_shared<vecReg<data_t> > (hypercube<data_t>(axis<data_t>(_par.nza,0,_model->getHyper()->getAxis(1).d),_model->getHyper()->getAxis(2), axis<data_t>(1+par.nt/par.sub,0,par.dt*par.sub)));
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        l_we_op_e::setDomainRange(domain, range);
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        return l_we_op_e::checkDomainRange(mod, dat);
    }
    l_we_op_ae * clone() const {
        param par = _par;
        std::shared_ptr<vecReg<data_t> > model = _model->clone();
        nl_we_op_e::convert_model(model->getVals(), model->getN123()/par.nmodels, false);
        l_we_op_ae * op = new l_we_op_ae(_domain,model,par);
        return op;
    }

    void setModel(std::shared_ptr<vecReg<data_t> > model){
        successCheck(*model->getHyper()==*_model->getHyper(),__FILE__,__LINE__,"Models must have the same hypercube\n");
        _model = model->clone();
        analyzeModel(_domain,_model,_par);
        nl_we_op_e::convert_model(_model->getVals(), _model->getN123()/_par.nmodels, true);
    }

    void forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<vecReg<data_t> > dat){
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        l_we_op_e::forward(add, mod, dat);
    }
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {
        l_we_op_e::apply_jacobianT(add, pmod, pmod0, pdat);
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){
        l_we_op_e::apply_forward(add, pmod, pdat);
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        l_we_op_e::apply_adjoint(add, pmod, pdat);
    }
};



// Acooustic (constant density) Born operator, taking reflectivity (slowness squared) and generating data
class born_op_a : virtual public loper {
public:
    param _par;
    std::shared_ptr<vecReg<data_t> > _allsrc;
    std::shared_ptr<vecReg<data_t> > _background_wfld;
    std::shared_ptr<vecReg<data_t> > _model; // background model

    born_op_a(){}
    virtual ~born_op_a(){}
    born_op_a(const std::shared_ptr<vecReg<data_t> > model, const std::shared_ptr<vecReg<data_t> > allsrc, const param &par){
        _model = model;
        _allsrc = allsrc;
        _par = par;
        _domain = hypercube<data_t>(model->getHyper()->getAxis(1),model->getHyper()->getAxis(2));
        axis<data_t> X(par.nr,0,1);
        axis<data_t> C(1,0,1);
        _range = hypercube<data_t>(allsrc->getHyper()->getAxis(1),X,C);
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        
        successCheck( domain.isCompatible(_domain),__FILE__,__LINE__,"New domain must be compatible with the original\n");
        successCheck( range.getAxis(1) == _range.getAxis(1),__FILE__,__LINE__,"New range must have the same time axis as the original\n");
        successCheck( range.getN123()/range.getAxis(1).n == _range.getN123()/_range.getAxis(1).n, __FILE__,__LINE__, "New range must have the same number of traces as the original\n");

        _domain=domain;
        _range=range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ( (!_domain.isCompatible(*mod->getHyper())) || ( dat->getHyper()->getAxis(1) != _range.getAxis(1)) || ( dat->getN123()/dat->getHyper()->getAxis(1).n != _range.getN123()/_range.getAxis(1).n) ) return false;
        else return true;
    }
    born_op_a * clone() const {
        born_op_a * op = new born_op_a(_model,_allsrc,_par);
        return op;
    }
    void compute_gradients(const data_t * model, const data_t * u_full, const data_t * curr, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const;
    void propagate(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const;

    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};