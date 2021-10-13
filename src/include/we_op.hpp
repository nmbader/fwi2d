#pragma once

#include "operator.hpp"
#include "injector.hpp"
#include "param.hpp"

// check how many wavelets are provided and how many are needed ; duplicate when necessary
std::shared_ptr<vecReg<data_t> > analyzeWavelet(std::shared_ptr<vecReg<data_t> > src, const param &par, bool verbose=true);

// Analyze model, sources,, and receiver geomtery, as well as boundary conditions
void analyzeGeometry(const hypercube<data_t> &model, param &par, bool verbose=true);

// Analyze model and modify as necessary
void analyzeModel(const hypercube<data_t> &domain, std::shared_ptr<vecReg<data_t> > model, param &par);

// non-linear elastic wave equation operator: taking an elastic model and computing data
class nl_we_op_e : virtual public nloper {
    
public:
    param _par;
    std::shared_ptr<vecReg<data_t> > _allsrc;
    std::shared_ptr<vecReg<data_t> > _full_wfld;

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
    nl_we_op_e * clone() const {
        param par = _par;
        nl_we_op_e * op = new nl_we_op_e(_domain,_allsrc,par);
        return op;
    }
    
    // convert Vp, Vs, rho to lambda, mu, rho and vice versa
    // mu = rho.vs2
    // lambda = rho.(vp2 - 2.vs2)
    void convert_model(data_t * m, int n, bool forward) const;
    void compute_gradients(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const;
    void propagate(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz) const;

    void forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        apply_forward(add, mod->getCVals(), dat->getVals());
    }
    virtual void jacobianT(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > mod0, const std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        successCheck(mod->getHyper()->isCompatible(*mod0->getHyper()),__FILE__,__LINE__,"Input and background models hypercubes are not compatible\n");
        apply_jacobianT(add, mod->getVals(), mod0->getCVals(), dat->getCVals());
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
};

// linear elastic wave equation operator: taking a source function and computing data
class l_we_op_e : public loper, public nl_we_op_e {
    
public:
    std::shared_ptr<vecReg<data_t> > _model;

    l_we_op_e(){}
    virtual ~l_we_op_e(){}
    l_we_op_e(const hypercube<data_t> &domain, std::shared_ptr<vecReg<data_t> > model, param &par){
        _model = model->clone();
        _par = par;
        convert_model(_model->getVals(), _model->getN123()/3, true);
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
        convert_model(model->getVals(), model->getN123()/3, false);
        l_we_op_e * op = new l_we_op_e(_domain,model,par);
        return op;
    }

    void setModel(std::shared_ptr<vecReg<data_t> > model){
        successCheck(*model->getHyper()==*_model->getHyper(),__FILE__,__LINE__,"Models must have the same hypercube\n");
        _model = model->clone();
        analyzeModel(_domain,_model,_par);
        convert_model(_model->getVals(), _model->getN123()/3, true);
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