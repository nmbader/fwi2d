#pragma once

#include "vecReg.hpp"
#include "misc.hpp"

// generic class for non-linear operators
class nloper {
protected:
    hypercube<data_t> _domain;
    hypercube<data_t> _range;
public:
    nloper(){}
    virtual ~nloper(){}
    nloper(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        _domain = domain;
        _range = range;
    }
    const hypercube<data_t> * getDomain() const {return &_domain;}
    const hypercube<data_t> * getRange() const {return &_range;}
    virtual void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        _domain = domain;
        _range = range;
    }
    virtual bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {return true;}
    bool checkSame(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((*mod->getHyper() != _domain) || (*dat->getHyper() != _range)) return false;
        else return true;
    }
    bool checkCompatible(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((mod->getHyper()->isCompatible(_domain)) && (dat->getHyper()->isCompatible(_range))) return true;
        else return false;
    }
    bool checkN123(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((mod->getN123() != _domain.getN123()) || (dat->getN123() != _range.getN123())) return false;
        else return true;
    }
    // clone the object
    virtual nloper * clone() const = 0;

    // d = F(m)
    virtual void forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        apply_forward(add, mod->getCVals(), dat->getVals());
    }
    // m = invF(d) if relevant
    virtual void inverse(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        apply_inverse(add, mod->getVals(), dat->getCVals());
    }
    // d = (dF/dm0)*m
    virtual void jacobian(bool add, const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > mod0, const std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        successCheck(mod->getHyper()->isCompatible(*mod0->getHyper()),__FILE__,__LINE__,"Input and background models hypercubes are not compatible\n");
        apply_jacobian(add, mod->getCVals(), mod0->getCVals(), dat->getVals());
    }
    //  m = (dF/dm0)' * d
    virtual void jacobianT(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > mod0, const std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        successCheck(mod->getHyper()->isCompatible(*mod0->getHyper()),__FILE__,__LINE__,"Input and background models hypercubes are not compatible\n");
        apply_jacobianT(add, mod->getVals(), mod0->getCVals(), dat->getCVals());
    }
    virtual void apply_forward(bool add, const data_t * pmod, data_t * pdat) = 0;
    virtual void apply_inverse(bool add, data_t * mod, const data_t * dat) {
        successCheck(false,__FILE__,__LINE__,"This function from the abstract non-linear operator must not be called\n");
    }
    virtual void apply_jacobian(bool add, const data_t * pmod, const data_t * pmod0, data_t * pdat) {
        successCheck(false,__FILE__,__LINE__,"This function from the abstract non-linear operator must not be called\n");
    }
    virtual void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {
        successCheck(false,__FILE__,__LINE__,"This function from the abstract non-linear operator must not be called\n");
    }

    // Approximate dot product test
    // Forward linear operator df/dm
    // Adjoint linear operator (df/dm)' implicit in the jacobian (df/dm)'.d
    // (F(m + eps.dm) - F(m))/eps ~= (dF/dm).dm     for eps << 1
    // Perform dot product < (dF/dm).dm , d > ~= < dm , (dF/dm)'.d >
    void dotProduct();
};

// generic class for linear operators
class loper : virtual public nloper {
protected:

public:
    loper(){}
    virtual ~loper(){}
    virtual void adjoint(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) {
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        apply_adjoint(add, mod->getVals(), dat->getCVals());
    }
    virtual void apply_adjoint(bool add, data_t * pmod, const data_t * pdat) {
        successCheck(false,__FILE__,__LINE__,"This function from the abstract linear operator must not be called.\n");
    }
    void apply_jacobian(bool add, const data_t * pmod, const data_t * pmod0, data_t * pdat) {
        apply_forward(add, pmod, pdat);
    }
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {
        apply_adjoint(add, pmod, pdat);
    }
    // dot product test
    void dotProduct();
};

// Chain of two non-linear (or linear) operators L(R(m))
class chainNLOper : public nloper {
protected:
    nloper * _L; // left operator
    nloper * _R; // right operator
public:
    chainNLOper(){}
    virtual ~chainNLOper(){}
    chainNLOper(nloper * L, nloper * R){
        successCheck(R->getRange()->getN123()==L->getDomain()->getN123(),__FILE__,__LINE__,"Number of samples in range of R must be the same as the one in domain of L\n");
        _domain = *R->getDomain();
        _range = *L->getRange();
        _L=L;
        _R=R;
    }
    nloper * clone() const {
        chainNLOper op(_L, _R);
        return &op;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat) {
        std::shared_ptr<vecReg<data_t> > v = std::make_shared<vecReg<data_t> >(*_R->getRange());
        v->zero();
        _R->apply_forward(false, pmod, v->getVals());
        _L->apply_forward(add, v->getVals(), pdat);
    }
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {
        std::shared_ptr<vecReg<data_t> > m0 = std::make_shared<vecReg<data_t> >(*_L->getRange());
        std::shared_ptr<vecReg<data_t> > v = std::make_shared<vecReg<data_t> >(*_L->getDomain());
        m0->zero();        
        v->zero();
        _R->apply_forward(false, pmod0, m0->getVals());
        _L->apply_jacobianT(false, v->getVals(), m0->getCVals(), pdat);
        _R->apply_jacobianT(add, pmod, pmod0, v->getCVals());
    }
};

// Chain of two linear operators L.R.m
class chainLOper : public loper {
protected:
    loper * _L; // left operator
    loper * _R; // right operator
public:
    chainLOper(){}
    virtual ~chainLOper(){}
    chainLOper(loper * L, loper * R){
        successCheck(R->getRange()->getN123()==L->getDomain()->getN123(),__FILE__,__LINE__,"Number of samples in range of R must be the same as the one in domain of L\n");
        _domain = *R->getDomain();
        _range = *L->getRange();
        _L=L;
        _R=R;
    }
    nloper * clone() const {
        chainLOper op(_L, _R);
        return &op;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat) {
        std::shared_ptr<vecReg<data_t> > v = std::make_shared<vecReg<data_t> >(*_R->getRange());
        v->zero();
        _R->apply_forward(false, pmod, v->getVals());
        _L->apply_forward(false, v->getVals(), pdat);
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat) {
        std::shared_ptr<vecReg<data_t> > v = std::make_shared<vecReg<data_t> >(*_R->getRange());
        v->zero();
        _L->apply_adjoint(false, v->getVals(), pdat);
        _R->apply_adjoint(add, pmod, v->getCVals());
    }
};

// operator to resampler data in time (along the fast axis)
class resampler : public loper {
public:
    resampler(){}
    virtual ~resampler(){}
    resampler(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        setDomainRange(domain,range);
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        
        successCheck( domain.getAxis(1).o == range.getAxis(1).o,__FILE__,__LINE__,"Domain and range must have the same time origin\n");
        successCheck( domain.getN123()/domain.getAxis(1).n == range.getN123()/range.getAxis(1).n, __FILE__,__LINE__, "Domain and range must have the same number of traces\n");

        _domain=domain;
        _range=range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((mod->getHyper()->getAxis(1).o != dat->getHyper()->getAxis(1).o) || (mod->getN123()/mod->getHyper()->getAxis(1).n != dat->getN123()/dat->getHyper()->getAxis(1).n)) return false;
        else return true;
    }
}; 

// linear resampler
class linear_resampler : public resampler{
public:
    // inherit constructor from the base class
    using resampler::resampler;
    ~linear_resampler(){}
    linear_resampler * clone() const{
        linear_resampler * resampler = new linear_resampler(_domain, _range);
        return resampler;
    }

    // forward operator (often interpolation)
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    // adjoint operator (often decimation)
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);

};

// sinc resampler
class sinc_resampler : public resampler{
protected:
    int _hl;

public:
    sinc_resampler(){}
    ~sinc_resampler(){}
    sinc_resampler(const hypercube<data_t> &domain, const hypercube<data_t> &range, int half_length):resampler(domain, range){
   
        successCheck(half_length>0,__FILE__,__LINE__,"The half length of the sinc filter must be larger than 0\n");
        _hl = half_length;
    }
    sinc_resampler * clone() const{
        sinc_resampler * resampler = new sinc_resampler(_domain, _range, _hl);
        return resampler;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// Matrix multiplication by a vector; the matrix is stored in row major
class matrix : public loper {
protected:
    std::shared_ptr<vecReg<data_t> > _mat;
public:
    matrix(){}
    ~matrix(){}
    matrix(const std::shared_ptr<vecReg<data_t> > mat, const hypercube<data_t> domain, const hypercube<data_t> range){
        successCheck(mat->getHyper()->getNdim()==2,__FILE__,__LINE__,"The matrix must be 2D\n");
        successCheck(mat->getHyper()->getAxis(1).n==domain.getN123(),__FILE__,__LINE__,"The domain dimension must match the matrix row size\n");
        successCheck(mat->getHyper()->getAxis(2).n==range.getN123(),__FILE__,__LINE__,"The range dimension must match the matrix column size\n");
        _domain = domain;
        _range = range;
        _mat=mat->clone();
    }
    matrix(const std::shared_ptr<vecReg<data_t> > mat){
        successCheck(mat->getHyper()->getNdim()==2,__FILE__,__LINE__,"The matrix must be 2D\n");
        _mat=mat->clone();
        _domain = hypercube<data_t>(mat->getHyper()->getAxis(1));
        _range = hypercube<data_t>(mat->getHyper()->getAxis(2));
    }
    matrix(int n, data_t val){
        _domain = hypercube<data_t>(n);
        _range = hypercube<data_t>(n);
        _mat = std::make_shared<vecReg<data_t> >(hypercube<data_t>(n,n));
        for (int i=0; i<n; i++) _mat->getVals()[i*n+i] = val;
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(_domain.getN123()==domain.getN123(),__FILE__,__LINE__,"The domain dimension must match the matrix row size\n");
        successCheck(_range.getN123()==range.getN123(),__FILE__,__LINE__,"The range dimension must match the matrix column size\n");
        _domain=domain;
        _range=range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        if ((_domain.getN123() != mod->getN123()) || (_range.getN123() != dat->getN123())) return false;
        else return true;
    }
    matrix * clone() const {
        matrix * mat = new matrix(_mat,_domain,_range);
        return mat;
    }
    std::shared_ptr<vecReg<data_t> > getMat() const {return _mat;}
    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){
        const data_t * pmat = _mat->getCVals();
        int n=_domain.getN123();
        int m=_range.getN123();
        data_t val;
        #pragma omp parallel for
        for (int i=0; i<m; i++){
            val=0;
            for (int j=0; j<n; j++){
                val += pmat[i*m+j]*pmod[j]; 
            }
            pdat[i] = add*pdat[i] + val;
        } 
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        const data_t * pmat = _mat->getCVals();
        int n=_domain.getN123();
        int m=_range.getN123();
        data_t val;
        #pragma omp parallel for
        for (int j=0; j<n; j++){
            val=0;
            for (int i=0; i<m; i++){
                val += pmat[i*m+j]*pdat[i]; 
            }
            pmod[j] = add*pmod[j] + val;
        }
    }
};