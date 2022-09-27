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
    virtual void apply_forward(bool add, const data_t * pmod, data_t * pdat) {
        successCheck(false,__FILE__,__LINE__,"This function from the abstract non-linear operator must not be called\n");
    }
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
    virtual loper * clone() const = 0;
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

// identity operator
class identity : public loper {
protected:

public:
    identity(){}
    ~identity(){}
    identity(const hypercube<data_t> &domain){
        _domain = domain;
        _range = domain;
    }
    identity * clone() const {
        identity * op = new identity(_domain);
        return op;
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain.getN123()==_domain.getN123(),__FILE__,__LINE__,"The new domain must have the same number of samples\n");
        successCheck(range.getN123()==_range.getN123(),__FILE__,__LINE__,"The new range must have the same number of samples\n");
        successCheck(domain.getN123()==range.getN123(),__FILE__,__LINE__,"The new domain and range must have the same number of samples\n");
        _domain = domain;
        _range = range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        bool ans1 = (mod->getN123()==_domain.getN123());
        bool ans2 = (dat->getN123()==_range.getN123());
        bool ans3 = (dat->getN123()==mod->getN123());
        return (ans1 && ans2 && ans3);
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat) {
        for (int i=0; i<_domain.getN123(); i++) pdat[i] = add*pdat[i] + pmod[i];
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat) {
        for (int i=0; i<_domain.getN123(); i++) pmod[i] = add*pmod[i] + pdat[i];
    }
};

// Chain of two non-linear (or linear) operators L(R(m))
class chainNLOper : public nloper {
protected:
    nloper * _L; // left operator
    nloper * _R; // right operator
    std::shared_ptr<vecReg<data_t> > _m0; // intermediate vectors used for chained forward and jacobian
    std::shared_ptr<vecReg<data_t> > _v;

public:
    chainNLOper(){}
    ~chainNLOper(){delete _L; delete _R;}
    chainNLOper(nloper * L, nloper * R){
        successCheck(R->getRange()->getN123()==L->getDomain()->getN123(),__FILE__,__LINE__,"Number of samples in range of R must be the same as the one in domain of L\n");
        _domain = *R->getDomain();
        _range = *L->getRange();
        _L=L->clone();
        _R=R->clone();
        _m0 = std::make_shared<vecReg<data_t> >(*_R->getRange());
        _v = std::make_shared<vecReg<data_t> >(*_R->getRange());
        _m0->zero();
        _v->zero();
    }
    chainNLOper * clone() const {
        chainNLOper * op = new chainNLOper(_L, _R);
        return op;
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain.getN123()==_domain.getN123(),__FILE__,__LINE__,"The new domain must have the same number of samples\n");
        successCheck(range.getN123()==_range.getN123(),__FILE__,__LINE__,"The new range must have the same number of samples\n");
        _domain = domain;
        _range = range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        bool ans1 = (mod->getN123()==_domain.getN123());
        bool ans2 = (dat->getN123()==_range.getN123());
        return (ans1 && ans2);
    }    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat) {
        _R->apply_forward(false, pmod, _v->getVals());
        _L->apply_forward(add, _v->getCVals(), pdat);
    }
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {
        _R->apply_forward(false, pmod0, _m0->getVals());
        _L->apply_jacobianT(false, _v->getVals(), _m0->getCVals(), pdat);
        _R->apply_jacobianT(add, pmod, pmod0, _v->getCVals());
    }
};

// Chain of two linear operators L.R.m
class chainLOper : public loper {
protected:
    loper * _L; // left operator
    loper * _R; // right operator
    std::shared_ptr<vecReg<data_t> > _v;
public:
    chainLOper(){}
    ~chainLOper(){delete _L; delete _R;}
    chainLOper(loper * L, loper * R){
        successCheck(R->getRange()->getN123()==L->getDomain()->getN123(),__FILE__,__LINE__,"Number of samples in range of R must be the same as the one in domain of L\n");
        _domain = *R->getDomain();
        _range = *L->getRange();
        _L=L->clone();
        _R=R->clone();
        _v = std::make_shared<vecReg<data_t> >(*_R->getRange());
        _v->zero();
    }
    chainLOper * clone() const {
        chainLOper * op = new chainLOper(_L, _R);
        return op;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat) {
        _R->apply_forward(false, pmod, _v->getVals());
        _L->apply_forward(add, _v->getCVals(), pdat);
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat) {
        _L->apply_adjoint(false, _v->getVals(), pdat);
        _R->apply_adjoint(add, pmod, _v->getCVals());
    }
};

// non-linear soft clip operator h(x) = ig(f(g(x)))
// g(x) = (2.(x-xmean)/Dx)^p    Dx = xmax-xmin  xmean = (xmax+xmin)/2
// f(g) = g-g^q/q if |g|<1 and = sign(g).(q-1)/q otherwise
// ig(f) = Dx/2.f^(1/p) + xmean
// h(x) = Dx/2.(g-g^q/q)^(1/p) + xmean if |g|<1  and = Dx/2.(sign(g).(q-1)/q)^(1/p) + xmean otherwise 
// h'(x) = (2.(x-xmean)/Dx)^(p-1).(1-g^(q-1)).(g-g^q/q)^(1/p-1) if |g|<1 and = 0 otherwise
// p and q are positive odd integers
// xmax and xmin are actually re-scaled so that the final clipping is within the provided range
class softClip : public nloper {
protected:
    data_t _xmin;
    data_t _xmax;
    int _p;
    int _q;
public:
    softClip(){}
    ~softClip(){}
    softClip(const hypercube<data_t> &domain, data_t xmin, data_t xmax, int p=1, int q=9){
        successCheck(xmax>xmin,__FILE__,__LINE__,"The upper bound must be strictly greater than the lower bound\n");
        successCheck((p%2==1) && (q%2==1),__FILE__,__LINE__,"The powers p and q must be odd numbers\n");
        _domain = domain;
        _range = domain;
        _xmin = xmin;
        _xmax = xmax;
        _p = p;
        _q = q;
    }
    softClip * clone() const {
        softClip * op = new softClip(_domain, _xmin, _xmax, _p, _q);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
};

// soft clip of the elastic model and the Vs/Vp ratio
class emodelSoftClip : public nloper {
    data_t _vpmin;
    data_t _vpmax;
    data_t _vsmin;
    data_t _vsmax;
    data_t _rhomin;
    data_t _rhomax;
    data_t _spratio;
    int _p;
    int _q;
public:
    emodelSoftClip(){}
    ~emodelSoftClip(){}
    emodelSoftClip(const hypercube<data_t> &domain, data_t vpmin, data_t vpmax, data_t vsmin, data_t vsmax, data_t rhomin, data_t rhomax, data_t spratio=1/sqrt(2.00001), int p=9, int q=9){
        successCheck((domain.getNdim()==3) && (domain.getAxis(3).n>=3),__FILE__,__LINE__,"The domain must be 3D with 3rd dimension containing at least 3 fields\n");
        _domain = domain;
        _range = domain;
        _vpmin=vpmin; _vpmax=vpmax; _vsmin=vsmin; _vsmax=vsmax; _rhomin=rhomin; _rhomax=rhomax; _spratio=spratio;
        _p = p;
        _q = q;
    }
    emodelSoftClip * clone() const {
        emodelSoftClip * op = new emodelSoftClip(_domain, _vpmin, _vpmax, _vsmin, _vsmax, _rhomin, _rhomax, _spratio, _p, _q);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
};

// Elastic model parameterization using P-impedance, S-impedance, density, then anisotropy if applicable
class pi_si_rho : public nloper {
public:
    pi_si_rho(){}
    ~pi_si_rho(){}
    pi_si_rho(const hypercube<data_t> &domain){
        successCheck((domain.getNdim()>=2) && (domain.getAxis(domain.getNdim()).n>=3),__FILE__,__LINE__,"The domain must be at least 2D with the last dimension containing at least 3 fields\n");
        _domain = domain;
        _range = domain;
    }
    pi_si_rho * clone() const {
        pi_si_rho * op = new pi_si_rho(_domain);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
    void apply_inverse(bool add, data_t * pmod, const data_t * pdat);
};

// Elastic model parameterization using log(Vs/Vs0), log(Vp/Vs - sqrt(2)), log(rho/rho0), then anisotropy if applicable
// Vs0 and rho0 are some fixed reference values
class vs_vpvs_rho : public nloper {
protected:
    data_t _vs0;
    data_t _rho0;
public:
    vs_vpvs_rho(){}
    ~vs_vpvs_rho(){}
    vs_vpvs_rho(const hypercube<data_t> &domain, data_t vs0=1.0, data_t rho0=1.0){
        successCheck((domain.getNdim()>=2) && (domain.getAxis(domain.getNdim()).n>=3),__FILE__,__LINE__,"The domain must be at least 2D with the last dimension containing at least 3 fields\n");
        _domain = domain;
        _range = domain;
        _vs0=vs0;
        _rho0=rho0;
    }
    vs_vpvs_rho * clone() const {
        vs_vpvs_rho * op = new vs_vpvs_rho(_domain,_vs0,_rho0);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
    void apply_inverse(bool add, data_t * pmod, const data_t * pdat);
};

// class to transform from x+iy to r.exp(theta) (all with real numbers)
class polar : public nloper {
public:
    polar(){}
    ~polar(){}
    polar(const hypercube<data_t> &domain){
        successCheck(2*round(domain.getN123()/2)==domain.getN123(),__FILE__,__LINE__,"The number of samples must be even\n");
        _domain = domain;
        _range = domain;
    }
    polar * clone() const {
        polar * op = new polar(_domain);
        return op;
    }   
    // parameter add is ignored
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
};

// class doing the inverse of polar
class ipolar : public nloper {
public:
    ipolar(){}
    ~ipolar(){}
    ipolar(const hypercube<data_t> &domain){
        successCheck(2*round(domain.getN123()/2)==domain.getN123(),__FILE__,__LINE__,"The number of samples must be even\n");
        _domain = domain;
        _range = domain;
    }
    ipolar * clone() const {
        ipolar * op = new ipolar(_domain);
        return op;
    }   
    // parameter add is ignored
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
};

// class to perform spectral division with white noise added
class sdivision : public nloper {
protected:
    data_t _fmin, _fmax, _eps;
public:
    sdivision(){}
    ~sdivision(){}
    sdivision(const hypercube<data_t> &domain, data_t fmin, data_t fmax, data_t eps){
        successCheck(fmax>fmin && fmin>=0 && fmax<=1,__FILE__,__LINE__,"The high and low cutoff must satisfy 0<=fmin<fmax<=1\n");
        successCheck(domain.getNdim()>=2,__FILE__,__LINE__,"The domain must contain at least 2 axes\n");
        _domain = domain;
        _range = domain;
        _fmin = fmin;
        _fmax = fmax;
        _eps = eps;
    }
    sdivision * clone() const {
        sdivision * op = new sdivision(_domain, _fmin, _fmax, _eps);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
};

// Non-linear operator to compute deconvolution between 2 consecutive traces
// The output is the matching filters (same size as the input)
// The deconvolution is performed in the frequency domain
class tr2trDecon : public nloper {
protected:
    data_t _fmin, _fmax; // low and high cut frequencies to be used in the deconvolution ; expressed as fraction of Nyquist
    data_t _eps; // white noise to be added to amplitude spectra before deconvolution
    int _hl; // half length smoothing window size applied to amplitude spectra after deconvolution
public:
    tr2trDecon(){}
    ~tr2trDecon(){}
    tr2trDecon(const hypercube<data_t> &domain, data_t fmin, data_t fmax, data_t eps, int window){
        successCheck(fmax>fmin && fmin>=0 && fmax<=1,__FILE__,__LINE__,"The high and low cutoff must satisfy 0<=fmin<fmax<=1\n");
        successCheck(window>=0,__FILE__,__LINE__,"The half length spectral smoothing window must be >= 0\n");
        successCheck(domain.getNdim()>=2,__FILE__,__LINE__,"The domain must contain at least 2 axes\n");
        _domain = domain;
        _range = domain;
        _fmin = fmin;
        _fmax = fmax;
        _eps = eps;
        _hl = window;
    }
    tr2trDecon * clone() const {
        tr2trDecon * op = new tr2trDecon(_domain, _fmin, _fmax, _eps, _hl);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat);
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
    data_t _alpha;

public:
    sinc_resampler(){}
    ~sinc_resampler(){}
    sinc_resampler(const hypercube<data_t> &domain, const hypercube<data_t> &range, int half_length, data_t alpha=0.5):resampler(domain, range){
   
        successCheck(half_length>0,__FILE__,__LINE__,"The half length of the sinc filter must be larger than 0\n");
        _hl = half_length;
        _alpha = alpha;
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
    matrix(const std::shared_ptr<vecReg<data_t> > mat, const hypercube<data_t> &domain, const hypercube<data_t> &range){
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

// Integral operator along the fast axis using the Trapezoidal quadrature
// The add option is inactive
class integral : public loper {
public:
    integral(){}
    ~integral(){}
    integral(const hypercube<data_t> &domain){
        _domain = domain;
        _range = domain;
    }
    integral * clone() const {
        integral * op = new integral(_domain);
        return op;
    }    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){
        int nt = _domain.getAxis(1).n;
        int nx = _domain.getN123() / nt;
        data_t dt = _domain.getAxis(1).d;
        #pragma omp parallel for
        for (int ix=0; ix<nx; ix++){
            int i1=ix*nt;
            pdat[i1] = 0.5*pmod[i1]*dt;
            for (int it=1; it<nt-1; it++){
                pdat[i1+it] = pdat[i1+it-1] + pmod[i1+it]*dt;
            }
            pdat[i1+nt-1] = pdat[i1+nt-2] + 0.5*pmod[i1+nt-1]*dt;
        }
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        int nt = _domain.getAxis(1).n;
        int nx = _domain.getN123() / nt;
        data_t dt = _domain.getAxis(1).d;
        for (int ix=0; ix<nx; ix++){
            int i1=ix*nt;
            pmod[i1+nt-1] = pdat[i1+nt-1]*dt;
            for (int it=nt-2; it>=0; it--){
                pmod[i1+it] = pmod[i1+it+1] + pdat[i1+it]*dt;
            }
            pmod[i1+nt-1] *= 0.5;
            pmod[i1] *= 0.5; 
        }
    }
};

// operator to take the difference of two consecutive traces in the second direction
// The add option is inactive
class xdifference : public loper {
public:
    xdifference(){}
    ~xdifference(){}
    xdifference(const hypercube<data_t> &domain){
        _domain = domain;
        _range = domain;
    }
    xdifference * clone() const {
        xdifference * op = new xdifference(_domain);
        return op;
    }    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){
        int nt = _domain.getAxis(1).n;
        int nx = _domain.getAxis(2).n;
        int ny = _domain.getN123() / (nt*nx);
        const data_t (*pm) [nx][nt] = (const data_t (*)[nx][nt]) pmod;
        data_t (*pd) [nx][nt] = (data_t (*)[nx][nt]) pdat;
        #pragma omp parallel for
        for (int iy=0; iy<ny; iy++){
            for (int ix=nx-1; ix>0; ix--){
                for (int it=0; it<nt; it++){
                    pd[iy][ix][it] = add*pd[iy][ix][it] + pm[iy][ix][it]-pm[iy][ix-1][it];
                }
            }
            for (int it=0; it<nt; it++) pd[iy][0][it] = 0;
        }
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        int nt = _domain.getAxis(1).n;
        int nx = _domain.getAxis(2).n;
        int ny = _domain.getN123() / (nt*nx);
        data_t (*pm) [nx][nt] = (data_t (*)[nx][nt]) pmod;
        const data_t (*pd) [nx][nt] = (const data_t (*)[nx][nt]) pdat;
        #pragma omp parallel for
        for (int iy=0; iy<ny; iy++){
            for (int it=0; it<nt; it++) {
                pm[iy][0][it] = - pd[iy][1][it];
            }
            for (int ix=1; ix<nx-1; ix++){
                for (int it=0; it<nt; it++){
                    pm[iy][ix][it] = pd[iy][ix][it]-pd[iy][ix+1][it];
                }
            }
            for (int it=0; it<nt; it++) {
                pm[iy][nx-1][it] = pd[iy][nx-1][it];
            }
        }
    }
};

// class to tranform from t-x (or z-x) to f-x space and vice versa
class fxTransform{
protected:
    hypercube<data_t> _domain;
    hypercube<data_t> _range;
    hypercube<data_t> _crange;

public:
	fxTransform(){}
    ~fxTransform(){}
    fxTransform(const hypercube<data_t> &domain){
        _domain = domain;
        std::vector<axis<data_t> > axes = domain.getAxes();
        axis<data_t> Z = axes[0];
        Z.n = Z.n/2 + 1;
        Z.d = 1.0/((domain.getAxis(1).n-1)*Z.d);
        Z.o = 0.0;
        axes[0] = Z;
        _range = hypercube<data_t>(axes);
        Z.n = domain.getAxis(1).n;
        Z.o = -(Z.n-1)/2 * Z.d;
        axes[0] = Z;
        _crange = hypercube<data_t> (axes);
    }
    fxTransform * clone() const{
        fxTransform * op = new fxTransform(_domain);
        return op;
    }
    const hypercube<data_t> * getDomain() const {return &_domain;}
    const hypercube<data_t> * getRange() const {return &_range;}
    const hypercube<data_t> * getCRange() const {return &_crange;}
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain.getN123()==_domain.getN123(),__FILE__,__LINE__,"The new domain must have the same number of samples\n");
        successCheck(range.getN123()==_range.getN123(),__FILE__,__LINE__,"The new range must have the same number of samples\n");
        successCheck(domain.getAxis(1)==_domain.getAxis(1),__FILE__,__LINE__,"The first axis in the domain must be the same\n");
        successCheck(range.getAxis(1)==_range.getAxis(1),__FILE__,__LINE__,"The first axis in the range must be the same\n");
        _domain = domain;
        _range = range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat) const {
        bool ans1 = (mod->getN123()==_domain.getN123());
        bool ans2 = (dat->getN123()==_range.getN123());
        bool ans3 = (mod->getHyper()->getAxis(1)==_domain.getAxis(1));
        bool ans4 = (dat->getHyper()->getAxis(1)==_range.getAxis(1));
        return (ans1 && ans2 && ans3 && ans4);
    }
    bool checkCDomainRange(const std::shared_ptr<cvecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat) const {
        bool ans1 = (mod->getN123()==_domain.getN123());
        bool ans2 = (dat->getN123()==_crange.getN123());
        bool ans3 = (mod->getHyper()->getAxis(1)==_domain.getAxis(1));
        bool ans4 = (dat->getHyper()->getAxis(1)==_crange.getAxis(1));
        return (ans1 && ans2 && ans3 && ans4);
    }
    void apply_forward(bool add, const data_t * pmod, std::complex<data_t> * pdat);
    void apply_adjoint(bool add, data_t * pmod, const std::complex<data_t> * pdat){
        apply_inverse(add, pmod, pdat);
    }
    void apply_inverse(bool add, data_t * pmod, const std::complex<data_t> * pdat);
    void forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat){
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        apply_forward(add, mod->getCVals(), dat->getVals());
    }
    void inverse(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat){
        successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");
        apply_inverse(add, mod->getVals(), dat->getCVals());
    }
    void adjoint(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat){
        inverse(add, mod, dat);
    }
    void cforward(bool add, const std::shared_ptr<cvecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat);
    void cinverse(bool add, std::shared_ptr<cvecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat);
    void cadjoint(bool add, std::shared_ptr<cvecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat){
        cinverse(add, mod, dat);
    }

    void testInverse();
    void testCInverse();
    void testAdjoint();
    void testCAdjoint();
};

// class to tranform from t-x (or z-x) to f-k space and vice versa
class fkTransform : public fxTransform{
protected:

public:
	fkTransform(){}
    ~fkTransform(){}
    fkTransform(const hypercube<data_t> &domain){
        successCheck(domain.getNdim()>1,__FILE__,__LINE__,"The domain must contain at least 2 dimensions\n");
        _domain = domain;
        std::vector<axis<data_t> > axes = domain.getAxes();
        axis<data_t> Z = axes[0];
        axis<data_t> X = axes[1];
        Z.n = Z.n/2 + 1;
        Z.d = 1.0/((domain.getAxis(1).n-1)*Z.d);
        Z.o = 0.0;
        X.d = 2*M_PI/((X.n-1)*X.d);
        X.o = - X.d * floor((X.n-1)/2);
        axes[0] = Z;
        axes[1] = X;
        _range = hypercube<data_t>(axes);
        Z.n = domain.getAxis(1).n;
        Z.o = -(Z.n-1)/2 * Z.d;
        axes[0] = Z;
        _crange = hypercube<data_t> (axes);
    }
    fkTransform * clone() const{
        fkTransform * op = new fkTransform(_domain);
        return op;
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain.getN123()==_domain.getN123(),__FILE__,__LINE__,"The new domain must have the same number of samples\n");
        successCheck(range.getN123()==_range.getN123(),__FILE__,__LINE__,"The new range must have the same number of samples\n");
        successCheck(domain.getAxis(1)==_domain.getAxis(1),__FILE__,__LINE__,"The first axis in the domain must be the same\n");
        successCheck(domain.getAxis(2)==_domain.getAxis(2),__FILE__,__LINE__,"The second axis in the domain must be the same\n");
        successCheck(range.getAxis(1)==_range.getAxis(1),__FILE__,__LINE__,"The first axis in the range must be the same\n");
        successCheck(range.getAxis(2)==_range.getAxis(2),__FILE__,__LINE__,"The second axis in the range must be the same\n");
        _domain = domain;
        _range = range;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat) const {
        bool ans1 = (mod->getN123()==_domain.getN123());
        bool ans2 = (dat->getN123()==_range.getN123());
        bool ans3 = (mod->getHyper()->getAxis(1)==_domain.getAxis(1));
        bool ans4 = (dat->getHyper()->getAxis(1)==_range.getAxis(1));
        bool ans5 = (mod->getHyper()->getAxis(2)==_domain.getAxis(2));
        bool ans6 = (dat->getHyper()->getAxis(2)==_range.getAxis(2));
        return (ans1 && ans2 && ans3 && ans4 && ans5 && ans6);
    }
    bool checkCDomainRange(const std::shared_ptr<cvecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat) const {
        bool ans1 = (mod->getN123()==_domain.getN123());
        bool ans2 = (dat->getN123()==_crange.getN123());
        bool ans3 = (mod->getHyper()->getAxis(1)==_domain.getAxis(1));
        bool ans4 = (dat->getHyper()->getAxis(1)==_crange.getAxis(1));
        bool ans5 = (mod->getHyper()->getAxis(2)==_domain.getAxis(2));
        bool ans6 = (dat->getHyper()->getAxis(2)==_crange.getAxis(2));
        return (ans1 && ans2 && ans3 && ans4 && ans5 && ans6);
    }
    void forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat);
    void inverse(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat);
    void cforward(bool add, const std::shared_ptr<cvecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat);
    void cinverse(bool add, std::shared_ptr<cvecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat);
};

// class to tranform from t-x (or z-x) to f-x space using real numbers only
// in forward mode, real parts are stored in even indices (0, 2, etc..) and imaginary parts in odd indices(1, 3, etc.)
// note that the adjoint of this operator is NOT the inverse Fourier transform
class fxTransformR : public loper {
public:
    fxTransformR(){}
    ~fxTransformR(){}
    fxTransformR(const hypercube<data_t> &domain){
        _domain = domain;
        std::vector<axis<data_t> > axes = domain.getAxes();
        axis<data_t> Z = axes[0];
        Z.n = Z.n/2 + 1;
        Z.n *=2;
        Z.d = 1.0/((domain.getAxis(1).n-1)*Z.d);
        Z.o = 0.0;
        axes[0] = Z;
        _range = hypercube<data_t>(axes);
    }
    fxTransformR * clone() const {
        fxTransformR * op = new fxTransformR(_domain);
        return op;
    } 
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// class to perform the inverse of fxTransformR
class ifxTransformR : public loper {
public:
    ifxTransformR(){}
    ~ifxTransformR(){}
    ifxTransformR(const hypercube<data_t> &range){
        _range = range;
        std::vector<axis<data_t> > axes = range.getAxes();
        axis<data_t> Z = axes[0];
        Z.n = Z.n/2 + 1;
        Z.n *=2;
        Z.d = 1.0/((range.getAxis(1).n-1)*Z.d);
        Z.o = 0.0;
        axes[0] = Z;
        _domain = hypercube<data_t>(axes);
    }
    ifxTransformR * clone() const {
        ifxTransformR * op = new ifxTransformR(_range);
        return op;
    } 
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// low order 2D gradient operator to be used in first order Tikhonov regularization
// the boundary gradient is discarded
class gradient2d : public loper {
protected:
    data_t _xw, _zw; // directional weights
public:
    gradient2d(){}
    ~gradient2d(){}
    gradient2d(const hypercube<data_t> &domain, data_t xw=1, data_t zw=1){
        successCheck(domain.getNdim()>=2,__FILE__,__LINE__,"The domain must contain at least 2 dimensions\n");
        _domain = domain;
        std::vector<axis<data_t> > axes = domain.getAxes();
        axes.push_back(axis<data_t>(2,0,1));
        _range = hypercube<data_t>(axes);
        _xw=xw;
        _zw=zw;
    }
    gradient2d * clone() const {
        gradient2d * op = new gradient2d(_domain);
        return op;
    }    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// low order 2D laplacian operator to be used in second order Tikhonov regularization
// the boundary laplacian is discarded
class laplacian2d : public loper {
protected:
    data_t _xw, _zw; // directional weights
public:
    laplacian2d(){}
    ~laplacian2d(){}
    laplacian2d(const hypercube<data_t> &domain, data_t xw=1, data_t zw=1){
        successCheck(domain.getNdim()>=2,__FILE__,__LINE__,"The domain must contain at least 2 dimensions\n");
        _domain = domain;
        _range = domain;
        _xw=xw;
        _zw=zw;
    }
    laplacian2d * clone() const {
        laplacian2d * op = new laplacian2d(_domain);
        return op;
    }    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// Model extrapolation from 1D to 2D along a given horizon in the second dimension
class extrapolator1d2d : public loper {
protected:
    std::shared_ptr<vecReg<data_t> > _hrz; // 1D horizon for extrapolation
public:
    extrapolator1d2d(){}
    ~extrapolator1d2d(){}
    extrapolator1d2d(const hypercube<data_t> &domain, std::shared_ptr<vecReg<data_t> > hrz){
        successCheck(hrz->getHyper()->getNdim()==1,__FILE__,__LINE__,"The horizon must be 1D\n");
        _domain = domain;
        std::vector<axis<data_t> > axes = domain.getAxes();
        axis<data_t> X = hrz->getHyper()->getAxis(1);
        axes.insert(axes.begin()+1, X);
        _range = hypercube<data_t>(axes);
        _hrz = hrz->clone();
    }
    extrapolator1d2d * clone() const {
        extrapolator1d2d * op = new extrapolator1d2d(_domain, _hrz);
        return op;
    }    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
    void apply_inverse(bool add, data_t * pmod, const data_t * pdat);
};

// class to perform spectral triangular smoothing
class ssmth : public loper {
protected:
    int _hl;
public:
    ssmth(){}
    ~ssmth(){}
    ssmth(const hypercube<data_t> &domain, int hl){
        successCheck(hl>=0,__FILE__,__LINE__,"The half length spectral smoothing window must be >= 0\n");
        _domain = domain;
        _range = domain;
        _hl=hl;
    }
    ssmth * clone() const {
        ssmth * op = new ssmth(_domain,_hl);
        return op;
    } 
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// Time domain Convolution of 1D filter with N-dimensional vector
class conv1dnd : public loper {
protected:
    std::shared_ptr<vecReg<data_t> > _f; // filter always treated as 1D
//    bool _frequency; // perform convolution in time or frequency domain
    bool _centered; // center the convolution

public:
    conv1dnd(){}
    ~conv1dnd(){}
    conv1dnd(const hypercube<data_t> domain, const std::shared_ptr<vecReg<data_t> > f, bool centered = true){
        successCheck(f->getHyper()->getAxis(1).d==domain.getAxis(1).d,__FILE__,__LINE__,"Filter must have the same sampling as the domain fast axis\n");
        _domain=domain;
        _range=domain;
        _f=f->clone();
        _centered = centered;
    }
    conv1dnd * clone() const {
        conv1dnd * op = new conv1dnd(_domain,_f,_centered);
        return op;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// Time domain Convolution of N-dimensional data with 1D filter
// domain always treated as 1D
class convnd1d : public loper {
protected:
    std::shared_ptr<vecReg<data_t> > _f; // N-dimensional data
//    bool _frequency; // perform convolution in time or frequency domain
    bool _centered; // center the convolution (in the case of time domain)

public:
    convnd1d(){}
    ~convnd1d(){}
    convnd1d(const hypercube<data_t> domain, const std::shared_ptr<vecReg<data_t> > f, bool centered = true){
        successCheck(f->getHyper()->getAxis(1).d==domain.getAxis(1).d,__FILE__,__LINE__,"Data must have the same sampling as the domain fast axis\n");
        _domain=domain;
        _range=*f->getHyper();
        _f=f->clone();
        _centered = centered;
    }
    convnd1d * clone() const {
        convnd1d * op = new convnd1d(_domain,_f,_centered);
        return op;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat);
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat);
};

// transform a filter to zero phase
std::shared_ptr<vecReg<data_t> > zero_phase(const std::shared_ptr<vecReg<data_t> > dat);
// transform a filter to minimum phase
std::shared_ptr<vecReg<data_t> > minimum_phase(const std::shared_ptr<vecReg<data_t> > dat, const data_t eps);
