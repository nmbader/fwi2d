#include "operator.hpp"
#include "fftw3.h"

#define EPS 0.0001

void nloper::dotProduct(){

    std::shared_ptr<vecReg<data_t> > m (new vecReg<data_t>(_domain));
    std::shared_ptr<vecReg<data_t> > dm (new vecReg<data_t>(_domain));
    std::shared_ptr<vecReg<data_t> > d (new vecReg<data_t>(_range));
    std::shared_ptr<vecReg<data_t> > m1;

    m->random(-1,1,123);
    d->random(-1,1,123);
    dm->random(-1,1,123);

    std::shared_ptr<vecReg<data_t> > mtilde = m->clone();
    std::shared_ptr<vecReg<data_t> > d1 = d->clone();
    std::shared_ptr<vecReg<data_t> > dtilde = d->clone();

    data_t eps = 0.1;

    m1 = dm->clone();
    m1->scaleAdd(m,eps,1); // m + eps.dm
    forward(false, m, d1); // f(m)
    forward(false,m1,dtilde); // f(m + eps.dm)
    dtilde->scaleAdd(d1,1,-1); // f(m + eps.dm) - f(m)
    dtilde->scale(1.0/eps); // (df/dm).dm
    jacobianT(false,mtilde,m,d); // (df/dm)'.d

    data_t sum1 = d->dot(dtilde);
    data_t sum2 = dm->dot(mtilde);

    fprintf(stderr,"Dot product with eps = %f:\n",eps);
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));

    eps /= 10;

    m1 = dm->clone();
    m1->scaleAdd(m,eps,1); // m + eps.dm
    forward(false,m1,dtilde); // f(m + eps.dm)
    dtilde->scaleAdd(d1,1,-1); // f(m + eps.dm) - f(m)
    dtilde->scale(1.0/eps); // (df/dm).dm
    jacobianT(false,mtilde,m,d); // (df/dm)'.d

    sum1 = d->dot(dtilde);
    sum2 = dm->dot(mtilde);

    fprintf(stderr,"Dot product with eps = %f:\n",eps);
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));

    eps /= 10;
    
    m1 = dm->clone();
    m1->scaleAdd(m,eps,1); // m + eps.dm
    forward(false,m1,dtilde); // f(m + eps.dm)
    dtilde->scaleAdd(d1,1,-1); // f(m + eps.dm) - f(m)
    dtilde->scale(1.0/eps); // (df/dm).dm
    jacobianT(false,mtilde,m,d); // (df/dm)'.d

    sum1 = d->dot(dtilde);
    sum2 = dm->dot(mtilde);

    fprintf(stderr,"Dot product with eps = %f:\n",eps);
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));   
}

void loper::dotProduct(){

    std::shared_ptr<vecReg<data_t> > m (new vecReg<data_t>(_domain));
    std::shared_ptr<vecReg<data_t> > d (new vecReg<data_t>(_range));

    //time_t timer; unsigned int seed=time(&timer);
    m->random(-1,1,123);
    d->random(-1,1,123);

    std::shared_ptr<vecReg<data_t> > mtild = m->clone();
    std::shared_ptr<vecReg<data_t> > dtild = d->clone();

    this->forward(false, m, dtild);
    this->adjoint(false, mtild, d);

    data_t sum1 = d->dot(dtild);
    data_t sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=false:\n");
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));

    this->forward(true, m, dtild);
    this->adjoint(true, mtild, d);

    sum1 = d->dot(dtild);
    sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=true:\n");
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));
}

void softClip::apply_forward(bool add, const data_t * pmod, data_t * pdat) {
    
    int n = _domain.getN123();
    data_t xbar = 0.5*(_xmax+_xmin);
    data_t dx = 2*(_xmax-xbar)*std::pow(1.0*_q/(_q-1),1.0/_p);
    data_t xmax = xbar + dx;
    data_t xmin = xbar - dx;
    data_t g = 0;

    #pragma omp parallel for private(g)
    for (int i=0; i<n; i++)
    {
        g = std::pow(2*(pmod[i]-xbar)/dx,_p);
        if (std::abs(g) < 1) {
            if (g>0) pdat[i] = add*pdat[i] + xbar + 0.5*dx*std::pow(g-std::pow(g,_q)/_q,1.0/_p);
            else if (g<0) pdat[i] = add*pdat[i] + xbar - 0.5*dx*std::pow(-g+std::pow(g,_q)/_q,1.0/_p);
            else pdat[i] = add*pdat[i] + xbar;
        }
        else if (g>=1) pdat[i] = add*pdat[i] + xbar + 0.5*dx*std::pow(1.0*(_q-1)/_q,1.0/_p);
        else pdat[i] = add*pdat[i] + xbar - 0.5*dx*std::pow(1.0*(_q-1)/_q,1.0/_p);
    }
}
void softClip::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat) {

    int n = _domain.getN123();
    data_t xbar = 0.5*(_xmax+_xmin);
    data_t dx = 2*(_xmax-xbar)*std::pow(1.0*_q/(_q-1),1.0/_p);
    data_t g = 0;
    #pragma omp parallel for private(g)
    for (int i=0; i<n; i++)
    {
        g = std::pow(2*(pmod0[i]-xbar)/dx,_p);
        if (std::abs(g) < 1 && g!=0) {
            if  (g>0) pmod[i] = add*pmod[i] + std::pow(2*(pmod0[i]-xbar)/dx,_p-1)*(1-std::pow(g,_q-1))*std::pow(g-std::pow(g,_q)/_q,1.0*(1-_p)/_p)*pdat[i];
            else pmod[i] = add*pmod[i] + std::pow(2*(pmod0[i]-xbar)/dx,_p-1)*(1-std::pow(g,_q-1))*std::pow(-g+std::pow(g,_q)/_q,1.0*(1-_p)/_p)*pdat[i];
        }
        else pmod[i] = add*pmod[i];
    }
}

void emodelSoftClip::apply_forward(bool add, const data_t * pmod, data_t * pdat){
    
    hypercube<data_t> hyper(_domain.getN123()/_domain.getAxis(_domain.getNdim()).n);
    softClip vpClip(hyper, _vpmin, _vpmax, _p, _q);
    softClip vsClip(hyper, _vsmin, _vsmax, _p, _q);
    softClip rhoClip(hyper, _rhomin, _rhomax, _p, _q);
    softClip spratioClip(hyper, -_spratio, _spratio, _p, _q);

    int n = hyper.getN123();

    data_t (* pm) [n] = (data_t (*) [n]) pmod;
    data_t (* pd) [n] = (data_t (*) [n]) pdat;

    data_t * v0;
    if (add) {
        v0 = new data_t [n];
        memcpy(v0, pd[1], sizeof(data_t)*n);
    }

    // vp, vs, rho clip
    vpClip.apply_forward(add, pm[0], pd[0]);
    vsClip.apply_forward(false, pm[1], pd[1]);
    rhoClip.apply_forward(add, pm[2], pd[2]);

    // vs/vp ratio clip (conditional)
    #pragma omp parallel for
    for (int i=0; i<n; i++) pd[1][i] = pd[1][i] / pd[0][i];

    if (_spratio < 1.0/sqrt(2)) spratioClip.apply_forward(false, pd[1], pd[1]);
    
    #pragma omp parallel for
    for (int i=0; i<n; i++) {
        pd[1][i] = pd[1][i] * pd[0][i];
    }

    if (add){
        #pragma omp parallel for
        for (int i=0; i<n; i++) {
            pd[1][i] += v0[i];
        }
        delete [] v0;
    }

    // copy other parameters if any
    #pragma omp parallel for
    for (int i=3*n; i<_domain.getN123(); i++){
        pdat[i] = add*pdat[i] + pmod[i];
    }
}

void emodelSoftClip::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat){

    hypercube<data_t> hyper(_domain.getN123()/_domain.getAxis(_domain.getNdim()).n);
    softClip vpClip(hyper, _vpmin, _vpmax, _p, _q);
    softClip vsClip(hyper, _vsmin, _vsmax, _p, _q);
    softClip rhoClip(hyper, _rhomin, _rhomax, _p, _q);
    softClip spratioClip(hyper, -_spratio, _spratio, _p, _q);

    int n = hyper.getN123();

    data_t (* pm) [n] = (data_t (*) [n]) pmod;
    const data_t (* pm0) [n] = (const data_t (*) [n]) pmod0;
    const data_t (* pd) [n] = (const data_t (*) [n]) pdat;

    // rho jacobian
    rhoClip.apply_jacobianT(add, pm[2], pm0[2], pd[2]);

    // adjoint for other parameters if any
    #pragma omp parallel for
    for (int i=3*n; i<_domain.getN123(); i++){
        pmod[i] = add*pmod[i] + pdat[i];
    }

    if (_spratio < 1.0/sqrt(2))
    {
        // forward vp,vs clip -> vp,spratio -> spratio clip
        data_t * v0 = new data_t [2*n];
        data_t * v1 = new data_t [n];
        data_t * m0 = new data_t [2*n];
        memset(v0, 0, sizeof(data_t)*2*n);
        memset(m0, 0, sizeof(data_t)*2*n);
        vpClip.apply_forward(false, pmod0, v0);
        vsClip.apply_forward(false, pmod0+n, v0+n);

        #pragma omp parallel for
        for (int i=0; i<n; i++) v0[n+i] = v0[n+i] / v0[i];

        memcpy(v1, v0+n, sizeof(data_t)*n);

        spratioClip.apply_forward(false, v1, v1);

        // jacobian of vp,spratio -> vp,vs
        #pragma omp parallel for
        for (int i=0; i<n; i++) {
            m0[i] = pd[0][i] + v1[i]*pd[1][i];
            m0[n+i] = v0[i]*pd[1][i];
        }

        // jacobian of spratio clip
        spratioClip.apply_jacobianT(false,m0+n, v0+n, m0+n);

        // jacobian of vp,vs -> vp,spratio
        #pragma omp parallel for
        for (int i=0; i<n; i++) {
            m0[i] = m0[i] - v0[n+i]/v0[i]*m0[n+i];
            m0[n+i] /= v0[i];
        }

        // jacobian of vp,vs clip
        vpClip.apply_jacobianT(add, pm[0], pm0[0], m0);
        vsClip.apply_jacobianT(add, pm[1], pm0[1], m0+n);

        delete [] v0;
        delete [] v1;
        delete [] m0;
    }
    else{
        // jacobian of vp,vs
        vpClip.apply_jacobianT(add, pm[0], pm0[0], pd[0]);
        vsClip.apply_jacobianT(add, pm[1], pm0[1], pd[1]);
    }
}

void lam_mu_rho::apply_forward(bool add, const data_t * pmod, data_t * pdat){
    
    int ncomp=1, n=1, nm=1;
    if (_domain.getNdim()==2) {
        ncomp=_domain.getAxis(2).n;
        n = _domain.getAxis(1).n;
    }
    else {
        ncomp=_domain.getAxis(3).n;
        n = _domain.getAxis(1).n*_domain.getAxis(2).n;
        if (_domain.getNdim()==4) nm = _domain.getAxis(4).n;
    }

    const data_t (* pm) [ncomp][n] = (const data_t (*) [ncomp][n]) pmod;
    data_t (*  pd) [ncomp][n] = (data_t (*) [ncomp][n]) pdat;

    for (int im=0; im<nm; im++){
        #pragma omp parallel for
        for (int i=0; i<n; i++){
            pd[im][0][i] = add*pd[im][0][i] + sqrt((pm[im][0][i]+2*pm[im][1][i])/pm[im][2][i]);
            pd[im][1][i] = add*pd[im][1][i] + sqrt(pm[im][1][i]/pm[im][2][i]);
            for (int ic=2; ic<ncomp; ic++) pd[im][ic][i] = add*pd[im][ic][i] + pm[im][ic][i];
        }
    }
}

void lam_mu_rho::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat){

    int ncomp=1, n=1, nm=1;
    if (_domain.getNdim()==2) {
        ncomp=_domain.getAxis(2).n;
        n = _domain.getAxis(1).n;
    }
    else {
        ncomp=_domain.getAxis(3).n;
        n = _domain.getAxis(1).n*_domain.getAxis(2).n;
        if (_domain.getNdim()==4) nm = _domain.getAxis(4).n;
    }

    const data_t (* pm0) [ncomp][n] = (const data_t (*) [ncomp][n]) pmod0;
    const data_t (* pd) [ncomp][n] = (const data_t (*) [ncomp][n]) pdat;
    data_t (*  pm) [ncomp][n] = (data_t (*) [ncomp][n]) pmod;

    for (int im=0; im<nm; im++){
        #pragma omp parallel for
        for (int i=0; i<n; i++){
            pm[im][0][i] = add*pm[im][0][i] + 0.5/sqrt(pm0[im][2][i]*(pm0[im][0][i]+2*pm0[im][1][i]))*pd[im][0][i];
            pm[im][1][i] = add*pm[im][1][i] + 1.0/sqrt(pm0[im][2][i]*(pm0[im][0][i]+2*pm0[im][1][i]))*pd[im][0][i]
                                    + 0.5/sqrt(pm0[im][2][i]*pm0[im][1][i])*pd[im][1][i];
            pm[im][2][i] = add*pm[im][2][i] - 0.5/sqrt(pm0[im][2][i]*pm0[im][2][i]*pm0[im][2][i])
                            * (sqrt(pm0[im][0][i]+2*pm0[im][1][i])*pd[im][0][i] 
                            + sqrt(pm0[im][1][i])*pd[im][1][i]) + pd[im][2][i];
            for (int ic=3; ic<ncomp; ic++) pm[im][ic][i] = add*pm[im][ic][i] + pd[im][ic][i];
        }
    }
}

void lam_mu_rho::apply_inverse(bool add, data_t * pmod, const data_t * pdat){

    int ncomp=1, n=1, nm=1;
    if (_domain.getNdim()==2) {
        ncomp=_domain.getAxis(2).n;
        n = _domain.getAxis(1).n;
    }
    else {
        ncomp=_domain.getAxis(3).n;
        n = _domain.getAxis(1).n*_domain.getAxis(2).n;
        if (_domain.getNdim()==4) nm = _domain.getAxis(4).n;
    }

    const data_t (* pd) [ncomp][n] = (const data_t (*) [ncomp][n]) pdat;
    data_t (*  pm) [ncomp][n] = (data_t (*) [ncomp][n]) pmod;

    for (int im=0; im<nm; im++){
        #pragma omp parallel for
        for (int i=0; i<n; i++){
            pm[im][0][i] = add*pm[im][0][i] + pd[im][2][i]*(pd[im][0][i]*pd[im][0][i] - 2*pd[im][1][i]*pd[im][1][i]);
            pm[im][1][i] = add*pm[im][1][i] + pd[im][2][i]*pd[im][1][i]*pd[im][1][i];
            for (int ic=2; ic<ncomp; ic++) pm[im][ic][i] = add*pm[im][ic][i] + pd[im][ic][i];
        }
    }
}

void ip_is_rho::apply_forward(bool add, const data_t * pmod, data_t * pdat){
    
    int ncomp=1, n=1, nm=1;
    if (_domain.getNdim()==2) {
        ncomp=_domain.getAxis(2).n;
        n = _domain.getAxis(1).n;
    }
    else {
        ncomp=_domain.getAxis(3).n;
        n = _domain.getAxis(1).n*_domain.getAxis(2).n;
        if (_domain.getNdim()==4) nm = _domain.getAxis(4).n;
    }

    const data_t (* pm) [ncomp][n] = (const data_t (*) [ncomp][n]) pmod;
    data_t (*  pd) [ncomp][n] = (data_t (*) [ncomp][n]) pdat;

    for (int im=0; im<nm; im++){
        #pragma omp parallel for
        for (int i=0; i<n; i++){
            pd[im][0][i] = add*pd[im][0][i] + pm[im][0][i]/pm[im][2][i];
            pd[im][1][i] = add*pd[im][1][i] + pm[im][1][i]/pm[im][2][i];
            for (int ic=2; ic<ncomp; ic++) pd[im][ic][i] = add*pd[im][ic][i] + pm[im][ic][i];
        }
    }
}

void ip_is_rho::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat){

    int ncomp=1, n=1, nm=1;
    if (_domain.getNdim()==2) {
        ncomp=_domain.getAxis(2).n;
        n = _domain.getAxis(1).n;
    }
    else {
        ncomp=_domain.getAxis(3).n;
        n = _domain.getAxis(1).n*_domain.getAxis(2).n;
        if (_domain.getNdim()==4) nm = _domain.getAxis(4).n;
    }

    const data_t (* pm0) [ncomp][n] = (const data_t (*) [ncomp][n]) pmod0;
    const data_t (* pd) [ncomp][n] = (const data_t (*) [ncomp][n]) pdat;
    data_t (*  pm) [ncomp][n] = (data_t (*) [ncomp][n]) pmod;

    for (int im=0; im<nm; im++){
        #pragma omp parallel for
        for (int i=0; i<n; i++){
            pm[im][0][i] = add*pm[im][0][i] + pd[im][0][i]/pm0[im][2][i];
            pm[im][1][i] = add*pm[im][1][i] + pd[im][1][i]/pm0[im][2][i];
            pm[im][2][i] = add*pm[im][2][i] - 1.0/(pm0[im][2][i]*pm0[im][2][i])*(pm0[im][0][i]*pd[im][0][i] + pm0[im][1][i]*pd[im][1][i]) + pd[im][2][i];
            for (int ic=3; ic<ncomp; ic++) pm[im][ic][i] = add*pm[im][ic][i] + pd[im][ic][i];
        }
    }
}

void ip_is_rho::apply_inverse(bool add, data_t * pmod, const data_t * pdat){

    int ncomp=1, n=1, nm=1;
    if (_domain.getNdim()==2) {
        ncomp=_domain.getAxis(2).n;
        n = _domain.getAxis(1).n;
    }
    else {
        ncomp=_domain.getAxis(3).n;
        n = _domain.getAxis(1).n*_domain.getAxis(2).n;
        if (_domain.getNdim()==4) nm = _domain.getAxis(4).n;
    }

    const data_t (* pd) [ncomp][n] = (const data_t (*) [ncomp][n]) pdat;
    data_t (*  pm) [ncomp][n] = (data_t (*) [ncomp][n]) pmod;

    for (int im=0; im<nm; im++){
        #pragma omp parallel for
        for (int i=0; i<n; i++){
            pm[im][0][i] = add*pm[im][0][i] + pd[im][0][i]*pd[im][2][i];
            pm[im][1][i] = add*pm[im][1][i] + pd[im][1][i]*pd[im][2][i];
            for (int ic=2; ic<ncomp; ic++) pm[im][ic][i] = add*pm[im][ic][i] + pd[im][ic][i];
        }
    }
}

void vs_vpvs_rho::apply_forward(bool add, const data_t * pmod, data_t * pdat){
    
    int ncomp=1, n=1, nm=1;
    if (_domain.getNdim()==2) {
        ncomp=_domain.getAxis(2).n;
        n = _domain.getAxis(1).n;
    }
    else {
        ncomp=_domain.getAxis(3).n;
        n = _domain.getAxis(1).n*_domain.getAxis(2).n;
        if (_domain.getNdim()==4) nm = _domain.getAxis(4).n;
    }

    const data_t (* pm) [ncomp][n] = (const data_t (*) [ncomp][n]) pmod;
    data_t (*  pd) [ncomp][n] = (data_t (*) [ncomp][n]) pdat;

    for (int im=0; im<nm; im++){
        #pragma omp parallel for
        for (int i=0; i<n; i++){
            pd[im][0][i] = add*pd[im][0][i] + _vs0*exp(pm[im][0][i])*(exp(pm[im][1][i]) + sqrt(2.0));
            pd[im][1][i] = add*pd[im][1][i] + _vs0*exp(pm[im][0][i]);
            pd[im][2][i] = add*pd[im][2][i] + _rho0*exp(pm[im][2][i]);
            for (int ic=3; ic<ncomp; ic++) pd[im][ic][i] = add*pd[im][ic][i] + pm[im][ic][i];
        }
    }
}

void vs_vpvs_rho::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat){

    int ncomp=1, n=1, nm=1;
    if (_domain.getNdim()==2) {
        ncomp=_domain.getAxis(2).n;
        n = _domain.getAxis(1).n;
    }
    else {
        ncomp=_domain.getAxis(3).n;
        n = _domain.getAxis(1).n*_domain.getAxis(2).n;
        if (_domain.getNdim()==4) nm = _domain.getAxis(4).n;
    }

    const data_t (* pm0) [ncomp][n] = (const data_t (*) [ncomp][n]) pmod0;
    const data_t (* pd) [ncomp][n] = (const data_t (*) [ncomp][n]) pdat;
    data_t (*  pm) [ncomp][n] = (data_t (*) [ncomp][n]) pmod;

    for (int im=0; im<nm; im++){
        #pragma omp parallel for
        for (int i=0; i<n; i++){
            pm[im][0][i] = add*pm[im][0][i] + _vs0*exp(pm0[im][0][i])*( (exp(pm0[im][1][i]) + sqrt(2.0))*pd[im][0][i] + pd[im][1][i] );
            pm[im][1][i] = add*pm[im][1][i] + _vs0*exp(pm0[im][0][i]+pm0[im][1][i])*pd[im][0][i];
            pm[im][2][i] = add*pm[im][2][i] + _rho0*exp(pm0[im][2][i])*pd[im][2][i];
            for (int ic=3; ic<ncomp; ic++) pm[im][ic][i] = add*pm[im][ic][i] + pd[im][ic][i];
        }
    }
}

void vs_vpvs_rho::apply_inverse(bool add, data_t * pmod, const data_t * pdat){

    int ncomp=1, n=1, nm=1;
    if (_domain.getNdim()==2) {
        ncomp=_domain.getAxis(2).n;
        n = _domain.getAxis(1).n;
    }
    else {
        ncomp=_domain.getAxis(3).n;
        n = _domain.getAxis(1).n*_domain.getAxis(2).n;
        if (_domain.getNdim()==4) nm = _domain.getAxis(4).n;
    }

    const data_t (* pd) [ncomp][n] = (const data_t (*) [ncomp][n]) pdat;
    data_t (*  pm) [ncomp][n] = (data_t (*) [ncomp][n]) pmod;

    for (int im=0; im<nm; im++){
        #pragma omp parallel for
        for (int i=0; i<n; i++){
            pm[im][0][i] = add*pm[im][0][i] + log(pd[im][1][i]/_vs0);
            pm[im][1][i] = add*pm[im][1][i] + log(std::max( (data_t)EPS, (data_t)(pd[im][0][i]/pd[im][1][i] - sqrt(2.0)) ) );
            pm[im][2][i] = add*pm[im][2][i] + log(pd[im][2][i]/_rho0 );
            for (int ic=3; ic<ncomp; ic++) pm[im][ic][i] = add*pm[im][ic][i] + pd[im][ic][i];
        }
    }
}

void cross_gradient::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    if (!add) memset(pdat,0,_range.getN123()*sizeof(data_t));

    int nz = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    data_t dz = 2*_domain.getAxis(1).d;
    data_t dx = 2*_domain.getAxis(2).d;

    const data_t (* pm)[nx][nz] = (const data_t (*) [nx][nz]) pmod;
    data_t (*  pd)[nz] = (data_t (*) [nz]) pdat;

    #pragma omp parallel for
    for (int ix=1; ix<nx-1; ix++){
        for (int iz=1; iz<nz-1; iz++){
            pd[ix][iz] += ( (pm[0][ix+1][iz] - pm[0][ix-1][iz])*(pm[1][ix][iz+1] - pm[1][ix][iz-1])
                          - (pm[0][ix][iz+1] - pm[0][ix][iz-1])*(pm[1][ix+1][iz] - pm[1][ix-1][iz]) )/(dx*dz);
        }
    }
}

void cross_gradient::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat){

    if (!add) memset(pmod,0,_domain.getN123()*sizeof(data_t));

    int nz = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    data_t dz = 2*_domain.getAxis(1).d;
    data_t dx = 2*_domain.getAxis(2).d;

    const data_t (* pm0)[nx][nz] = (const data_t (*) [nx][nz]) pmod0;
    data_t (* pm)[nx][nz] = (data_t (*) [nx][nz]) pmod;
    const data_t (* pd) [nz] = (const data_t (*) [nz]) pdat;

    for (int ix=1; ix<nx-1; ix++){
        for (int iz=1; iz<nz-1; iz++){
            data_t g1x = (pm0[0][ix+1][iz] - pm0[0][ix-1][iz]) / dx * pd[ix][iz];
            data_t g1z = (pm0[0][ix][iz+1] - pm0[0][ix][iz-1]) / dz * pd[ix][iz];
            data_t g2x = (pm0[1][ix+1][iz] - pm0[1][ix-1][iz]) / dx * pd[ix][iz];
            data_t g2z = (pm0[1][ix][iz+1] - pm0[1][ix][iz-1]) / dz * pd[ix][iz];

            pm[0][ix+1][iz] += g2z / dx;
            pm[0][ix-1][iz] -= g2z / dx;
            pm[0][ix][iz+1] -= g2x / dz;
            pm[0][ix][iz-1] += g2x / dz;

            pm[1][ix+1][iz] -= g1z / dx;
            pm[1][ix-1][iz] += g1z / dx;
            pm[1][ix][iz+1] += g1x / dz;
            pm[1][ix][iz-1] -= g1x / dz;
        }
    }
}

void polar::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    #pragma omp parallel for
    for (int i=0; i<_domain.getN123()/2; i++){
        data_t val = pmod[2*i];
        pdat[2*i]= sqrt(pmod[2*i]*pmod[2*i] + pmod[2*i+1]*pmod[2*i+1]);
        pdat[2*i+1]= atan2(pmod[2*i+1], val);
    }
}

void polar::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat){

    #pragma omp parallel for
    for (int i=0; i<_domain.getN123()/2; i++) {
        data_t r = sqrt(pmod0[2*i]*pmod0[2*i] + pmod0[2*i+1]*pmod0[2*i+1]);
        data_t x = pmod0[2*i];
        data_t y = pmod0[2*i+1];
        data_t val = pdat[2*i];
        pmod[2*i]= x/r*pdat[2*i] - y/r*pdat[2*i+1];
        pmod[2*i+1]= y/r*val + x/r*pdat[2*i+1];
    }
}

void ipolar::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    #pragma omp parallel for
    for (int i=0; i<_domain.getN123()/2; i++){
        data_t val = pmod[2*i];
        pdat[2*i] = pmod[2*i]*cos(pmod[2*i+1]);
        pdat[2*i+1] = val*sin(pmod[2*i+1]);
    }
}
void ipolar::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat){

    #pragma omp parallel for
    for (int i=0; i<_domain.getN123()/2; i++) {
        data_t r = pmod0[2*i];
        data_t phi = pmod0[2*i+1];
        data_t val = pdat[2*i];
        pmod[2*i] = cos(phi)*pdat[2*i] + sin(phi)*pdat[2*i+1];
        pmod[2*i+1] = - r*sin(phi)*val + r*cos(phi)*pdat[2*i+1];
    }
}

void sdivision::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    if (!add) memset(pdat,0,_domain.getN123()*sizeof(data_t));

    int nf=_domain.getAxis(1).n/2;
    int ntr=_domain.getAxis(2).n;
    int ny = _domain.getN123()/(2*nf*ntr);
    int fmin=round(_fmin*nf);
    int fmax=round(_fmax*nf);

    const data_t (*pv1) [ntr][2*nf] = (const data_t (*)[ntr][2*nf]) pmod;
    data_t (*pv2) [ntr][2*nf] = (data_t (*)[ntr][2*nf]) pdat;

    for (int iy=0; iy<ny; iy++){
        #pragma omp parallel for
        for (int itr=1; itr<ntr; itr++){
            for (int f=fmin; f<fmax; f++){
                pv2[iy][itr][2*f] += (pv1[iy][itr][2*f]+_eps) / (pv1[iy][itr-1][2*f]+_eps);
                pv2[iy][itr][2*f+1] += pv1[iy][itr][2*f+1] - pv1[iy][itr-1][2*f+1];
            }
        }
    }
}
void sdivision::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat){

    if (!add) memset(pmod,0,_domain.getN123()*sizeof(data_t));

    int nf=_domain.getAxis(1).n/2;
    int ntr=_domain.getAxis(2).n;
    int ny = _domain.getN123()/(2*nf*ntr);
    int fmin=round(_fmin*nf);
    int fmax=round(_fmax*nf);

    data_t (*pv2) [ntr][2*nf] = (data_t (*)[ntr][2*nf]) pmod;
    const data_t (*pv0) [ntr][2*nf] = (const data_t (*)[ntr][2*nf]) pmod0;
    const data_t (*pv1) [ntr][2*nf] = (const data_t (*)[ntr][2*nf]) pdat;

    for (int iy=0; iy<ny; iy++){
        #pragma omp parallel for
        for (int itr=1; itr<ntr-1; itr++){
            for (int f=fmin; f<fmax; f++){
                pv2[iy][itr][2*f] += pv1[iy][itr][2*f]/(pv0[iy][itr-1][2*f]+_eps) - pv1[iy][itr+1][2*f]*(pv0[iy][itr+1][2*f]+_eps)/((pv0[iy][itr][2*f]+_eps)*(pv0[iy][itr][2*f]+_eps));
                pv2[iy][itr][2*f+1] += pv1[iy][itr][2*f+1] - pv1[iy][itr+1][2*f+1];
            }
        }
        for (int f=fmin; f<fmax; f++){
            pv2[iy][0][2*f] += -pv1[iy][1][2*f]*(pv0[iy][1][2*f]+_eps)/((pv0[iy][0][2*f]+_eps)*(pv0[iy][0][2*f]+_eps));
            pv2[iy][0][2*f+1] += -pv1[iy][1][2*f+1];
            pv2[iy][ntr-1][2*f] += pv1[iy][ntr-1][2*f]/(pv0[iy][ntr-2][2*f]+_eps);
            pv2[iy][ntr-1][2*f+1] += pv1[iy][ntr-1][2*f+1];
        }
    }
}

void tr2trDecon::apply_forward(bool add, const data_t * pmod, data_t * pdat) {

    // This operator can be written as: d = F^-1.P^-1(S.D(P(F.m))))
    // F is FT transform (real to real vectors)
    // P is polar coordinates transform (non-linear)
    // D is spectral division in a given frequency band between traces with added white noise: r_i,r_i+1,r_i+2,phi_i,phi_i+1,phi_i+2 --> 0,(r_i+1+e)/(r_i+e),(r_i+2+e)/(r_i+1+e),0,phi_i+1-phi_i,phi_i+2-phi_i+1
    // S is a triangular smoothing applied to r
    // F^-1 is inverse FT transform (real to real vectors)

    // F.m
    fxTransformR fx(_domain);
    std::shared_ptr<vecReg<data_t> > vec1 = std::make_shared<vecReg<data_t> >(*fx.getRange());
    vec1->zero();
    data_t * pvec1=vec1->getVals();
    fx.apply_forward(false,pmod,pvec1);

    // P(F.m)
    polar P(*vec1->getHyper());
    P.apply_forward(false,pvec1,pvec1);

    // D(P(F.m))
    sdivision D(*vec1->getHyper(),_fmin,_fmax,_eps);
    std::shared_ptr<vecReg<data_t> > vec2 = std::make_shared<vecReg<data_t> >(*D.getRange());
    vec2->zero();
    data_t * pvec2=vec2->getVals();
    D.apply_forward(false, pvec1, pvec2);

    // S.D(P(F.m))
    ssmth S(*vec2->getHyper(),_hl);
    S.apply_forward(false, pvec2, pvec1);

    // P^-1(S.D(P(F.m)))
    ipolar iP(*vec1->getHyper());
    iP.apply_forward(false,pvec1,pvec1);

    // F^-1.P^-1(S.D(P(F.m)))
    ifxTransformR ifx(_range);
    ifx.apply_forward(add,pvec1,pdat);
}

void tr2trDecon::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat){

    // d = F^-1.P^-1(S.D(P(F.m))))
    // The adjoint of the above operator is: m = F'.(dP/dF)'.(dD/dP)'.S'.(dP^-1/dS)'.F'^-1.d
    // dP/dF is evaluated at F_0=(F.m_0)
    // dD/dP is evaluated at P_0=P(F.m_0)
    // dP^-1/dS is evaluated at S_0=S.D((P(F.m_0))

    // F'^-1.d
    ifxTransformR ifx(_range);
    std::shared_ptr<vecReg<data_t> > vec1 = std::make_shared<vecReg<data_t> >(*ifx.getDomain());
    vec1->zero();
    data_t * pvec1=vec1->getVals();
    ifx.apply_adjoint(false,pvec1,pdat);

    // S_0
    std::shared_ptr<vecReg<data_t> > vec0 = std::make_shared<vecReg<data_t> >(_range);
    vec0->zero();
    apply_forward(false, pmod0, vec0->getVals());
    fxTransformR fx(_domain);
    std::shared_ptr<vecReg<data_t> > vec2 = std::make_shared<vecReg<data_t> >(*fx.getRange());
    vec2->zero();
    data_t * pvec2=vec2->getVals();
    fx.apply_forward(false,vec0->getVals(),pvec2);
    polar P(*vec2->getHyper());
    P.apply_forward(false,pvec2,pvec2);

    // (dP^-1/dS)'.F'^-1.d
    ipolar iP(*vec1->getHyper());
    iP.apply_jacobianT(false,pvec1,pvec2,pvec1);

    // S'.(dP^-1/dS)'.F'^-1.d
    ssmth S(*vec1->getHyper(),_hl);
    S.apply_adjoint(false, pvec2, pvec1);

    // P_0
    fx.apply_forward(false,pmod0,pvec1);
    P.apply_forward(false,pvec1,pvec1);

    // (dD/dP)'.S'.(dP^-1/dS)'.F'^-1.d
    sdivision D(*vec2->getHyper(),_fmin,_fmax,_eps);
    std::shared_ptr<vecReg<data_t> > vec3 = std::make_shared<vecReg<data_t> >(*D.getDomain());
    vec3->zero();
    data_t * pvec3=vec3->getVals();
    D.apply_jacobianT(false,pvec3,pvec1,pvec2);

    // F_0
    iP.apply_forward(false,pvec1,pvec1);

    // (dP/dF)'.(dD/dP)'.S'.(dP^-1/dS)'.F'^-1.d
    P.apply_jacobianT(false,pvec3,pvec1,pvec3);

    // F'.(dP/dF)'.(dD/dP)'.S'.(dP^-1/dS)'.F'^-1.d
    fx.apply_adjoint(add, pmod, pvec3);
}

inline data_t sinc(data_t x){
    if (x==0.0) return 1.0;
    else return sin(x)/x;
}

void linear_resampler::apply_forward(bool add, const data_t * pmod, data_t * pdat) {

    axis<data_t> Tm = _domain.getAxis(1);
    axis<data_t> Td = _range.getAxis(1);
    int nx = _domain.getN123()/Tm.n;

    int i, j;
    data_t wi, wj;
    #pragma omp parallel for private(i,j,wi,wj)
    for (int ix=0; ix<nx; ix++){
        for (int it=0; it<Td.n; it++){

            i = std::min(Tm.n - 1, (int) floor(it*Td.d / Tm.d));
            j = std::min(Tm.n - 1, (int) ceil(it*Td.d / Tm.d));
            wj = (it*Td.d - i*Tm.d)/Tm.d;
            wi = 1 - wj;
            pdat[ix*Td.n+it] = add*pdat[ix*Td.n+it] + wi*pmod[ix*Tm.n+i] + wj*pmod[ix*Tm.n+j];
        }
    }

}

void linear_resampler::apply_adjoint(bool add, data_t * pmod, const data_t * pdat) {

    axis<data_t> Tm = _domain.getAxis(1);
    axis<data_t> Td = _range.getAxis(1);
    int nx = _domain.getN123()/Tm.n;

    if (add == false) {
        for (int i=0; i<Tm.n*nx; i++) pmod[i] = 0;
    }

    int i, j;
    data_t wi, wj;
    #pragma omp parallel for private(i,j,wi,wj)
    for (int ix=0; ix<nx; ix++){
        for (int it=0; it<Td.n; it++){

            i = std::min(Tm.n - 1, (int) floor(it*Td.d / Tm.d));
            j = std::min(Tm.n - 1, (int) ceil(it*Td.d / Tm.d));
            wj = (it*Td.d - i*Tm.d)/Tm.d;
            wi = 1 - wj;
            pmod[ix*Tm.n+i] += wi*pdat[ix*Td.n+it];
            pmod[ix*Tm.n+j] += wj*pdat[ix*Td.n+it];
        }
    }
}

void sinc_resampler::apply_forward(bool add, const data_t * pmod, data_t * pdat) {

    axis<data_t> Tm = _domain.getAxis(1);
    axis<data_t> Td = _range.getAxis(1);
    int nx = _domain.getN123()/Tm.n;

    int itm0;
    data_t t;
    data_t a = _alpha; // 0 < a < 1 for the windowing of the cosine taper

    #pragma omp parallel for private(itm0,t)
    for (int ix=0; ix<nx; ix++){
        for (int itd=0; itd<Td.n; itd++){
            t = itd * Td.d;
            itm0 = floor(t / Tm.d);
            pdat[ix*Td.n+itd] = add*pdat[ix*Td.n+itd];
            for (int itm = std::max(0, itm0 -_hl + 1) ; itm <= std::min(Tm.n-1, itm0 + _hl); itm++){
                pdat[ix*Td.n+itd] += pmod[ix*Tm.n+itm]*sinc(M_PI * (t-itm*Tm.d) / Tm.d)*(a +(1-a)*cos(M_PI*(itm-t/Tm.d)/_hl));
            }
        }
    }

}

void sinc_resampler::apply_adjoint(bool add, data_t * pmod, const data_t * pdat) {

    axis<data_t> Tm = _domain.getAxis(1);
    axis<data_t> Td = _range.getAxis(1);
    int nx = _domain.getN123()/Tm.n;

    if (add == false) {
       for (int i=0; i<Tm.n*nx; i++) pmod[i] = 0;
    }

    int itm0;
    data_t t;
    data_t a = _alpha; // 0 < a < 1 for the windowing of the cosine taper
    
    #pragma omp parallel for private(itm0,t)
    for (int ix=0; ix<nx; ix++){
        for (int itd=0; itd<Td.n; itd++){
            t = itd * Td.d;
            itm0 = floor(t / Tm.d);
            for (int itm = std::max(0, itm0 -_hl + 1) ; itm <= std::min(Tm.n-1, itm0 + _hl); itm++){
                pmod[ix*Tm.n+itm] += pdat[ix*Td.n+itd]*sinc(M_PI * (t-itm*Tm.d) / Tm.d)*(a +(1-a)*cos(M_PI*(itm-t/Tm.d)/_hl));
            }
        }
    }
}

void fxTransform::apply_forward(bool add, const data_t * pmod, std::complex<data_t> * pdat)
{
    int nt = _domain.getAxis(1).n;
    int ntr = _domain.getN123()/nt;
    int nf = _range.getAxis(1).n;
    const data_t (*m) [nt] = (const data_t (*) [nt]) pmod;
    std::complex<data_t> (*d) [nf] = (std::complex<data_t> (*) [nf]) pdat;

    // double precision
    if ((int)sizeof(data_t)==8){

        double * in = new double[nt];
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf);
        fftw_plan p = fftw_plan_dft_r2c_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt);

            // forward transform
            fftw_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) {
                d[itr][f].real(add*d[itr][f].real() + out[f][0]);
                d[itr][f].imag(add*d[itr][f].imag() + out[f][1]);
            }
        }

        delete [] in;
        fftw_destroy_plan(p);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        float * in = new float[nt];
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf);
        fftwf_plan p = fftwf_plan_dft_r2c_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt);

            // forward transform
            fftwf_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) {
                d[itr][f].real(add*d[itr][f].real() + out[f][0]);
                d[itr][f].imag(add*d[itr][f].imag() + out[f][1]);
            }
        }

        delete [] in;
        fftwf_destroy_plan(p);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void fxTransform::apply_inverse(bool add, data_t * pmod, const std::complex<data_t> * pdat)
{

    int nt = _domain.getAxis(1).n;
    int ntr = _domain.getN123()/nt;
    int nf = _range.getAxis(1).n;
    data_t (*m) [nt] = (data_t (*) [nt]) pmod;
    const std::complex<data_t> (*d) [nf] = (const std::complex<data_t> (*) [nf]) pdat;

    // double precision
    if ((int)sizeof(data_t)==8){

        double * out = new double[nt];
        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf);
        fftw_plan p = fftw_plan_dft_c2r_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, d[itr], sizeof(fftw_complex)*nf);

            // inverse transform
            fftw_execute(p);

            // copy to output
            for (int it=0; it<nt; it++) m[itr][it]= add*m[itr][it] + out[it]/nt;
        }

        delete [] out;
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_cleanup();
    }
    else
    {
        float * out = new float[nt];
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf);
        fftwf_plan p = fftwf_plan_dft_c2r_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, d[itr], sizeof(fftwf_complex)*nf);

            // inverse transform
            fftwf_execute(p);

            // copy to output
            for (int it=0; it<nt; it++) m[itr][it]= add*m[itr][it] + out[it]/nt;
        }

        delete [] out;
        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_cleanup();
    }
}

void fxTransform::cforward(bool add, const std::shared_ptr<cvecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkCDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int ntr = _domain.getN123()/nt;
    int nf = _crange.getAxis(1).n;
    const std::complex<data_t> (*m) [nt] = (const std::complex<data_t> (*) [nt]) mod->getCVals();
    std::complex<data_t> (*d) [nf] = (std::complex<data_t> (*) [nf]) dat->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){

        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt);
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt);
        fftw_plan p = fftw_plan_dft_1d(nt, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(fftw_complex)*nt);

            // forward transform
            fftw_execute(p);

            // copy to output
            for (int f=0; f<(nf-1)/2; f++) {
                d[itr][f].real(add*d[itr][f].real() + out[nf-(nf-1)/2+f][0]);
                d[itr][f].imag(add*d[itr][f].imag() + out[nf-(nf-1)/2+f][1]);
            }
            for (int f=(nf-1)/2; f<nf; f++) {
                d[itr][f].real(add*d[itr][f].real() + out[f-(nf-1)/2][0]);
                d[itr][f].imag(add*d[itr][f].imag() + out[f-(nf-1)/2][1]);
            }
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
        fftwf_plan p = fftwf_plan_dft_1d(nt, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(fftwf_complex)*nt);

            // forward transform
            fftwf_execute(p);

            // copy to output
            for (int f=0; f<(nf-1)/2; f++) {
                d[itr][f].real(add*d[itr][f].real() + out[nf-(nf-1)/2+f][0]);
                d[itr][f].imag(add*d[itr][f].imag() + out[nf-(nf-1)/2+f][1]);
            }
            for (int f=(nf-1)/2; f<nf; f++) {
                d[itr][f].real(add*d[itr][f].real() + out[f-(nf-1)/2][0]);
                d[itr][f].imag(add*d[itr][f].imag() + out[f-(nf-1)/2][1]);
            }
        }

        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void fxTransform::cinverse(bool add, std::shared_ptr<cvecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkCDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int ntr = _domain.getN123()/nt;
    int nf = _crange.getAxis(1).n;
    std::complex<data_t> (*m) [nt] = (std::complex<data_t> (*) [nt]) mod->getCVals();
    const std::complex<data_t> (*d) [nf] = (const std::complex<data_t> (*) [nf]) dat->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){

        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt);
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt);
        fftw_plan p = fftw_plan_dft_1d(nt, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in+nt-(nt-1)/2, d[itr], sizeof(fftw_complex)*(nt-1)/2);
            memcpy(in, d[itr]+(nt-1)/2, sizeof(fftw_complex)*(nt-(nt-1)/2));

            // forward transform
            fftw_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) {
                m[itr][f].real(add*m[itr][f].real() + out[f][0]/nt);
                m[itr][f].imag(add*m[itr][f].imag() + out[f][1]/nt);
            }
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt);
        fftwf_plan p = fftwf_plan_dft_1d(nt, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in+nt-(nt-1)/2, d[itr], sizeof(fftwf_complex)*(nt-1)/2);
            memcpy(in, d[itr]+(nt-1)/2, sizeof(fftwf_complex)*(nt-(nt-1)/2));

            // forward transform
            fftwf_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) {
                m[itr][f].real(add*m[itr][f].real() + out[f][0]/nt);
                m[itr][f].imag(add*m[itr][f].imag() + out[f][1]/nt);
            }
        }

        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void fxTransform::testInverse() {

    std::shared_ptr<vecReg<data_t> > mod1 (new vecReg<data_t> (_domain));
    std::shared_ptr<vecReg<data_t> > mod1b (new vecReg<data_t> (_domain));

    mod1->random();

    std::shared_ptr<cvecReg<data_t> > dat1 (new cvecReg<data_t> (_range));
    dat1->zero();
    mod1b->zero();

    this->forward(false,mod1,dat1);
    this->inverse(false,mod1b,dat1);

    data_t trueNorm1, invNorm1, error1;
    trueNorm1 = mod1->norm2();
    invNorm1 = mod1b->norm2();
    mod1b->scaleAdd(mod1,1,-1);
    error1 = mod1b->norm2();

    fprintf(stderr,"input = %f output = %f inverse = %f relative error: %f\n", trueNorm1, dat1->norm2(), invNorm1, error1/trueNorm1);
}

void fxTransform::testCInverse() {

    std::shared_ptr<cvecReg<data_t> > mod1 (new cvecReg<data_t> (_domain));
    std::shared_ptr<cvecReg<data_t> > mod1b (new cvecReg<data_t> (_domain));

    mod1->random();
    
    std::shared_ptr<cvecReg<data_t> > dat1 (new cvecReg<data_t> (_crange));
    dat1->zero();
    mod1b->zero();

    this->cforward(false,mod1,dat1);
    this->cinverse(false,mod1b,dat1);
    
    data_t trueNorm1, invNorm1, error1;
    trueNorm1 = mod1->norm2();
    invNorm1 = mod1b->norm2();
    mod1b->scaleAdd(mod1,1,-1);
    error1 = mod1b->norm2();

    fprintf(stderr,"input = %f output = %f inverse = %f relative error: %f\n", trueNorm1, dat1->norm2(), invNorm1, error1/trueNorm1);
}

void fxTransform::testAdjoint(){

    std::shared_ptr<vecReg<data_t> > m (new vecReg<data_t>(_domain));
    std::shared_ptr<cvecReg<data_t> > d (new cvecReg<data_t>(_range));

    //time_t timer; unsigned int seed=time(&timer);
    m->random();
    d->random();
    
    int nf = _range.getAxis(1).n;
    int nx = _range.getN123()/nf;
    for (int i=0; i<nx; i++) d->getVals()[i*nf].imag(0);

    std::shared_ptr<vecReg<data_t> > mtild = m->clone();
    std::shared_ptr<cvecReg<data_t> > dtild = d->clone();

    this->forward(false, m, dtild);
    this->adjoint(false, mtild, d);

    data_t sum1 = dtild->real()->dot(d->real()) + dtild->imag()->dot(d->imag());
    data_t sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=false:\n");
    fprintf(stderr,"sum1 = %f\nsum2 = %f\n",sum1, sum2);

    this->forward(true, m, dtild);
    this->adjoint(true, mtild, d);

    sum1 = dtild->real()->dot(d->real()) + dtild->imag()->dot(d->imag());
    sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=true:\n");
    fprintf(stderr,"sum1 = %f\nsum2 = %f\n",sum1, sum2);
}

void fxTransform::testCAdjoint(){

    std::shared_ptr<cvecReg<data_t> > m (new cvecReg<data_t>(_domain));
    std::shared_ptr<cvecReg<data_t> > d (new cvecReg<data_t>(_crange));

    //time_t timer; unsigned int seed=time(&timer);
    m->random();
    d->random();

    std::shared_ptr<cvecReg<data_t> > mtild = m->clone();
    std::shared_ptr<cvecReg<data_t> > dtild = d->clone();

    this->cforward(false, m, dtild);
    dtild->scale(1.0/sqrt(_domain.getAxis(1).n));
    this->cadjoint(false, mtild, d);
    mtild->scale(sqrt(_domain.getAxis(1).n));

    std::complex<data_t> sum1 = dtild->dot(d);
    std::complex<data_t> sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=false:\n");
    fprintf(stderr,"sum1 = %f + i%f\nsum2 = %f + i%f\n",sum1.real(), sum1.imag(), sum2.real(), sum2.imag());

    m->scale(1.0/sqrt(_domain.getAxis(1).n));
    this->cforward(true, m, dtild);
    m->scale(sqrt(_domain.getAxis(1).n));
    d->scale(sqrt(_domain.getAxis(1).n));
    this->cadjoint(true, mtild, d);
    d->scale(1.0/sqrt(_domain.getAxis(1).n));

    sum1 = dtild->dot(d);
    sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=true:\n");
    fprintf(stderr,"sum1 = %f + i%f\nsum2 = %f + i%f\n",sum1.real(), sum1.imag(), sum2.real(), sum2.imag());
}

void fkTransform::forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ntr = _domain.getN123()/(nt*nx);
    int nf = _range.getAxis(1).n;
    const data_t (*m) [nx][nt] = (const data_t (*) [nx][nt]) mod->getCVals();
    std::complex<data_t> (*d) [nx][nf] = (std::complex<data_t> (*) [nx][nf]) dat->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){

        double * in = new double[nt*nx];
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*nx);
        fftw_plan p = fftw_plan_dft_r2c_2d(nx, nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt*nx);

            // forward transform
            fftw_execute(p);

            // copy to output
            int x=0, f=0;
            for (int i=0; i<nx*nf; i++){
                x = i/nf;
                f = i - x*nf;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                d[itr][nx-x-1][f].real(add*d[itr][nx-x-1][f].real() + out[i][0]);
                d[itr][nx-x-1][f].imag(add*d[itr][nx-x-1][f].imag() + out[i][1]);
            }
        }

        delete [] in;
        fftw_destroy_plan(p);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        float * in = new float[nt*nx];
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf*nx);
        fftwf_plan p = fftwf_plan_dft_r2c_2d(nx, nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt*nx);

            // forward transform
            fftwf_execute(p);

            // copy to output
            int x=0, f=0;
            for (int i=0; i<nx*nf; i++){
                x = i/nf;
                f = i - x*nf;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                d[itr][nx-x-1][f].real(add*d[itr][nx-x-1][f].real() + out[i][0]);
                d[itr][nx-x-1][f].imag(add*d[itr][nx-x-1][f].imag() + out[i][1]);
            }
        }

        delete [] in;
        fftwf_destroy_plan(p);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void fkTransform::inverse(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ntr = _domain.getN123()/(nt*nx);
    int nf = _range.getAxis(1).n;
    data_t (*m) [nx][nt] = (data_t (*) [nx][nt]) mod->getCVals();
    const std::complex<data_t> (*d) [nx][nf] = (const std::complex<data_t> (*) [nx][nf]) dat->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){

        double * out = new double[nt*nx];
        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*nx);
        fftw_plan p = fftw_plan_dft_c2r_2d(nx, nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            int x=0, f=0;
            for (int i=0; i<nx*nf; i++){
                x = i/nf;
                f = i - x*nf;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                in[i][0] = d[itr][nx-1-x][f].real();
                in[i][1] = d[itr][nx-1-x][f].imag();
            }

            // inverse transform
            fftw_execute(p);

            // copy to output
            int it=0;
            for (int i=0; i<nx*nt; i++){
                x = i/nt;
                it = i - x*nt;
                m[itr][x][it]= add*m[itr][x][it] + out[i]/(nt*nx);
            }
        }

        delete [] out;
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_cleanup();
    }
    else
    {
        float * out = new float[nt*nx];
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf*nx);
        fftwf_plan p = fftwf_plan_dft_c2r_2d(nx, nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            int x=0, f=0;
            for (int i=0; i<nx*nf; i++){
                x = i/nf;
                f = i - x*nf;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                in[i][0] = d[itr][nx-1-x][f].real();
                in[i][1] = d[itr][nx-1-x][f].imag();
            }

            // inverse transform
            fftwf_execute(p);

            // copy to output
            int it=0;
            for (int i=0; i<nx*nt; i++){
                x = i/nt;
                it = i - x*nt;
                m[itr][x][it]= add*m[itr][x][it] + out[i]/(nt*nx);
            }
        }

        delete [] out;
        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_cleanup();
    }
}

void fkTransform::cforward(bool add, const std::shared_ptr<cvecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkCDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ntr = _domain.getN123()/(nt*nx);
    int nf = _crange.getAxis(1).n;
    const data_t (*m) [nx][nt] = (const data_t (*) [nx][nt]) mod->getCVals();
    std::complex<data_t> (*d) [nx][nf] = (std::complex<data_t> (*) [nx][nf]) dat->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){

        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*nx);
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*nx);
        fftw_plan p = fftw_plan_dft_2d(nx, nt, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt*nx);

            // forward transform
            fftw_execute(p);

            // copy to output
            int x=0, f=0;
            for (int i=0; i<nx*nf; i++){
                x = i/nt;
                f = i - x*nt;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                f += (nt-1)/2;
                f = f - nt*(f/nt);
                d[itr][x][f].real(add*d[itr][x][f].real() + out[i][0]);
                d[itr][x][f].imag(add*d[itr][x][f].imag() + out[i][1]);
            }
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf*nx);
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf*nx);
        fftwf_plan p = fftwf_plan_dft_2d(nx, nt, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt*nx);

            // forward transform
            fftwf_execute(p);

            // copy to output
            int x=0, f=0;
            for (int i=0; i<nx*nf; i++){
                x = i/nt;
                f = i - x*nt;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                f += (nt-1)/2;
                f = f - nt*(f/nt);
                d[itr][x][f].real(add*d[itr][x][f].real() + out[i][0]);
                d[itr][x][f].imag(add*d[itr][x][f].imag() + out[i][1]);
            }
        }

        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void fkTransform::cinverse(bool add, std::shared_ptr<cvecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkCDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ntr = _domain.getN123()/(nt*nx);
    int nf = _crange.getAxis(1).n;
    std::complex<data_t> (*m) [nx][nt] = (std::complex<data_t> (*) [nx][nt]) mod->getCVals();
    const std::complex<data_t> (*d) [nx][nf] = (const std::complex<data_t> (*) [nx][nf]) dat->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){

        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt*nx);
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt*nx);
        fftw_plan p = fftw_plan_dft_2d(nx, nt, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            int x=0, f=0;
            for (int i=0; i<nx*nt; i++){
                x = i/nt;
                f = i - x*nt;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                f += (nt-1)/2;
                f = f - nt*(f/nt);
                in[i][0] = d[itr][x][f].real();
                in[i][1] = d[itr][x][f].imag();
            }

            // forward transform
            fftw_execute(p);

            // copy to output
            int it=0;
            for (int i=0; i<nx*nt; i++){
                x = i/nt;
                it = i - x*nt;
                m[itr][x][f].real(add*m[itr][x][f].real() + out[i][0]/(nt*nx));
                m[itr][x][f].imag(add*m[itr][x][f].imag() + out[i][1]/(nt*nx));
            }
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nt*nx);
        fftwf_plan p = fftwf_plan_dft_2d(nx, nt, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            int x=0, f=0;
            for (int i=0; i<nx*nt; i++){
                x = i/nt;
                f = i - x*nt;
                x += (nx-1)/2;
                x = x - nx*(x/nx);
                f += (nt-1)/2;
                f = f - nt*(f/nt);
                in[i][0] = d[itr][x][f].real();
                in[i][1] = d[itr][x][f].imag();
            }

            // forward transform
            fftwf_execute(p);

            // copy to output
            int it=0;
            for (int i=0; i<nx*nt; i++){
                x = i/nt;
                it = i - x*nt;
                m[itr][x][f].real(add*m[itr][x][f].real() + out[i][0]/(nt*nx));
                m[itr][x][f].imag(add*m[itr][x][f].imag() + out[i][1]/(nt*nx));
            }
        }

        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void fxTransformR::apply_forward(bool add, const data_t * pmod, data_t * pdat)
{   
    int nt = _domain.getAxis(1).n;
    int ntr = _domain.getN123()/nt;
    int nf = _range.getAxis(1).n/2;
    const data_t (*m) [nt] = (const data_t (*) [nt]) pmod;
    data_t (*d) [2*nf] = (data_t (*) [2*nf]) pdat;

    // double precision
    if ((int)sizeof(data_t)==8){

        double * in = new double[nt];
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf);
        fftw_plan p = fftw_plan_dft_r2c_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt);

            // forward transform
            fftw_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) {
                d[itr][2*f] = add*d[itr][2*f] + out[f][0]/sqrt(nt); // real part
                d[itr][2*f+1] = add*d[itr][2*f+1] + out[f][1]/sqrt(nt); // imaginary part
            }
        }
        delete [] in;
        fftw_destroy_plan(p);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        float * in = new float[nt];
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf);
        fftwf_plan p = fftwf_plan_dft_r2c_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt);

            // forward transform
            fftwf_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) {
                d[itr][2*f] = add*d[itr][2*f] + out[f][0]/sqrt(nt); // real part
                d[itr][2*f+1] = add*d[itr][2*f+1] + out[f][1]/sqrt(nt); // imaginary part
            }
        }
        delete [] in;
        fftwf_destroy_plan(p);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void fxTransformR::apply_adjoint(bool add, data_t * pmod, const data_t * pdat)
{
    int nt = _domain.getAxis(1).n;
    int ntr = _domain.getN123()/nt;
    int nf = _range.getAxis(1).n/2;
    data_t (*m) [nt] = (data_t (*) [nt]) pmod;
    const data_t (*d) [2*nf] = (const data_t (*) [2*nf]) pdat;

    if (!add) memset(pmod, 0, sizeof(data_t)*_domain.getN123());

    // double precision
    if ((int)sizeof(data_t)==8){

        double * in = new double[nt];
        memset(in, 0, nt*sizeof(double));
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf);
        fftw_plan p = fftw_plan_dft_r2c_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer real input
            for (int f=0; f<nf; f++) in[f] = d[itr][2*f];

            // forward transform
            fftw_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) m[itr][f] +=  out[f][0]/sqrt(nt);
            for (int f=nf; f<nt; f++) m[itr][f] += out[2*nf-f-1][0]/sqrt(nt);

            // transfer imaginary input
            for (int f=0; f<nf; f++) in[f] = d[itr][2*f+1];

            // forward transform
            fftw_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) m[itr][f] +=  out[f][1]/sqrt(nt);
            for (int f=nf; f<nt; f++) m[itr][f] -= out[2*nf-f-1][1]/sqrt(nt);
        }
        delete [] in;
        fftw_destroy_plan(p);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        float * in = new float[nt];
        memset(in, 0, nt*sizeof(float));
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf);
        fftwf_plan p = fftwf_plan_dft_r2c_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer real input
            for (int f=0; f<nf; f++) in[f] = d[itr][2*f];

            // forward transform
            fftwf_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) m[itr][f] += out[f][0]/sqrt(nt);
            for (int f=nf; f<nt; f++) m[itr][f] += out[2*nf-f-1][0]/sqrt(nt);

            // transfer imaginary input
            for (int f=0; f<nf; f++) in[f] = d[itr][2*f+1];

            // forward transform
            fftwf_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) m[itr][f] += out[f][1]/sqrt(nt);
            for (int f=nf; f<nt; f++) m[itr][f] -= out[2*nf-f-1][1]/sqrt(nt);
        }
        delete [] in;
        fftwf_destroy_plan(p);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void ifxTransformR::apply_forward(bool add, const data_t * pmod, data_t * pdat)
{
    int nt = _range.getAxis(1).n;
    int ntr = _range.getN123()/nt;
    int nf = _domain.getAxis(1).n/2;
    data_t (*m) [nt] = (data_t (*) [nt]) pdat;
    const data_t (*d) [2*nf] = (const data_t (*) [2*nf]) pmod;

    // double precision
    if ((int)sizeof(data_t)==8){

        double * out = new double[nt];
        fftw_complex * in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf);
        fftw_plan p = fftw_plan_dft_c2r_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, d[itr], sizeof(fftw_complex)*nf);

            // inverse transform
            fftw_execute(p);

            // copy to output
            for (int it=0; it<nt; it++) m[itr][it]= add*m[itr][it] + out[it]/sqrt(nt);
        }
        delete [] out;
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_cleanup();
    }
    else
    {
        float * out = new float[nt];
        fftwf_complex * in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf);
        fftwf_plan p = fftwf_plan_dft_c2r_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, d[itr], sizeof(fftwf_complex)*nf);

            // inverse transform
            fftwf_execute(p);

            // copy to output
            for (int it=0; it<nt; it++) m[itr][it]= add*m[itr][it] + out[it]/sqrt(nt);
        }

        delete [] out;
        fftwf_destroy_plan(p);
        fftwf_free(in);
        fftwf_cleanup();
    }
}

void ifxTransformR::apply_adjoint(bool add, data_t * pmod, const data_t * pdat)
{   
    int nt = _range.getAxis(1).n;
    int ntr = _range.getN123()/nt;
    int nf = _domain.getAxis(1).n/2;
    const data_t (*m) [nt] = (const data_t (*) [nt]) pdat;
    data_t (*d) [2*nf] = (data_t (*) [2*nf]) pmod;

    // double precision
    if ((int)sizeof(data_t)==8){

        double * in = new double[nt];
        fftw_complex * out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf);
        fftw_plan p = fftw_plan_dft_r2c_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt);

            // forward transform
            fftw_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) {
                d[itr][2*f] = add*d[itr][2*f] + 2*out[f][0]/sqrt(nt); // real part
                d[itr][2*f+1] = add*d[itr][2*f+1] + 2*out[f][1]/sqrt(nt) ; // imaginary part
            }
            d[itr][0] -= out[0][0]/sqrt(nt);
        }
        delete [] in;
        fftw_destroy_plan(p);
        fftw_free(out);
        fftw_cleanup();
    }
    else
    {
        float * in = new float[nt];
        fftwf_complex * out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf);
        fftwf_plan p = fftwf_plan_dft_r2c_1d(nt, in, out, FFTW_ESTIMATE);
        
        for (int itr=0; itr<ntr; itr++){

            // transfer input
            memcpy(in, m[itr], sizeof(data_t)*nt);

            // forward transform
            fftwf_execute(p);

            // copy to output
            for (int f=0; f<nf; f++) {
                d[itr][2*f] = add*d[itr][2*f] + 2*out[f][0]/sqrt(nt); // real part
                d[itr][2*f+1] = add*d[itr][2*f+1] + 2*out[f][1]/sqrt(nt) ; // imaginary part
            }
            d[itr][0] -= out[0][0]/sqrt(nt);
        }
        delete [] in;
        fftwf_destroy_plan(p);
        fftwf_free(out);
        fftwf_cleanup();
    }
}

void gradient2d::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    if (!add) memset(pdat,0,_range.getN123()*sizeof(data_t));

    int nz = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ny = _domain.getN123() / (nz*nx);
    data_t dz = 2*_domain.getAxis(1).d;
    data_t dx = 2*_domain.getAxis(2).d;

    data_t (* pm)[nx][nz] = (data_t (*) [nx][nz]) pmod;
    data_t (* pd)[ny][nx][nz] = (data_t (*) [ny][nx][nz]) pdat;

    for (int iy=0; iy<ny; iy++)
    {
        #pragma omp parallel for
        for (int ix=1; ix<nx-1; ix++){
            for (int iz=1; iz<nz-1; iz++){
                pd[0][iy][ix][iz] += _xw*(pm[iy][ix+1][iz] - pm[iy][ix-1][iz])/dx;
                pd[1][iy][ix][iz] += _zw*(pm[iy][ix][iz+1] - pm[iy][ix][iz-1])/dz;
            }
        }
    }
}

void gradient2d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    if (!add) memset(pmod,0,_domain.getN123()*sizeof(data_t));

    int nz = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ny = _domain.getN123() / (nz*nx);
    data_t dz = 2*_domain.getAxis(1).d;
    data_t dx = 2*_domain.getAxis(2).d;

    data_t (* pm)[nx][nz] = (data_t (*) [nx][nz]) pmod;
    data_t (* pd)[ny][nx][nz] = (data_t (*) [ny][nx][nz]) pdat;

    for (int iy=0; iy<ny; iy++)
    {
        for (int ix=1; ix<nx-1; ix++){
            for (int iz=1; iz<nz-1; iz++){
                pm[iy][ix+1][iz] += _xw*pd[0][iy][ix][iz] / dx;
                pm[iy][ix-1][iz] -= _xw*pd[0][iy][ix][iz] / dx;
                pm[iy][ix][iz+1] += _zw*pd[1][iy][ix][iz] / dz;
                pm[iy][ix][iz-1] -= _zw*pd[1][iy][ix][iz] / dz;
            }
        }
    }
}

void laplacian2d::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    if (!add) memset(pdat,0,_range.getN123()*sizeof(data_t));

    int nz = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ny = _domain.getN123() / (nz*nx);
    data_t dz = _domain.getAxis(1).d;
    data_t dx = _domain.getAxis(2).d;
    dz *=dz;
    dx *=dx;

    data_t (* pm)[nx][nz] = (data_t (*) [nx][nz]) pmod;
    data_t (* pd)[nx][nz] = (data_t (*) [nx][nz]) pdat;

    for (int iy=0; iy<ny; iy++)
    {
        #pragma omp parallel for
        for (int ix=1; ix<nx-1; ix++){
            for (int iz=1; iz<nz-1; iz++){
                pd[iy][ix][iz] += _xw*(pm[iy][ix+1][iz] -2*pm[iy][ix][iz] + pm[iy][ix-1][iz])/dx + _zw*(pm[iy][ix][iz+1] -2*pm[iy][ix][iz] + pm[iy][ix][iz-1])/dz;
            }
        }
    }
}

void laplacian2d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    if (!add) memset(pmod,0,_domain.getN123()*sizeof(data_t));

    int nz = _domain.getAxis(1).n;
    int nx = _domain.getAxis(2).n;
    int ny = _domain.getN123() / (nz*nx);
    data_t dz = _domain.getAxis(1).d;
    data_t dx = _domain.getAxis(2).d;
    dz *=dz;
    dx *=dx;

    data_t (* pm)[nx][nz] = (data_t (*) [nx][nz]) pmod;
    data_t (* pd)[nx][nz] = (data_t (*) [nx][nz]) pdat;

    for (int iy=0; iy<ny; iy++)
    {
        for (int ix=1; ix<nx-1; ix++){
            for (int iz=1; iz<nz-1; iz++){
                pm[iy][ix][iz] -= 2*pd[iy][ix][iz]*(_xw/dx + _zw/dz);
                pm[iy][ix+1][iz] += _xw*pd[iy][ix][iz] / dx;
                pm[iy][ix-1][iz] += _xw*pd[iy][ix][iz] / dx;
                pm[iy][ix][iz+1] += _zw*pd[iy][ix][iz] / dz;
                pm[iy][ix][iz-1] += _zw*pd[iy][ix][iz] / dz;
            }
        }
    }
}

void extrapolator1d2d::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    int nz = _range.getAxis(1).n;
    int nx = _range.getAxis(2).n;
    int ny = _range.getN123() / (nz*nx);
    data_t oz = _range.getAxis(1).o;
    data_t dz = _range.getAxis(1).d;

    data_t (* pm)[nz] = (data_t (*) [nz]) pmod;
    data_t (* pd)[nx][nz] = (data_t (*) [nx][nz]) pdat;
    const data_t * ph = _hrz->getCVals();
    for (int iy=0; iy<ny; iy++)
    {
        #pragma omp parallel for
        for (int ix=0; ix<nx; ix++){
            data_t z0=oz+ph[ix]-ph[0];
            if (z0<=oz){
                int iz0min = floor((oz-z0)/dz);
                data_t wmin=((iz0min+1)*dz+z0-oz)/dz;
                for (int iz=0; iz<nz-iz0min-1; iz++) pd[iy][ix][iz] = add*pd[iy][ix][iz] + wmin*pm[iy][iz+iz0min] + (1-wmin)*pm[iy][iz+iz0min+1];
                for (int iz=nz-iz0min-1; iz<nz; iz++) pd[iy][ix][iz] = add*pd[iy][ix][iz] + pm[iy][nz-1];
            }
            else{
                int iz0min = floor((z0-oz)/dz);
                data_t wmin=((iz0min+1)*dz+oz-z0)/dz;
                for (int iz=0; iz<iz0min+1; iz++) pd[iy][ix][iz] = add*pd[iy][ix][iz] + pm[iy][0];
                for (int iz=iz0min+1; iz<nz; iz++) pd[iy][ix][iz] = add*pd[iy][ix][iz] + (1-wmin)*pm[iy][iz-iz0min-1] + wmin*pm[iy][iz-iz0min];
            }
        }
    }
}

void extrapolator1d2d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    if (!add) memset(pmod,0,_domain.getN123()*sizeof(data_t));

    int nz = _range.getAxis(1).n;
    int nx = _range.getAxis(2).n;
    int ny = _range.getN123() / (nz*nx);
    data_t oz = _range.getAxis(1).o;
    data_t dz = _range.getAxis(1).d;

    data_t (* pm)[nz] = (data_t (*) [nz]) pmod;
    data_t (* pd)[nx][nz] = (data_t (*) [nx][nz]) pdat;
    const data_t * ph = _hrz->getCVals();

    for (int iy=0; iy<ny; iy++)
    {
        for (int ix=0; ix<nx; ix++){
            data_t z0=oz+ph[ix]-ph[0];
            if (z0<=oz){
                int iz0min = floor((oz-z0)/dz);
                data_t wmin=((iz0min+1)*dz+z0-oz)/dz;
                for (int iz=0; iz<nz-iz0min-1; iz++) {
                    pm[iy][iz+iz0min] += wmin*pd[iy][ix][iz];
                    pm[iy][iz+iz0min+1] += (1-wmin)*pd[iy][ix][iz];
                }
                for (int iz=nz-iz0min-1; iz<nz; iz++) pm[iy][nz-1] += pd[iy][ix][iz];
            }
            else{
                int iz0min = floor((z0-oz)/dz);
                data_t wmin=((iz0min+1)*dz+oz-z0)/dz;
                for (int iz=0; iz<iz0min+1; iz++) pm[iy][0] += pd[iy][ix][iz];
                for (int iz=iz0min+1; iz<nz; iz++) {
                    pm[iy][iz-iz0min-1] += (1-wmin)*pd[iy][ix][iz];
                    pm[iy][iz-iz0min] += wmin*pd[iy][ix][iz];
                }
            }
        }
    }
}

void extrapolator1d2d::apply_inverse(bool add, data_t * pmod, const data_t * pdat){


    int nz = _range.getAxis(1).n;
    int nx = _range.getAxis(2).n;
    int ny = _range.getN123() / (nz*nx);

    data_t (* pm)[nz] = (data_t (*) [nz]) pmod;
    data_t (* pd)[nx][nz] = (data_t (*) [nx][nz]) pdat;

    for (int iy=0; iy<ny; iy++)
    {
        for (int iz=0; iz<nz; iz++) pm[iy][iz] = add*pm[iy][iz] + pd[iy][0][iz]; 
    }
}

void ssmth::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    if (!add) memset(pdat,0,_domain.getN123()*sizeof(data_t));

    int nf=_domain.getAxis(1).n/2;
    int ntr = _domain.getN123()/(2*nf);

    data_t * pvec = new data_t[_domain.getN123()];
    memset(pvec,0,_domain.getN123()*sizeof(data_t));

    const data_t (*pv1) [2*nf] = (const data_t (*)[2*nf]) pmod;
    data_t (*pv0) [2*nf] = (data_t (*)[2*nf]) pvec;
    data_t (*pv2) [2*nf] = (data_t (*)[2*nf]) pdat;

    #pragma omp parallel for
    for (int itr=0; itr<ntr; itr++){
        for (int f=0; f<nf; f++){
            data_t sc = std::min(nf,f+_hl+1)-std::max(0,f-_hl);
            for (int j=std::max(0,f-_hl); j<std::min(nf,f+_hl+1); j++) pv0[itr][2*f] += pv1[itr][2*j]/sc;
            pv0[itr][2*f+1] = pv1[itr][2*f+1];
        }
    }
    #pragma omp parallel for
    for (int itr=0; itr<ntr; itr++){
        for (int f=0; f<nf; f++){
            data_t sc = std::min(nf,f+_hl+1)-std::max(0,f-_hl);
            for (int j=std::max(0,f-_hl); j<std::min(nf,f+_hl+1); j++) pv2[itr][2*f] += pv0[itr][2*j]/sc;
            pv2[itr][2*f+1] += pv0[itr][2*f+1];
        }
    }
    delete [] pvec;
}
void ssmth::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    if (!add) memset(pmod,0,_domain.getN123()*sizeof(data_t));

    int nf=_domain.getAxis(1).n/2;
    int ntr = _domain.getN123()/(2*nf);

    data_t * pvec = new data_t[_domain.getN123()];
    memset(pvec,0,_domain.getN123()*sizeof(data_t));

    data_t (*pv1) [2*nf] = (data_t (*)[2*nf]) pmod;
    data_t (*pv0) [2*nf] = (data_t (*)[2*nf]) pvec;
    const data_t (*pv2) [2*nf] = (const data_t (*)[2*nf]) pdat;

    #pragma omp parallel for
    for (int itr=0; itr<ntr; itr++){
        for (int f=0; f<nf; f++){
            data_t sc = std::min(nf,f+_hl+1)-std::max(0,f-_hl);
            for (int j=std::max(0,f-_hl); j<std::min(nf,f+_hl+1); j++) pv0[itr][2*j] += pv2[itr][2*f]/sc;
            pv0[itr][2*f+1] = pv2[itr][2*f+1];
        }
    }

    #pragma omp parallel for
    for (int itr=0; itr<ntr; itr++){
        for (int f=0; f<nf; f++){
            data_t sc = std::min(nf,f+_hl+1)-std::max(0,f-_hl);
            for (int j=std::max(0,f-_hl); j<std::min(nf,f+_hl+1); j++) pv1[itr][2*j] += pv0[itr][2*f]/sc;
            pv1[itr][2*f+1] += pv0[itr][2*f+1];
        }
    }
    delete [] pvec;
}

void conv1dnd::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    if (!add) memset(pdat,0,_range.getN123()*sizeof(data_t));

    int nt = _domain.getAxis(1).n;
    int nx = _domain.getN123()/nt;
    int nf = _f->getN123();

    data_t * pf = _f->getVals();
    data_t scale = 1;
    //data_t scale = nt;

    int sh = 0;
    if (_centered) sh = std::floor((nf+1)/2)-1;

    #pragma omp parallel for
    for (int ix=0; ix<nx; ix++){
        for (int it=0; it<nt; it++){
            for (int j=std::max(0,it-nt+1+sh); j<std::min(nf,it+1+sh); j++) pdat[ix*nt+it] += pf[j]*pmod[ix*nt+it-j+sh]/scale;
        }
    }
}

void conv1dnd::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    if (!add) memset(pmod,0,_domain.getN123()*sizeof(data_t));

    int nt = _domain.getAxis(1).n;
    int nx = _domain.getN123()/nt;
    int nf = _f->getN123();

    data_t * pf = _f->getVals();
    data_t scale = 1;
    //data_t scale = nt;

    int sh = 0;
    if (_centered) sh = std::floor((nf+1)/2)-1;

    #pragma omp parallel for
    for (int ix=0; ix<nx; ix++){
        for (int it=0; it<nt; it++){
            for (int j=std::max(0,it-nt+1+sh); j<std::min(nf,it+1+sh); j++) pmod[ix*nt+it-j+sh] += pf[j]*pdat[ix*nt+it]/scale;
        }
    }
}

void convnd1d::apply_forward(bool add, const data_t * pmod, data_t * pdat){

    if (!add) memset(pdat,0,_range.getN123()*sizeof(data_t));

    int nt = _range.getAxis(1).n;
    int nx = _range.getN123()/nt;
    int nf = _domain.getN123();

    data_t * pf = _f->getVals();
    data_t scale = 1;
    //data_t scale = nt;

    int sh = 0;
    if (_centered) sh = std::floor((nf+1)/2)-1;
    
    #pragma omp parallel for
    for (int ix=0; ix<nx; ix++){
        for (int it=0; it<nt; it++){
            for (int j=std::max(0,it-nt+1+sh); j<std::min(nf,it+1+sh); j++) pdat[ix*nt+it] += pmod[j]*pf[ix*nt+it-j+sh]/scale;
        }
    }
}

void convnd1d::apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

    if (!add) memset(pmod,0,_domain.getN123()*sizeof(data_t));

    int nt = _range.getAxis(1).n;
    int nx = _range.getN123()/nt;
    int nf = _domain.getN123();

    data_t * pf = _f->getVals();
    data_t scale = 1; 
    //data_t scale = nt; 

    int sh = 0;
    if (_centered) sh = std::floor((nf+1)/2)-1;

    if (!_centered){
        #pragma omp parallel for
        for (int j=0; j<nf; j++){
            for (int ix=0; ix<nx; ix++){
                for (int it=j; it<nt; it++) pmod[j] += pdat[ix*nt+it]*pf[ix*nt+it-j]/scale;
            }
        }
    }
    else{
        for (int ix=0; ix<nx; ix++){
            for (int it=0; it<nt; it++){
                for (int j=std::max(0,it-nt+1+sh); j<std::min(nf,it+1+sh); j++) pmod[j] += pdat[ix*nt+it]*pf[ix*nt+it-j+sh]/scale;
            }
        }
    }
}





std::shared_ptr<vecReg<data_t> > zero_phase(const std::shared_ptr<vecReg<data_t> > dat){

    axis<data_t>T = dat->getHyper()->getAxis(1);
    axis<data_t>X(1,0,1);
    
    // Make a new vector with odd length
    std::shared_ptr<vecReg<data_t> > vec0;
	if (2*(T.n/2)==T.n) {
		T.n+=1;
        vec0 = std::make_shared<vecReg<data_t> >(hypercube<data_t>(T));

        for (int iz=0; iz<T.n-1; iz++){
        vec0->getVals()[iz] = dat->getVals()[iz];
        }
        vec0->getVals()[T.n-1] = 0;
    }
	
    else {
        vec0 = dat;
        vec0->setHyper((hypercube<data_t>(T,X)));
    }

    // Fourier transform the vector
    fxTransform fx(hypercube<data_t>(T,X));
    axis<data_t>F;
    F.d=fx.getRange()->getAxis(1).d;
    F.o=0;
    F.n=T.n/2 + 1;
    std::shared_ptr<cvecReg<data_t> > fxvec = std::make_shared<cvecReg<data_t> >(hypercube<data_t>(F,X));
    fx.forward(false,vec0,fxvec);

    int N = F.n;
    data_t dw = 2*M_PI/T.n;
    data_t r;
    std::complex<data_t> z (0,0);

    // Modify the FT by setting a linear phase while keeping the same amplitude spectrum
	for (int i=0; i<N; i++){

        r = abs(fxvec->getVals()[i]);

        // add linear phase shift (* exp (-j.omega.dt.(Nt-1)/2)))
        z.real(r*cos(-i*dw*(T.n - 1)/2));
        z.imag(r*sin(-i*dw*(T.n - 1)/2));

        fxvec->getVals()[i] = z;
	}

    // inverse FT and modify the original vector
    fx.inverse(false,vec0,fxvec);
    T.o = (int)(-T.n/2) * T.d; 
    vec0->setHyper(hypercube<data_t>(T));
    
    return vec0;
}
std::shared_ptr<vecReg<data_t> > minimum_phase(const std::shared_ptr<vecReg<data_t> > dat, const data_t eps){
    
    axis<data_t>T = dat->getHyper()->getAxis(1);
    axis<data_t>X(1,0,1);
    
    // Make a new vector with odd length
    std::shared_ptr<vecReg<data_t> > vec0;
	if (2*(T.n/2)==T.n) {
		T.n+=1;
        vec0 = std::make_shared<vecReg<data_t> >(hypercube<data_t>(T));

        for (int iz=0; iz<T.n-1; iz++){
        vec0->getVals()[iz] = dat->getVals()[iz];
        }
        vec0->getVals()[T.n-1] = 0;
    }
	
    else {
        vec0 = dat;
        vec0->setHyper((hypercube<data_t>(T,X)));
    }

    // Fourier transform the vector
    fxTransform fx(hypercube<data_t>(T,X));
    axis<data_t>F;
    F.d=fx.getRange()->getAxis(1).d;
    F.o=0;
    F.n=T.n/2 + 1;
    std::shared_ptr<cvecReg<data_t> > fxvec = std::make_shared<cvecReg<data_t> >(hypercube<data_t>(F,X));
    fx.forward(false,vec0,fxvec);

    int N = F.n;
    data_t r;
    std::complex<data_t> z (0,0);

    // Modify the FT by making it real and equal to the logarithm of the power spectrum
	for (int i=0; i<N; i++){

        r = sqrt(norm(fxvec->getVals()[i]));

        // set FT = log(r+eps)
        z.real(log(r+eps));
        z.imag(0);

        fxvec->getVals()[i] = z;
	}

    // inverse FT and keep the causal part
    fx.inverse(false,vec0,fxvec);
    for (int iz=1; iz<T.n/2+1; iz++){
        vec0->getVals()[iz] *= 2;
    }
    for (int iz=T.n/2+1; iz<T.n; iz++){
        vec0->getVals()[iz] = 0;
    }

    // FT again and take the exponent
    fx.forward(false,vec0,fxvec);
    for (int iz=0; iz<N; iz++){
        fxvec->getVals()[iz] = exp(fxvec->getVals()[iz]);
    }

    // inverse FT 
    fx.inverse(false,vec0,fxvec);
    T.o = 0; 
    vec0->setHyper(hypercube<data_t>(T));
    
    return vec0;
}

#undef EPS