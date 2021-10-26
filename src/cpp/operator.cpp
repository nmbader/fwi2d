#include "operator.hpp"
#include <fftw3.h>

void nloper::dotProduct(){

    std::shared_ptr<vecReg<data_t> > m (new vecReg<data_t>(_domain));
    std::shared_ptr<vecReg<data_t> > dm (new vecReg<data_t>(_domain));
    std::shared_ptr<vecReg<data_t> > d (new vecReg<data_t>(_range));
    std::shared_ptr<vecReg<data_t> > m1;

    m->random(-10,10);
    d->random(-10,10);
    dm = m->clone();

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
    jacobian(false,mtilde,m,d); // (df/dm)'.d

    data_t sum1 = d->dot(dtilde);
    data_t sum2 = dm->dot(mtilde);

    fprintf(stderr,"Dot product with eps = %f:\n",eps);
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));

    eps = 0.01;

    m1 = dm->clone();
    m1->scaleAdd(m,eps,1); // m + eps.dm
    forward(false,m1,dtilde); // f(m + eps.dm)
    dtilde->scaleAdd(d1,1,-1); // f(m + eps.dm) - f(m)
    dtilde->scale(1.0/eps); // (df/dm).dm
    jacobian(false,mtilde,m,d); // (df/dm)'.d

    sum1 = d->dot(dtilde);
    sum2 = dm->dot(mtilde);

    fprintf(stderr,"Dot product with eps = %f:\n",eps);
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));

    eps = 0.001;
    
    m1 = dm->clone();
    m1->scaleAdd(m,eps,1); // m + eps.dm
    forward(false,m1,dtilde); // f(m + eps.dm)
    dtilde->scaleAdd(d1,1,-1); // f(m + eps.dm) - f(m)
    dtilde->scale(1.0/eps); // (df/dm).dm
    jacobian(false,mtilde,m,d); // (df/dm)'.d

    sum1 = d->dot(dtilde);
    sum2 = dm->dot(mtilde);

    fprintf(stderr,"Dot product with eps = %f:\n",eps);
    fprintf(stderr,"sum1 = %f\nsum2 = %f\ndiff (x1000) = %f\n",sum1,sum2,1000*(sum1 - sum2));   
}

void loper::dotProduct(){

    std::shared_ptr<vecReg<data_t> > m (new vecReg<data_t>(_domain));
    std::shared_ptr<vecReg<data_t> > d (new vecReg<data_t>(_range));

    //time_t timer; unsigned int seed=time(&timer);
    m->random();
    d->random();

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

void fxTransform::forward(bool add, const std::shared_ptr<vecReg<data_t> > mod, std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int ntr = _domain.getN123()/nt;
    int nf = _range.getAxis(1).n;
    const data_t (*m) [nt] = (const data_t (*) [nt]) mod->getCVals();
    std::complex<data_t> (*d) [nf] = (std::complex<data_t> (*) [nf]) dat->getVals();

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

void fxTransform::inverse(bool add, std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<cvecReg<data_t> > dat)
{
    successCheck(checkDomainRange(mod,dat),__FILE__,__LINE__,"Vectors hypercube do not match the operator domain and range\n");

    int nt = _domain.getAxis(1).n;
    int ntr = _domain.getN123()/nt;
    int nf = _range.getAxis(1).n;
    data_t (*m) [nt] = (data_t (*) [nt]) mod->getCVals();
    const std::complex<data_t> (*d) [nf] = (const std::complex<data_t> (*) [nf]) dat->getVals();

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

    std::shared_ptr<vecReg<data_t> > mtild = m->clone();
    std::shared_ptr<cvecReg<data_t> > dtild = d->clone();

    this->forward(false, m, dtild);
    this->adjoint(false, mtild, d);

    data_t sum1 = (dtild->dot(d)).real();
    data_t sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=false:\n");
    fprintf(stderr,"sum1 = %f\nsum2 = %f\n",sum1, sum2);

    this->forward(true, m, dtild);
    this->adjoint(true, mtild, d);

    sum1 = (dtild->dot(d)).real();
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
    this->cadjoint(false, mtild, d);

    std::complex<data_t> sum1 = dtild->dot(d);
    std::complex<data_t> sum2 = m->dot(mtild);

    fprintf(stderr,"Dot product with add=false:\n");
    fprintf(stderr,"sum1 = %f + i%f\nsum2 = %f + i%f\n",sum1.real(), sum1.imag(), sum2.real(), sum2.imag());

    this->cforward(true, m, dtild);
    this->cadjoint(true, mtild, d);

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