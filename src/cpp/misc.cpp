#include "misc.hpp"
#include <fftw3.h>

void successCheck(bool success, std::string file, int line, std::string msg) {
    if (!success) {
        fprintf(stderr,"\n========ERROR in file %s, line %d========\n",file.c_str(),line);
        throw std::logic_error(msg.c_str());
    }
}

data_t getHcoef(data_t * coef, int size, int n, int i){
    
    if (i<size) return coef[i];
    else if (i>= n-size) return coef[n-1-i]; 
    else return 1;
}

void dipole_to_strain(bool adj, data_t * in, const data_t * dip, int ntr, int nt, int itrmin, int itrmax){

    data_t (*p) [ntr][nt] = (data_t (*) [ntr][nt]) in;
    if (!adj){
        #pragma omp parallel for
        for (int ix=itrmin; ix<itrmax; ix++){
            data_t c = cos(dip[ix]);
            data_t s = sin(dip[ix]);
            for (int it=0; it<nt; it++){
                p[0][ix][it] = c*p[0][ix][it] + s*p[1][ix][it];
                p[1][ix][it] = 0;
            }
        }
    }
    else{
        #pragma omp parallel for
        for (int ix=itrmin; ix<itrmax; ix++){
            data_t c = cos(dip[ix]);
            data_t s = sin(dip[ix]);
            for (int it=0; it<nt; it++){
                p[1][ix][it] = s*p[0][ix][it];
                p[0][ix][it] *= c;
            }
        }
    }
}

void applyHt(bool inv, bool add, const data_t * in, data_t * out, int nx, int nt, data_t dt, int ixmin, int ixmax){

    if (!inv){
        #pragma omp parallel for
        for (int ix=ixmin; ix<ixmax; ix++){
            int i1=ix*nt;
            out[i1] = add*out[i1] + in[i1]*0.5*dt;
            for (int it=1; it<nt-1; it++){
                out[i1+it] = add*out[i1+it] + in[i1+it]*dt;
            }
            out[i1+nt-1] = add*out[i1+nt-1] + in[i1+nt-1]*0.5*dt;
        }
    }
    else{
        #pragma omp parallel for
        for (int ix=ixmin; ix<ixmax; ix++){
            int i1=ix*nt;
            out[i1] = add*out[i1] + in[i1]/(0.5*dt);
            for (int it=1; it<nt-1; it++){
                out[i1+it] = add*out[i1+it] + in[i1+it]/dt;
            }
            out[i1+nt-1] = add*out[i1+nt-1] + in[i1+nt-1]/(0.5*dt);
        }
    }
}

void Dt(bool adj, bool add, const data_t * in, data_t * out, int nx, int nt, data_t dt, int ixmin, int ixmax){

    if (!adj){
        #pragma omp parallel for
        for (int ix=ixmin; ix<ixmax; ix++){
            int i1=ix*nt;
            out[i1] = add*out[i1] + (in[i1+1]-in[i1])/dt;
            for (int it=1; it<nt-1; it++){
                out[i1+it] = add*out[i1+it] + (in[i1+it+1]-in[i1+it-1])/(2*dt);
            }
            out[i1+nt-1] = add*out[i1+nt-1] + (in[i1+nt-1]-in[i1+nt-2])/dt;
        }
    }
    else{
        #pragma omp parallel for
        for (int ix=ixmin; ix<ixmax; ix++){
            int i1=ix*nt;
            out[i1] = add*out[i1] + (in[i1+1]+in[i1])/dt;
            for (int it=1; it<nt-1; it++){
                out[i1+it] = add*out[i1+it] + (in[i1+it+1]-in[i1+it-1])/(2*dt);
            }
            out[i1+nt-1] = add*out[i1+nt-1] + (-in[i1+nt-1]-in[i1+nt-2])/dt;
        }
    }
}

void ttnormalize(data_t * d, data_t * norms, int nt, int ntr, data_t zero){
    #pragma omp parallel for
    for (int ix=0; ix<ntr; ix++){
        norms[ix]=0;
        for (int it=0; it<nt; it++){
            norms[ix] += d[ix*nt+it]*d[ix*nt+it];
        }
        if (norms[ix] <= zero) norms[ix]=1;
        norms[ix] = sqrt(norms[ix]);
        for (int it=0; it<nt; it++){
            d[ix*nt+it] /= norms[ix];
        }
    }
}

void hilbert(std::shared_ptr<vecReg<data_t> > input){
    int n123 = input->getN123();
    std::vector<axis<data_t> > axes = input->getHyper()->getAxes();
    axis<data_t> T = axes[0];
    int nx = round(n123/T.n);
    axis<data_t> X(nx,0,1);
    data_t * pin = input->getVals();

    // double precision
    if ((int)sizeof(data_t)==8){
        fftw_complex *in;
        fftw_plan p1, p2;

        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * T.n);
        p1 = fftw_plan_dft_1d(T.n, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
        p2 = fftw_plan_dft_1d(T.n, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        for (int ix=0; ix<nx; ix++){
            // transfer input to a complex vector
            for (int i=0; i<T.n; i++) {in[i][0]=pin[ix*T.n+i]; in[i][1]=0;}

            // forward transform
            fftw_execute(p1);

            // zero the negative frequencies
            for (int i=(T.n-1)/2; i<T.n; i++) {in[i][0]=0; in[i][1]=0;}

            // inverse transform
            fftw_execute(p2);

            // keep the imaginary part and multiply by 2
            for (int i=0; i<T.n; i++) pin[ix*T.n+i] = 2*in[i][1]/T.n;
        }

        fftw_destroy_plan(p1);
        fftw_destroy_plan(p2);
        fftw_free(in);
        fftw_cleanup();
    }

    // single precision
    else {
        fftwf_complex *in;
        fftwf_plan p1, p2;

        in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * T.n);
        p1 = fftwf_plan_dft_1d(T.n, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
        p2 = fftwf_plan_dft_1d(T.n, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);

        for (int ix=0; ix<nx; ix++){
            // transfer input to a complex vector
            for (int i=0; i<T.n; i++) {in[i][0]=pin[ix*T.n+i]; in[i][1]=0;}

            // forward transform
            fftwf_execute(p1);

            // zero the negative frequencies
            for (int i=(T.n-1)/2; i<T.n; i++) {in[i][0]=0; in[i][1]=0;}

            // inverse transform
            fftwf_execute(p2);

            // keep the imaginary part and multiply by 2
            for (int i=0; i<T.n; i++) pin[ix*T.n+i] = 2*in[i][1]/T.n;
        }

        fftwf_destroy_plan(p1);
        fftwf_destroy_plan(p2);
        fftwf_free(in);
        fftwf_cleanup();
    }
}

void envelop1(std::shared_ptr<vecReg<data_t> > in){
    std::shared_ptr<vecReg<data_t> > hil = in->clone();
    hilbert(hil);
    data_t * pin = in->getVals();
    data_t * phil = hil->getVals();
    for (int i=0; i<in->getN123(); i++) pin[i] = sqrt(phil[i]*phil[i] + pin[i]*pin[i]);
}

void envelop2(std::shared_ptr<vecReg<data_t> > in){
    std::shared_ptr<vecReg<data_t> > hil = in->clone();
    hilbert(hil);
    data_t * pin = in->getVals();
    data_t * phil = hil->getVals();
    for (int i=0; i<in->getN123(); i++) pin[i] = phil[i]*phil[i] + pin[i]*pin[i];
}