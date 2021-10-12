#include "misc.hpp"

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
