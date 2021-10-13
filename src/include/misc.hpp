#pragma once

#include "vecReg.hpp"

void successCheck(bool success, std::string file, int line, std::string msg);

// get the coefficient of the spatial quadrature operator given the index in the array
// the coefficients do not include the spatial sampling
data_t getHcoef(data_t * coef, int size, int n, int i);

// convert 2-components dipole data u(gl/2) - u(-gl/2) to strain. The strain is stored in the first component
// the adjoint performs the opposite transformation
// forward: xout = cos(dip).xin + sin(dip).zin  ; zout=0
// adjoint: xout = cos(dip).xin ; zout = sin(dip).xin
void dipole_to_strain(bool adj, data_t * in, const data_t * dip, int ntr, int nt, int itrmin, int itrmax);

// trapezoidal time quadrature operator Ht (or its inverse)
void applyHt(bool inv, bool add, const data_t * in, data_t * out, int nx, int nt, data_t dt, int ixmin, int ixmax);

// First order time derivative operator (consistent with Ht trapezoidal quadrature)
// The adjoint here is actually Ht^-1.Dt^T.Ht in reverse time, to be used in adjoint wave propagation only
void Dt(bool adj, bool add, const data_t * in, data_t * out, int nx, int nt, data_t dt, int ixmin, int ixmax);

// normalize a data trace by trace, store the norms in a vector
void ttnormalize(data_t * d, data_t * norms, int nt, int ntr, data_t zero=1e-16);

// compute the Hilbert transform of a group of traces
void hilbert(std::shared_ptr<vecReg<data_t> > input);

// compute the envelop of a signal
void envelop1(std::shared_ptr<vecReg<data_t> > in);

// compute the envelop squared of a signal
void envelop2(std::shared_ptr<vecReg<data_t> > in);

// Perform in-place LU factorization without pivoting
// The size of the matrix A is n x n, stored in column major
// L has 1s on the diagonal
template<typename T>
void LU(T * A, int n){
    
    // loop over columns
    for (int k=0; k<n; k++){
        if (A[k*n+k] == 0) {
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::runtime_error("Zero pivot has been encountered during LU factorization.\n");
        }

        // Fill L below the diagonal
        for (int i=k+1; i<n; i++){
            A[k*n+i] /= A[k*n+k];
        }

        // subtract from the lower right block the outer product between column k and row k
        for (int i=k+1; i<n; i++){
            for (int j=k+1; j<n; j++){
                A[j*n+i] -= A[k*n+i]*A[j*n+k];
            }
        }
    }
}

// Solve Ax=b using the in-place LU factorization
// n is the size of the vectors
template<typename T>
void solve_Axb(T* A, T*x, T*b, int n){
    
    // factorize A
    LU(A,n);

    // solve Lz = b
    std::vector<T> z(n);
    memcpy(x,b,sizeof(T)*n);
    for (int j=0; j<n; j++){
        for (int i=j+1; i<n; i++){
            x[i] -= A[j*n+i] * x[j];
        }
    }

    // solve Ux = z
    for (int j=n-1; j>=0; j--){
        x[j] = x[j] / A[j*n+j];
        for (int i=j-1; i>=0; i--){
            x[i] -= A[j*n+i] * x[j];
        }
    }
}