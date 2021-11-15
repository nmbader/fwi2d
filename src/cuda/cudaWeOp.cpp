#include "cudaMisc.h"
#include "cudaCppWrappers.h"
#include "we_op.hpp"
#include "misc.hpp"
#include "spatial_operators.hpp"

static inline data_t expr1(const data_t ** par, int i){return par[0][i];} // e.g. = lambda
static inline data_t expr2(const data_t ** par, int i){return par[1][i];} // e.g. = mu
static inline data_t expr3(const data_t ** par, int i){return par[2][i];} // e.g. = rho
static inline data_t expr4a(const data_t ** par, int i){return par[0][i] + 2*par[1][i];} // e.g. = lambda + 2 mu
static inline data_t expr4b(const data_t ** par, int i){return 2*par[1][i];} // e.g. = 2 mu
static inline data_t expr5(const data_t ** par, int i){return 0.5*sqrt(par[1][i]*par[2][i]);} // e.g. = 1/2 sqrt(rho.mu)
static inline data_t expr6(const data_t ** par, int i){return 0.5*sqrt((par[0][i]+2*par[1][i])*par[2][i]);} // e.g. = 1/2 sqrt(rho.(lambda+2 mu))

void nl_we_op_e::propagate_gpu(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz) const
{
    int nxz=nx*nz;
    
    // wavefields allocations and pointers
    data_t *dev_u_prev, *dev_u_curr, *dev_u_next, *dev_dux, *dev_duz, *dev_bucket, *dev_tmp, *dev_grad;
    cudaCheckError( cudaMalloc((void**)&dev_u_prev, 2*nxz*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_u_curr, 2*nxz*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_u_next, 2*nxz*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_dux, 2*nxz*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_duz, 2*nxz*sizeof(data_t)) );
    cudaCheckError( cudaMemset(dev_u_prev, 0, 2*nxz*sizeof(data_t)) );
    cudaCheckError( cudaMemset(dev_u_curr, 0, 2*nxz*sizeof(data_t)) );
    cudaCheckError( cudaMemset(dev_u_next, 0, 2*nxz*sizeof(data_t)) );
    cudaCheckError( cudaMemset(dev_dux, 0, 2*nxz*sizeof(data_t)) );
    cudaCheckError( cudaMemset(dev_duz, 0, 2*nxz*sizeof(data_t)) );

    data_t (*u_full) [2][nxz];
    data_t *dev_u_for;

    if (par.sub>0) u_full = (data_t (*) [2][nxz]) full_wfld;
    if (grad != nullptr) 
    {
        cudaCheckError( cudaMalloc((void**)&dev_tmp, 4*nxz*sizeof(data_t)) );
        cudaCheckError( cudaMalloc((void**)&dev_grad, 3*nxz*sizeof(data_t)) );
        cudaCheckError( cudaMalloc((void**)&dev_u_for, 6*nxz*sizeof(data_t)) );
        cudaCheckError( cudaMemset(dev_tmp, 0, 4*nxz*sizeof(data_t)) );
        cudaCheckError( cudaMemset(dev_grad, 0, 3*nxz*sizeof(data_t)) );
    }

    // copy model to device
	const data_t* mod[3] = {model, model+nxz, model+2*nxz};
    data_t *dev_model;
    cudaCheckError( cudaMalloc((void**)&dev_model, 3*nxz*sizeof(data_t)) );
    cudaCheckError( cudaMemcpy((data_t *)dev_model, model, 3*nxz*sizeof(data_t), cudaMemcpyHostToDevice) );

    // compute P- and S-impedances divided by 2 (1/2 sqrt(rho.(lambda+2 mu)) and sqrt(1/2 rho.mu)) at the boundaries and copy them to device for absorbing BC
    data_t *dev_tb_imp, *dev_lr_imp;
    cudaCheckError( cudaMalloc((void**)&dev_tb_imp, 4*nx*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_lr_imp, 4*nz*sizeof(data_t)) );
    {
        data_t tb_imp[2][2][nx];
        data_t lr_imp[2][2][nz];
        #pragma omp parallel for 
        for (int i=0; i<nx; i++){
            tb_imp[0][0][i] = 0.5*sqrt((mod[0][i*nz]+2*mod[1][i*nz])*mod[2][i*nz]); // top P-impedance
            tb_imp[0][1][i] = 0.5*sqrt(mod[1][i*nz]*mod[2][i*nz]); // top S-impedance
            tb_imp[1][0][i] = 0.5*sqrt((mod[0][i*nz+nz-1]+2*mod[1][i*nz+nz-1])*mod[2][i*nz+nz-1]); // bottom P-impedance
            tb_imp[1][1][i] = 0.5*sqrt(mod[1][i*nz+nz-1]*mod[2][i*nz+nz-1]); // bottom S-impedance
        }
        #pragma omp parallel for 
        for (int i=0; i<nz; i++){
            lr_imp[0][0][i] = 0.5*sqrt((mod[0][i]+2*mod[1][i])*mod[2][i]); // left P-impedance
            lr_imp[0][1][i] = 0.5*sqrt(mod[1][i]*mod[2][i]); // left S-impedance
            lr_imp[1][0][i] = 0.5*sqrt((mod[0][(nx-1)*nz+i]+2*mod[1][(nx-1)*nz+i])*mod[2][(nx-1)*nz+i]); // right P-impedance
            lr_imp[1][1][i] = 0.5*sqrt(mod[1][(nx-1)*nz+i]*mod[2][(nx-1)*nz+i]); // right S-impedance
        }
        
        cudaCheckError( cudaMemcpy(dev_tb_imp, tb_imp, 4*nx*sizeof(data_t), cudaMemcpyHostToDevice) );
        cudaCheckError( cudaMemcpy(dev_lr_imp, lr_imp, 4*nz*sizeof(data_t), cudaMemcpyHostToDevice) );
    }

    // source and receiver components for injection/extraction
    int ns=par.ns, ntr_s=1;
    int nr=par.nr, ntr_r=ext->_ntr;
    if (adj)
    {
        ns = par.nr; ntr_s=inj->_ntr;
        nr = par.ns; ntr_r=1;
    }
    data_t *dev_srcx1, *dev_srcx2, *dev_srcz1, *dev_srcz2;
    data_t *dev_rcvx1, *dev_rcvz1, *dev_rcvx2, *dev_rcvz2;
    data_t **dev_srcx, **dev_srcz;
    data_t **dev_rcvx, **dev_rcvz;
    cudaCheckError( cudaMalloc((void**)&dev_srcx1, par.nt*ntr_s*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_srcz1, par.nt*ntr_s*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_rcvx1, par.nt*ntr_r*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_rcvz1, par.nt*ntr_r*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_srcx, 2*sizeof(data_t*)) );
    cudaCheckError( cudaMalloc((void**)&dev_srcz, 2*sizeof(data_t*)) );
    cudaCheckError( cudaMalloc((void**)&dev_rcvx, 2*sizeof(data_t*)) );
    cudaCheckError( cudaMalloc((void**)&dev_rcvz, 2*sizeof(data_t*)) );
    
    cudaCheckError( cudaMemcpy(dev_srcx1, allsrc, par.nt*ntr_s*sizeof(data_t), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_srcz1, allsrc+ns*par.nt, par.nt*ntr_s*sizeof(data_t), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_rcvx1, allrcv, par.nt*ntr_r*sizeof(data_t), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_rcvz1, allrcv+nr*par.nt, par.nt*ntr_r*sizeof(data_t), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_srcx, &dev_srcx1, sizeof(data_t*), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_srcz, &dev_srcz1, sizeof(data_t*), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_rcvx, &dev_rcvx1, sizeof(data_t*), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_rcvz, &dev_rcvz1, sizeof(data_t*), cudaMemcpyHostToDevice) );

    const data_t * srcx[2] = {allsrc, nullptr};
    const data_t * srcz[2] = {allsrc + ns*par.nt, nullptr};
    data_t * rcvx[2] = {allrcv, nullptr};
    data_t * rcvz[2] = {allrcv + nr*par.nt, nullptr};
    if (par.mt==true)
    {
        if (!adj)
        {
            cudaCheckError( cudaMalloc((void**)&dev_srcx2, par.nt*ntr_s*sizeof(data_t)) );
            cudaCheckError( cudaMalloc((void**)&dev_srcz2, par.nt*ntr_s*sizeof(data_t)) );
            cudaCheckError( cudaMemcpy(dev_srcx2, allsrc+2*ns*par.nt, par.nt*ntr_s*sizeof(data_t), cudaMemcpyHostToDevice) );
            cudaCheckError( cudaMemcpy(dev_srcz1, allsrc+2*ns*par.nt, par.nt*ntr_s*sizeof(data_t), cudaMemcpyHostToDevice) );
            cudaCheckError( cudaMemcpy(dev_srcz2, allsrc+ns*par.nt, par.nt*ntr_s*sizeof(data_t), cudaMemcpyHostToDevice) );
            
            cudaCheckError( cudaMemcpy(dev_srcx+1, &dev_srcx2, sizeof(data_t*), cudaMemcpyHostToDevice) );
            cudaCheckError( cudaMemcpy(dev_srcz+1, &dev_srcz2, sizeof(data_t*), cudaMemcpyHostToDevice) );
            
            srcx[1] = allsrc + 2*ns*par.nt;
            srcz[0] = allsrc + 2*ns*par.nt;
            srcz[1] = allsrc + ns*par.nt;
        }
        else
        {
            cudaCheckError( cudaMalloc((void**)&dev_rcvx2, par.nt*ntr_r*sizeof(data_t)) );
            cudaCheckError( cudaMalloc((void**)&dev_rcvz2, par.nt*ntr_r*sizeof(data_t)) );
            cudaCheckError( cudaMemcpy(dev_rcvx2, allrcv+2*nr*par.nt, par.nt*ntr_r*sizeof(data_t), cudaMemcpyHostToDevice) );
            cudaCheckError( cudaMemcpy(dev_rcvz1, allrcv+2*nr*par.nt, par.nt*ntr_r*sizeof(data_t), cudaMemcpyHostToDevice) );
            cudaCheckError( cudaMemcpy(dev_rcvz2, allrcv+nr*par.nt, par.nt*ntr_r*sizeof(data_t), cudaMemcpyHostToDevice) );
            
            cudaCheckError( cudaMemcpy(dev_rcvx+1, &dev_rcvx2, sizeof(data_t*), cudaMemcpyHostToDevice) );
            cudaCheckError( cudaMemcpy(dev_rcvz+1, &dev_rcvz2, sizeof(data_t*), cudaMemcpyHostToDevice) );
            
            rcvx[1] = allrcv + 2*nr*par.nt;
            rcvz[0] = allrcv + 2*nr*par.nt;
            rcvz[1] = allrcv + nr*par.nt;
        }
    }

    // copy sources/receivers locations and weights to device
    data_t *dev_inj_xw, *dev_inj_zw, *dev_ext_xw, *dev_ext_zw;
    int *dev_inj_xind, *dev_inj_zind, *dev_ext_xind, *dev_ext_zind;
    cudaCheckError( cudaMalloc((void**)&dev_inj_xw, inj->_xw.size()*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_inj_zw, inj->_zw.size()*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_ext_xw, ext->_xw.size()*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_ext_zw, ext->_zw.size()*sizeof(data_t)) );
    cudaCheckError( cudaMalloc((void**)&dev_inj_xind, inj->_xind.size()*sizeof(int)) );
    cudaCheckError( cudaMalloc((void**)&dev_inj_zind, inj->_zind.size()*sizeof(int)) );
    cudaCheckError( cudaMalloc((void**)&dev_ext_xind, ext->_xind.size()*sizeof(int)) );
    cudaCheckError( cudaMalloc((void**)&dev_ext_zind, ext->_zind.size()*sizeof(int)) );
    cudaCheckError( cudaMemcpy(dev_inj_xw, inj->_xw.data(), inj->_xw.size()*sizeof(data_t), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_inj_zw, inj->_zw.data(), inj->_zw.size()*sizeof(data_t), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_ext_xw, ext->_xw.data(), ext->_xw.size()*sizeof(data_t), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_ext_zw, ext->_zw.data(), ext->_zw.size()*sizeof(data_t), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_inj_xind, inj->_xind.data(), inj->_xind.size()*sizeof(int), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_inj_zind, inj->_zind.data(), inj->_zind.size()*sizeof(int), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_ext_xind, ext->_xind.data(), ext->_xind.size()*sizeof(int), cudaMemcpyHostToDevice) );
    cudaCheckError( cudaMemcpy(dev_ext_zind, ext->_zind.data(), ext->_zind.size()*sizeof(int), cudaMemcpyHostToDevice) );


    // prev = 1/(2 * rho) * dt2 * src
    data_t * wfld = new data_t[2*nxz];
    memset(wfld,0,2*nxz*sizeof(data_t));
    inj->inject(false, srcx, wfld, nx, nz, par.nt, inj->_ntr, 0, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
    inj->inject(false, srcz, wfld+nxz, nx, nz, par.nt, inj->_ntr, 0, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
    #pragma omp parallel for
    for (int i=0; i<nxz; i++)
    {
        wfld[i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
        wfld[nxz+i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
    }
    cudaCheckError( cudaMemcpy(dev_u_prev, wfld, 2*nxz*sizeof(data_t), cudaMemcpyHostToDevice) );
    delete [] wfld;

    int pct10 = round(par.nt/10);

    for (int it=0; it<par.nt-1; it++)
    {
        // copy the current wfld to the full wfld vector
        if ((par.sub>0) && (it%par.sub==0) && (grad==nullptr)) cudaCheckError( cudaMemcpy(u_full[it/par.sub], dev_u_curr, 2*nxz*sizeof(data_t), cudaMemcpyDeviceToHost) );

        // extract receivers
        if (grad == nullptr)
        {
            ext->extract_gpu(true, dev_u_curr, dev_rcvx, nx, nz, par.nt, ext->_ntr, it, 0, ext->_ntr, dev_ext_xind, dev_ext_zind, dev_ext_xw, dev_ext_zw);
            ext->extract_gpu(true, dev_u_curr+nxz, dev_rcvz, nx, nz, par.nt, ext->_ntr, it, 0, ext->_ntr, dev_ext_xind, dev_ext_zind, dev_ext_xw, dev_ext_zw);
        }

        // compute FWI gradients except for first and last time samples
        if ((grad != nullptr) && (it%par.sub==0) && it!=0) {
            cudaCheckError( cudaMemcpy(dev_u_for, u_full[par.nt/par.sub+1-it], 6*nxz*sizeof(data_t), cudaMemcpyHostToDevice) );
            compute_gradients_gpu(model, dev_u_for, dev_u_curr, dev_dux, dev_duz, dev_tmp, dev_grad, par, nx, nz, it/par.sub, dx, dz, par.sub*par.dt);
        }

        // apply spatial SBP operators
        Dxx_var_gpu(false, dev_u_curr, dev_u_next, nx, nz, dx, dev_model+nxz , 2);
        Dzz_var_gpu(false, dev_u_curr+nxz, dev_u_next+nxz, nx, nz, dz, dev_model+nxz , 2);
        cudaDeviceSynchronize();
        
        Dxx_var_gpu(true, dev_u_curr+nxz, dev_u_next+nxz, nx, nz, dx, dev_model+nxz , 1);
        Dzz_var_gpu(true, dev_u_curr, dev_u_next, nx, nz, dz, dev_model+nxz , 1);
        cudaDeviceSynchronize();

        Dx_gpu(false, dev_u_curr, dev_dux, nx, nz, dx);
        Dz_gpu(false, dev_u_curr+nxz, dev_duz+nxz, nx, nz, dz);
        cudaDeviceSynchronize();
        
        mult_Dx_gpu(true, dev_dux, dev_u_next, nx, nz, dx, dev_model, 1.0);
        mult_Dz_gpu(true, dev_duz+nxz, dev_u_next+nxz, nx, nz, dz, dev_model, 1.0);
        cudaDeviceSynchronize();
        
        mult_Dx_gpu(true, dev_duz+nxz, dev_u_next, nx, nz, dx, dev_model, 1.0);
        mult_Dz_gpu(true, dev_dux, dev_u_next+nxz, nx, nz, dz, dev_model, 1.0);
        cudaDeviceSynchronize();

        Dx_gpu(false, dev_u_curr+nxz, dev_dux+nxz, nx, nz, dx);
        Dz_gpu(false, dev_u_curr, dev_duz, nx, nz, dz);
        cudaDeviceSynchronize();
        
        mult_Dx_gpu(true, dev_duz, dev_u_next+nxz, nx, nz, dx, dev_model+nxz, 1.0);
        mult_Dz_gpu(true, dev_dux+nxz, dev_u_next, nx, nz, dz, dev_model+nxz, 1.0);
        cudaDeviceSynchronize();

        // inject sources
        inj->inject_gpu(true, (const data_t**)dev_srcx, dev_u_next, nx, nz, par.nt, inj->_ntr, it, 0, inj->_ntr, dev_inj_xind, dev_inj_zind, dev_inj_xw, dev_inj_zw);
        inj->inject_gpu(true, (const data_t**)dev_srcz, dev_u_next+nxz, nx, nz, par.nt, inj->_ntr, it, 0, inj->_ntr, dev_inj_xind, dev_inj_zind, dev_inj_xw, dev_inj_zw);

        // apply boundary conditions
        if (par.bc_top==1)
        {
            esat_neumann_top_gpu(true, dev_u_curr+nxz, dev_u_curr, dev_u_next, nx, nz, dx, dz, dev_model+nxz, dev_model+nxz, 1, 1);
            esat_neumann_top_gpu(true, dev_u_curr, dev_u_curr+nxz, dev_u_next+nxz, nx, nz, dx, dz, dev_model, dev_model+nxz, 1, 2);
            esat_Dz_top_gpu(true, dev_u_curr+nxz, dev_u_next+nxz, nx, nz, dz, dev_model, 1.0);
        }
        else if (par.bc_top==2)
        {
            esat_absorbing_top_gpu(true, dev_u_curr+nxz, dev_u_curr, dev_u_prev, dev_u_next, nx, nz, dx, dz, par.dt, dev_model+nxz, dev_model+nxz, dev_tb_imp+nx, 1, 1, 1);
            esat_absorbing_top_gpu(true, dev_u_curr, dev_u_curr+nxz, dev_u_prev+nxz, dev_u_next+nxz, nx, nz, dx, dz, par.dt, dev_model, dev_model+nxz, dev_tb_imp, 1, 2, 1);
            esat_Dz_top_gpu(true, dev_u_curr+nxz, dev_u_next+nxz, nx, nz, dz, dev_model, 1.0);
        }
        if (par.bc_bottom==1)
        {
            esat_neumann_bottom_gpu(true, dev_u_curr+nxz, dev_u_curr, dev_u_next, nx, nz, dx, dz, dev_model+nxz, dev_model+nxz, 1, 1);
            esat_neumann_bottom_gpu(true, dev_u_curr, dev_u_curr+nxz, dev_u_next+nxz, nx, nz, dx, dz, dev_model, dev_model+nxz, 1, 2);
            esat_Dz_bottom_gpu(true, dev_u_curr+nxz, dev_u_next+nxz, nx, nz, dz, dev_model, 1.0);
        }
        else if (par.bc_bottom==2)
        {
            esat_absorbing_bottom_gpu(true, dev_u_curr+nxz, dev_u_curr, dev_u_prev, dev_u_next, nx, nz, dx, dz, par.dt, dev_model+nxz, dev_model+nxz, dev_tb_imp+3*nx, 1, 1, 1);
            esat_absorbing_bottom_gpu(true, dev_u_curr, dev_u_curr+nxz, dev_u_prev+nxz, dev_u_next+nxz, nx, nz, dx, dz, par.dt, dev_model, dev_model+nxz, dev_tb_imp+2*nx, 1, 2, 1);
            esat_Dz_bottom_gpu(true, dev_u_curr+nxz, dev_u_next+nxz, nx, nz, dz, dev_model, 1.0);
        }
        if (par.bc_left==1)
        {
            esat_neumann_left_gpu(true, dev_u_curr+nxz, dev_u_curr, dev_u_next, nx, nz, dx, dz, dev_model, dev_model+nxz, 1, 2);
            esat_Dx_left_gpu(true, dev_u_curr, dev_u_next, nx, nz, dx, dev_model, 1.0);
            esat_neumann_left_gpu(true, dev_u_curr, dev_u_curr+nxz, dev_u_next+nxz, nx, nz, dx, dz, dev_model+nxz, dev_model+nxz, 1, 1);
        }
        else if (par.bc_left==2)
        {
            esat_absorbing_left_gpu(true, dev_u_curr+nxz, dev_u_curr, dev_u_prev, dev_u_next, nx, nz, dx, dz, par.dt, dev_model, dev_model+nxz, dev_lr_imp, 1, 2, 1);
            esat_Dx_left_gpu(true, dev_u_curr, dev_u_next, nx, nz, dx, dev_model, 1.0);
            esat_absorbing_left_gpu(true, dev_u_curr, dev_u_curr+nxz, dev_u_prev+nxz, dev_u_next+nxz, nx, nz, dx, dz, par.dt, dev_model+nxz, dev_model+nxz, dev_lr_imp+nz, 1, 1, 1);
        }
        if (par.bc_right==1)
        {
            esat_neumann_right_gpu(true, dev_u_curr+nxz, dev_u_curr, dev_u_next, nx, nz, dx, dz, dev_model, dev_model+nxz, 1, 2);
            esat_Dx_right_gpu(true, dev_u_curr, dev_u_next, nx, nz, dx, dev_model, 1.0);
            esat_neumann_right_gpu(true, dev_u_curr, dev_u_curr+nxz, dev_u_next+nxz, nx, nz, dx, dz, dev_model+nxz, dev_model+nxz, 1, 1);
        }
        else if (par.bc_right==2)
        {
            esat_absorbing_right_gpu(true, dev_u_curr+nxz, dev_u_curr, dev_u_prev, dev_u_next, nx, nz, dx, dz, par.dt, dev_model, dev_model+nxz, dev_lr_imp+2*nz, 1, 2, 1);
            esat_Dx_right_gpu(true, dev_u_curr, dev_u_next, nx, nz, dx, dev_model, 1.0);
            esat_absorbing_right_gpu(true, dev_u_curr, dev_u_curr+nxz, dev_u_prev+nxz, dev_u_next+nxz, nx, nz, dx, dz, par.dt, dev_model+nxz, dev_model+nxz, dev_lr_imp+3*nz, 1, 1, 1);
        }

        // update wfld with the 2 steps time recursion
        time_step_gpu(dev_u_prev, dev_u_curr, dev_u_next, dev_model+2*nxz, nx, nz, par.dt);
        time_step_gpu(dev_u_prev+nxz, dev_u_curr+nxz, dev_u_next+nxz, dev_model+2*nxz, nx, nz, par.dt);

        // scale boundaries when relevant (for locally absorbing BC only)
        esat_scale_boundaries_gpu(dev_u_next, nx, nz, dx, dz, dev_model, par.dt, par.bc_top==2, par.bc_bottom==2, par.bc_left==2, par.bc_right==2);
        
        // apply taper when relevant
        if (par.taper_top>0) {
            taper_top_gpu(dev_u_curr, nx, nz, par.taper_top, par.taper_strength);
            taper_top_gpu(dev_u_curr+nxz, nx, nz, par.taper_top, par.taper_strength);
            taper_top_gpu(dev_u_next, nx, nz, par.taper_top, par.taper_strength);
            taper_top_gpu(dev_u_next+nxz, nx, nz, par.taper_top, par.taper_strength);
        }
        if (par.taper_bottom>0) {
            taper_bottom_gpu(dev_u_curr, nx, nz, par.taper_bottom, par.taper_strength);
            taper_bottom_gpu(dev_u_curr+nxz, nx, nz, par.taper_bottom, par.taper_strength);
            taper_bottom_gpu(dev_u_next, nx, nz, par.taper_bottom, par.taper_strength);
            taper_bottom_gpu(dev_u_next+nxz, nx, nz, par.taper_bottom, par.taper_strength);
        }
        if (par.taper_left>0) {
            taper_left_gpu(dev_u_curr, nx, nz, par.taper_left, par.taper_strength);
            taper_left_gpu(dev_u_curr+nxz, nx, nz, par.taper_left, par.taper_strength);
            taper_left_gpu(dev_u_next, nx, nz, par.taper_left, par.taper_strength);
            taper_left_gpu(dev_u_next+nxz, nx, nz, par.taper_left, par.taper_strength);
        }
        if (par.taper_right>0) {
            taper_right_gpu(dev_u_curr, nx, nz, par.taper_right, par.taper_strength);
            taper_right_gpu(dev_u_curr+nxz, nx, nz, par.taper_right, par.taper_strength);
            taper_right_gpu(dev_u_next, nx, nz, par.taper_right, par.taper_strength);
            taper_right_gpu(dev_u_next+nxz, nx, nz, par.taper_right, par.taper_strength);
        }

        dev_bucket=dev_u_prev;
        dev_u_prev=dev_u_curr;
        dev_u_curr=dev_u_next;
        dev_u_next=dev_bucket;

        if ((it+1) % pct10 == 0 && par.verbose>2) fprintf(stderr,"Propagation progress = %d\%\n",10*(it+1)/pct10);
    }

    // copy the last wfld to the full wfld vector
    if ((par.sub>0) && ((par.nt-1)%par.sub==0) && (grad==nullptr)) cudaCheckError( cudaMemcpy(u_full[(par.nt-1)/par.sub], dev_u_curr, 2*nxz*sizeof(data_t), cudaMemcpyDeviceToHost) );

    // extract receivers last sample then copy back to host
    if (grad == nullptr)
    {
        ext->extract_gpu(true, dev_u_curr, dev_rcvx, nx, nz, par.nt, ext->_ntr, par.nt-1, 0, ext->_ntr, dev_ext_xind, dev_ext_zind, dev_ext_xw, dev_ext_zw);
        ext->extract_gpu(true, dev_u_curr+nxz, dev_rcvz, nx, nz, par.nt, ext->_ntr, par.nt-1, 0, ext->_ntr, dev_ext_xind, dev_ext_zind, dev_ext_xw, dev_ext_zw);
        cudaCheckError( cudaMemcpy(allrcv, dev_rcvx1, par.nt*ntr_r*sizeof(data_t), cudaMemcpyDeviceToHost) );
        cudaCheckError( cudaMemcpy(allrcv+nr*par.nt, dev_rcvz1, par.nt*ntr_r*sizeof(data_t), cudaMemcpyDeviceToHost) );
        if (par.mt==true &&  adj==true){
            cudaCheckError( cudaMemcpy(allrcv+2*nr*par.nt, dev_rcvx2, par.nt*ntr_r*sizeof(data_t), cudaMemcpyDeviceToHost) );
            cudaCheckError( cudaMemcpy(allrcv+2*nr*par.nt, dev_rcvz1, par.nt*ntr_r*sizeof(data_t), cudaMemcpyDeviceToHost) );
            cudaCheckError( cudaMemcpy(allrcv+nr*par.nt, dev_rcvz2, par.nt*ntr_r*sizeof(data_t), cudaMemcpyDeviceToHost) );
        }
    }
    else
    {
        cudaCheckError( cudaMemcpy(grad, dev_grad, 3*nxz*sizeof(data_t), cudaMemcpyDeviceToHost) );
    }
    
    cudaCheckError( cudaFree(dev_u_prev) );
    cudaCheckError( cudaFree(dev_u_curr) );
    cudaCheckError( cudaFree(dev_u_next) );
    cudaCheckError( cudaFree(dev_dux) );
    cudaCheckError( cudaFree(dev_duz) );
    cudaCheckError( cudaFree(dev_model) );
    cudaCheckError( cudaFree(dev_tb_imp) );
    cudaCheckError( cudaFree(dev_lr_imp) );

    cudaCheckError( cudaFree(dev_srcx1) );
    cudaCheckError( cudaFree(dev_srcz1) );
    cudaCheckError( cudaFree(dev_rcvx1) );
    cudaCheckError( cudaFree(dev_rcvz1) );
    cudaCheckError( cudaFree(dev_srcx) );
    cudaCheckError( cudaFree(dev_srcz) );
    cudaCheckError( cudaFree(dev_rcvx) );
    cudaCheckError( cudaFree(dev_rcvz) );
    
    cudaCheckError( cudaFree(dev_inj_xw) );
    cudaCheckError( cudaFree(dev_inj_zw) );
    cudaCheckError( cudaFree(dev_ext_xw) );
    cudaCheckError( cudaFree(dev_ext_zw) );
    cudaCheckError( cudaFree(dev_inj_xind) );
    cudaCheckError( cudaFree(dev_inj_zind) );
    cudaCheckError( cudaFree(dev_ext_xind) );
    cudaCheckError( cudaFree(dev_ext_zind) );
    
    if (par.mt==true){
        if (!adj){
            cudaCheckError( cudaFree(dev_srcx2) );
            cudaCheckError( cudaFree(dev_srcz2) );
        }
        else{
            cudaCheckError( cudaFree(dev_rcvx2) );
            cudaCheckError( cudaFree(dev_rcvz2) );
        }
    }
    if (grad != nullptr) {
        cudaCheckError( cudaFree(dev_tmp) );
        cudaCheckError( cudaFree(dev_grad) );
        cudaCheckError( cudaFree(dev_u_for) );
    }
}