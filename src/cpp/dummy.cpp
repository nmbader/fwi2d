#pragma once

#include "injector.hpp"
#include "we_op.hpp"

void delta_m3::inject_gpu(bool add, const data_t ** in, data_t * out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw) const{}
void delta_m3::extract_gpu(bool add, const data_t * in, data_t ** out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw) const{}
void ddelta_m3::inject_gpu(bool add, const data_t ** in, data_t * out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw) const{}
void ddelta_m3::extract_gpu(bool add, const data_t * in, data_t ** out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw) const{}
void dipole_m3::inject_gpu(bool add, const data_t ** in, data_t * out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw) const{}
void dipole_m3::extract_gpu(bool add, const data_t * in, data_t ** out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw) const{}

void nl_we_op_e::compute_gradients_gpu(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const{}
void nl_we_op_e::propagate_gpu(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz) const{}

void nl_we_op_vti::compute_gradients_gpu(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const{}
void nl_we_op_vti::propagate_gpu(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz) const{}
