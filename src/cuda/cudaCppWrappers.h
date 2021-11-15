#pragma once

#include "vecReg.hpp"

void Dz_gpu(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d);
void Dx_gpu(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d);
void mult_Dz_gpu(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t * par, data_t a);
void mult_Dx_gpu(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t * par, data_t a);
void Dzz_var_gpu(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t * par, data_t a);
void Dxx_var_gpu(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t * par, data_t a);
void esat_neumann_top_gpu(bool add, const data_t* in0, const data_t*in1, data_t* out, int nx, int nz, data_t dx, data_t dz, const data_t * par0, const data_t * par1, data_t a0, data_t a1);
void esat_neumann_bottom_gpu(bool add, const data_t* in0, const data_t*in1, data_t* out, int nx, int nz, data_t dx, data_t dz, const data_t * par0, const data_t * par1, data_t a0, data_t a1);
void esat_neumann_left_gpu(bool add, const data_t* in0, const data_t*in1, data_t* out, int nx, int nz, data_t dx, data_t dz, const data_t * par0, const data_t * par1, data_t a0, data_t a1);
void esat_neumann_right_gpu(bool add, const data_t* in0, const data_t*in1, data_t* out, int nx, int nz, data_t dx, data_t dz, const data_t * par0, const data_t * par1, data_t a0, data_t a1);
void esat_absorbing_top_gpu(bool add, const data_t* in0, const data_t* in1, const data_t* in2, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, const data_t * par0,  const data_t * par1, const data_t * par2, data_t a0, data_t a1, data_t a2);
void esat_absorbing_bottom_gpu(bool add, const data_t* in0, const data_t* in1, const data_t* in2, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, const data_t * par0,  const data_t * par1, const data_t * par2, data_t a0, data_t a1, data_t a2);
void esat_absorbing_left_gpu(bool add, const data_t* in0, const data_t* in1, const data_t* in2, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, const data_t * par0,  const data_t * par1, const data_t * par2, data_t a0, data_t a1, data_t a2);
void esat_absorbing_right_gpu(bool add, const data_t* in0, const data_t* in1, const data_t* in2, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, const data_t * par0,  const data_t * par1, const data_t * par2, data_t a0, data_t a1, data_t a2);
void esat_Dz_top_gpu(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, const data_t * par, data_t a);
void esat_Dz_bottom_gpu(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, const data_t * par, data_t a);
void esat_Dx_left_gpu(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, const data_t * par, data_t a);
void esat_Dx_right_gpu(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, const data_t * par, data_t a);
void esat_scale_boundaries_gpu(data_t* in, int nx, int nz, data_t dx, data_t dz, const data_t* par, data_t dt, bool top, bool bottom, bool left, bool right);
void vtisat_scale_boundaries_gpu(data_t* in, int nx, int nz, data_t dx, data_t dz, const data_t* par, data_t dt, bool top, bool bottom, bool left, bool right);
void taper_top_gpu(data_t* in, int nx, int nz, int taper, data_t a);
void taper_bottom_gpu(data_t* in, int nx, int nz, int taper, data_t a);
void taper_left_gpu(data_t* in, int nx, int nz, int taper, data_t a);
void taper_right_gpu(data_t* in, int nx, int nz, int taper, data_t a);
void time_step_gpu(const data_t * prev, const data_t * curr, data_t * next, const data_t * par, int nx, int nz, data_t dt);