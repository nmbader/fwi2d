#ifdef DOUBLE_PRECISION
    typedef double data_t;
#else
    typedef float data_t;
#endif

export void Dz_ispc(uniform int add, uniform data_t pmod[], uniform data_t  pdat[], uniform data_t scale, uniform int i1, uniform int imin, uniform int imax){

    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    
    foreach ( i=imin ... imax){
         pdat[ i1+ i] = add*pdat[ i1+ i] + scale * (coef0 * (pmod[ i1+ i+1] - pmod[ i1+ i-1]) + coef1 * (pmod[ i1+ i+2] - pmod[ i1+ i-2]));
    }
}

export void Dx_ispc(uniform int add, uniform data_t pmod[], uniform data_t  pdat[], uniform data_t scale, uniform int i1, uniform int n, uniform int imin, uniform int imax){

    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    
    foreach ( i=imin ... imax){
         pdat[ i1*n+ i] = add*pdat[ i1*n+ i] + scale * (coef0 * (pmod[(i1+1)*n+i] - pmod[(i1-1)*n+i]) + coef1 * (pmod[(i1+2)*n+i] - pmod[(i1-2)*n+i]));
    }
}

export void mult_Dz_ispc(uniform int add, uniform data_t pmod[], uniform data_t  pdat[], uniform data_t scale, uniform int i1, uniform int imin, uniform int imax, uniform data_t  par[]){

    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    
    foreach ( i=imin ... imax){
         pdat[ i1+ i] = add*pdat[ i1+ i] + scale * (coef0 * (pmod[i1+i+1]*par[i1+i+1] - pmod[i1+i-1]*par[i1+i-1]) + coef1 * (pmod[i1+i+2]*par[i1+i+2] - pmod[i1+i-2]*par[i1+i-2]));
    }
}

export void mult_Dx_ispc(uniform int add, uniform data_t pmod[], uniform data_t  pdat[], uniform data_t scale, uniform int i1, uniform int n, uniform int imin, uniform int imax, uniform data_t  par[]){

    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    
    foreach ( i=imin ... imax){
         pdat[ i1*n+ i] = add*pdat[ i1*n+ i] + scale * (coef0 * (pmod[(i1+1)*n+i]*par[(i1+1)*n+i] - pmod[(i1-1)*n+i]*par[(i1-1)*n+i]) + coef1 * (pmod[(i1+2)*n+i]*par[(i1+2)*n+i] - pmod[(i1-2)*n+i]*par[(i1-2)*n+i]));
    }
}

export void taperz_ispc(uniform data_t  pdat[], uniform int i1, uniform int istart, uniform int iend, uniform data_t a){ 
    
    data_t val;
    if (iend>istart){
        foreach ( i=istart ... iend){
            val = cos(a*0.5*PI*(i-istart)/(iend-istart));
            pdat[ i1+ i] *= val*val;
        }
    }
    else{
        foreach ( i=iend ... istart){
            val = cos(a*0.5*PI*(i+1-istart)/(iend-istart));
            pdat[ i1+ i] *= val*val;
        }
    }
}

export void taperx_ispc(uniform data_t pdat[], uniform int i1,  uniform int n, uniform int izmin, uniform int izmax, uniform data_t a){ 
    
    foreach ( i=izmin ... izmax){
        pdat[i1+i] = a*a*pdat[i1+i];
    }
}