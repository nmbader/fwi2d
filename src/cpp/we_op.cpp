#include "we_op.hpp"
#include "misc.hpp"
#include "spatial_operators.hpp"
#include "IO.hpp"

void analyzeGeometry(const hypercube<data_t> &model, param &par, bool verbose)
{
    successCheck(model.getNdim()==3,__FILE__,__LINE__,"The model must have three axes\n");
    par.nmodels = model.getAxis(3).n;
    successCheck(par.nmodels==2 || par.nmodels==3 || par.nmodels==5,__FILE__,__LINE__,"The model must contain 2, 3, or 5 parameters\n");
    if (par.acoustic_elastic) successCheck(par.nmodels==3,__FILE__,__LINE__,"For acoustic-elastic option, the model must contain 3 parameters\n");

    if (par.nmodels==2){
        par.mt=false;
        par.gl=0;
        par.free_surface_stiffness=std::max(par.free_surface_stiffness,(data_t)1.0);
        par.soft_clip=0;
        par.model_parameterization=1;
        par.pml=false;
    }
    if (par.acoustic_elastic && par.acoustic_source) par.mt=false;

    if (verbose) fprintf(stderr,"\n==========================\n Subsurface model geometry\n==========================\n");
    axis<data_t> Z = model.getAxis(1);
    axis<data_t> X = model.getAxis(2);
    if (verbose){
        fprintf(stderr,"xmin=%.5f km, xmax=%.5f km, dx=%0.5f km, nx=%d\n",X.o,(X.n-1)*X.d+X.o,X.d,X.n);
        fprintf(stderr,"zmin=%.5f km, zmax=%.5f km, dz=%0.5f km, nz=%d\n",Z.o,(Z.n-1)*Z.d+Z.o,Z.d,Z.n);
    }
    if (par.acoustic_elastic) {
        successCheck(Z.o>0 && std::abs(Z.o/Z.d-round(Z.o/Z.d))<Z.d/100.,__FILE__,__LINE__,"For acoustic-elastic option, the z-origin must be positive multiple of the sampling\n");
        par.nza = floor(Z.o/Z.d) + 1;
        if (verbose) fprintf(stderr,"A fluid layer is implicitly added on top of the model, all the way up to z=0 km with nz=%d samples v=%.2f km/s and rho=%.2f g/cc\n",par.nza,par.water_velocity,par.water_density);
    }

    if (verbose) fprintf(stderr,"\n==========================\n Boundary conditions\n==========================\n");
    par.pml_T=false; par.pml_B=false; par.pml_L=false; par.pml_R=false;
    if (par.pml==true)
    {
        int max_taper = 0;
        max_taper = std::max(par.taper_bottom, par.taper_top);
        max_taper = std::max(max_taper, par.taper_left);
        max_taper = std::max(max_taper, par.taper_right);
        if (max_taper == 0) par.pml=false;
        else
        {
            if (par.taper_top>0) {par.bc_top=0; par.taper_top=max_taper; par.pml_T=true;} else par.bc_top = std::max(1, par.bc_top); 
            if (par.taper_bottom>0) {par.bc_bottom=0; par.taper_bottom=max_taper; par.pml_B=true;} else par.bc_bottom = std::max(1, par.bc_bottom);
            if (par.taper_left>0) {par.bc_left=0; par.taper_left=max_taper; par.pml_L=true;} else par.bc_left = std::max(1, par.bc_left);
            if (par.taper_right>0) {par.bc_right=0; par.taper_right=max_taper; par.pml_R=true;} else par.bc_right = std::max(1, par.bc_right);
        }
    }
    if (par.pml==false){
        par.bc_top = std::max(1, par.bc_top);
        par.bc_bottom = std::max(1, par.bc_bottom);
        par.bc_left = std::max(1, par.bc_left);
        par.bc_right = std::max(1, par.bc_right);
    }
    std::string bc[6] = {"none","free surface","locally absorbing","rigid-wall","mixed rigid-absorbing","mixed rigid-free"};
    std::string taper_type[3] = {"none","cosine squared","PML"};
    if (verbose){
        fprintf(stderr,"Top boundary condition = %s\t taper size = %d\t taper type = %s\n",bc[par.bc_top].c_str(), par.taper_top, taper_type[(par.taper_top>0)*(1+par.pml)].c_str()); 
        fprintf(stderr,"Bottom boundary condition = %s\t taper size = %d\t taper type = %s\n",bc[par.bc_bottom].c_str(), par.taper_bottom, taper_type[(par.taper_bottom>0)*(1+par.pml)].c_str());
        fprintf(stderr,"Left boundary condition = %s\t taper size = %d\t taper type = %s\n",bc[par.bc_left].c_str(), par.taper_left, taper_type[(par.taper_left>0)*(1+par.pml)].c_str()); 
        fprintf(stderr,"Right boundary condition = %s\t taper size = %d\t taper type = %s\n",bc[par.bc_right].c_str(), par.taper_right, taper_type[(par.taper_right>0)*(1+par.pml)].c_str());
    }

    if (verbose) {
        fprintf(stderr,"\n==========================\n Sources' and receivers' geometry\n==========================\n");
        if (par.srcoord_from_file==true) {
            fprintf(stderr,"Sources' and receivers' coordinates are read from file %s\n",par.srcoord.c_str());
        }
        else {
            fprintf(stderr,"Sources' and receivers' coordinates are read from parameters list\n");
        }
        std::string seismotype1[2] = {"displacement","velocity"};
        std::string seismotype2[2] = {"strain","strain rate"};
        std::string seismotype3[2] = {" "," derivative of "};
        std::string seismotype4[2] = {"momentum potential","hydraulic pressure"};
        if (par.gl<=0) {
            if (par.nmodels>=3 && !par.acoustic_elastic) fprintf(stderr,"Receivers are point measurement of type particle %s\n",seismotype1[par.seismotype].c_str());
            else if (par.acoustic_elastic) fprintf(stderr,"Receivers are dual point measurement of type %s and particle %s\n",seismotype4[par.seismotype].c_str(),seismotype1[par.seismotype].c_str());
            else fprintf(stderr,"Receivers are point measurement of type%shydraulic pressure\n",seismotype3[par.seismotype].c_str());
        }
        else {
            if (!par.acoustic_elastic) fprintf(stderr,"Receivers are DAS measurement of type %s with gauge length = %f km\n",seismotype2[par.seismotype].c_str(), par.gl);
            else fprintf(stderr,"Receivers are dual point measurement of type %s and DAS of type %s with gauge length = %f km\n",seismotype4[par.seismotype].c_str(),seismotype2[par.seismotype].c_str(), par.gl);
        }
    }

    if (verbose) fprintf(stderr,"Number of sources = %d\n",par.ns);
    if (par.acoustic_elastic && par.acoustic_source) fprintf(stderr,"All sources must fall within the fluid layer\n");
    else if (par.acoustic_elastic && !par.acoustic_source) fprintf(stderr,"All sources must fall within the solid layer\n");
    bool check=true;
    int s=0;
    par.nr=0;
    data_t precision = std::min(X.d,Z.d)*1e-3;
    if (par.acoustic_elastic) par.rxza = par.rxz;
    while (s<par.ns)
    {
        int nr = par.rxz[s].size();
        par.nr+=nr;
        // check source inside model boundaries
        data_t oz = Z.o; int nz=Z.n;
        data_t ozs = Z.o; int nzs=Z.n;
        if (par.acoustic_elastic) {
            oz=0; nz=Z.n+par.nza-1;
            if (par.acoustic_source){ozs=0; nzs=par.nza;}
        }
        if (verbose) fprintf(stderr,"Shot %d located at x=%.5f km, z=%.5f km, has %d receivers\n",s,par.sxz[s][0],par.sxz[s][1],nr);
        check = ((par.sxz[s][0] - X.o >= -precision) && (-precision <= -par.sxz[s][0] + X.o+(X.n-1)*X.d) && (par.sxz[s][1] -ozs >= -precision) && (-precision <= -par.sxz[s][1] + ozs+(nzs-1)*Z.d));
        successCheck(check,__FILE__,__LINE__,"Shot falls outside subsurface model boundaries\n");

        // check receivers inside model boundaries
        int r=0;
        while (r<nr)
        {
            data_t rx=par.rxz[s][r][0];
            data_t rz=par.rxz[s][r][1];
            data_t dip=par.rxz[s][r][2];

            // project receivers to fluid-solid interface when relevant
            if (par.acoustic_elastic){
                if (rz<Z.o && rz>=0) {par.rxz[s][r][1]=Z.o; rz=Z.o;}
                else if (rz>=Z.o) par.rxza[s][r][1]=Z.o;
            }
            if (par.gl<=0)
            {
                check = ((rx -X.o >= -precision) && (-precision <= -rx + X.o+(X.n-1)*X.d) && (rz -Z.o >= -precision) && (-precision <= -rz + Z.o+(Z.n-1)*Z.d));
                std::string msg = "Receiver "+std::to_string(r)+" at x="+std::to_string(rx)+" km, z="+std::to_string(rz)+" km, falls outside subsurface model boundaries\n";
                successCheck(check,__FILE__,__LINE__,msg);
            }
            else
            {
                check=
                    !((rx - par.gl/2*cos(dip) < X.o) 
                    || (rx - par.gl/2*cos(dip) > X.o+(X.n-1)*X.d)
                    || (rx + par.gl/2*cos(dip) < X.o)
                    || (rx + par.gl/2*cos(dip) > X.o+(X.n-1)*X.d)
                    || (rz - par.gl/2*sin(dip) < Z.o) 
                    || (rz - par.gl/2*sin(dip) > Z.o+(Z.n-1)*Z.d)
                    || (rz + par.gl/2*sin(dip) < Z.o)
                    || (rz + par.gl/2*sin(dip) > Z.o+(Z.n-1)*Z.d)
                    );
                std::string msg = "DAS channel "+std::to_string(r)+" at x="+std::to_string(rx)+" km, z="+std::to_string(rz)+" km, dip="+std::to_string(dip)+" rad, must be away from the subsurface model boundaries by half of the gauge length along the fiber\n";
                successCheck(check,__FILE__,__LINE__,msg);
            }
            r++;
        }
        s++;
    }
    if (verbose){
        if (par.gl<=0) {
            if (par.nmodels>=3 && !par.acoustic_elastic) fprintf(stderr,"Total number of 2-components receivers to be modeled = %d\n",par.nr);
            else if (par.nmodels==3 && par.acoustic_elastic) fprintf(stderr,"Total number of 3-components receivers to be modeled = %d\n",par.nr);
            else fprintf(stderr,"Total number of hydrophones to be modeled = %d\n",par.nr);
        }
        else {
            if (!par.acoustic_elastic) fprintf(stderr,"Total number of DAS channels to be modeled = %d\n",par.nr);
            else fprintf(stderr,"Total number of hydrophones+DAS channels to be modeled = %d\n",par.nr);
        }
    }
    if (par.gl>0)
    {
        par.rdip.clear();
        for (int s=0; s<par.ns; s++){
            for (int r=0; r<par.rxz[s].size(); r++){
                par.rdip.push_back(par.rxz[s][r][2]);
            }
        }
    }
}

std::shared_ptr<vecReg<data_t> > analyzeWavelet(std::shared_ptr<vecReg<data_t> > src, const param &par, bool verbose)
{
    if (verbose) fprintf(stderr,"\n==========================\n Source wavelet\n==========================\n");
    if (src->getHyper()->getNdim()==1)
    {
        if (par.nmodels==2 || (par.acoustic_elastic && par.acoustic_source))
        {
            if (verbose) {
                fprintf(stderr,"Acoustic source (pressure or momentum potential) is assumed\n");
                fprintf(stderr,"Input wavelet will be duplicated %d times\n",par.ns);
                fprintf(stderr,"Parameters fangle, mt, mxx, mzz, mxz are ignored\n");
            }
            axis<data_t> T = src->getHyper()->getAxis(1);
            axis<data_t> S(par.ns,0,1);
            axis<data_t> C(1,0,1);
            hypercube<data_t> hyp(T,S,C);
            std::shared_ptr<vecReg<data_t> > allsrc = std::make_shared<vecReg<data_t> > (hyp);
            const data_t * p = src->getCVals();
            data_t (* pall) [T.n] = (data_t (*) [T.n]) allsrc->getVals();
            for (int its=0; its<S.n; its++) memcpy(pall[its], p, T.n*sizeof(data_t));
            return allsrc;
        }
        else if (par.mt==false)
        {
            if (verbose) {
                fprintf(stderr,"Vector force is assumed with an angle of %f rad clockwise wrt to the horizontal\n",par.fangle);
                fprintf(stderr,"Input wavelet will be duplicated 2 x %d times with the appropriate scaling\n",par.ns);
                fprintf(stderr,"Parameters mxx, mzz, mxz are ignored\n");
            }
            axis<data_t> T = src->getHyper()->getAxis(1);
            axis<data_t> S(par.ns,0,1);
            axis<data_t> C(2,0,1);
            hypercube<data_t> hyp(T,S,C);
            std::shared_ptr<vecReg<data_t> > allsrc = std::make_shared<vecReg<data_t> > (hyp);
            const data_t * p = src->getCVals();
            data_t (* pall) [S.n][T.n] = (data_t (*) [S.n][T.n]) allsrc->getVals();
            for (int its=0; its<S.n; its++){
                for (int it=0; it<T.n; it++){
                    pall[0][its][it] = p[it]*cos(par.fangle); // fx component
                    pall[1][its][it] = p[it]*sin(par.fangle); // fz component
                }
            }
            return allsrc;
        } 
        else 
        {
            if (verbose){
                fprintf(stderr,"Moment tensor is assumed with mxx=%f, mzz=%f, mxz=%f\n",par.mxx,par.mzz,par.mxz);
                fprintf(stderr,"Input wavelet will be duplicated 3 x %d times with the appropriate scaling\n",par.ns);
                fprintf(stderr,"Parameter fangle is ignored\n");
            }
            axis<data_t> T = src->getHyper()->getAxis(1);
            axis<data_t> S(par.ns,0,1);
            axis<data_t> C(3,0,1);
            hypercube<data_t> hyp(T,S,C);
            std::shared_ptr<vecReg<data_t> > allsrc = std::make_shared<vecReg<data_t> > (hyp);
            const data_t * p = src->getCVals();
            data_t (* pall) [S.n][T.n] = (data_t (*) [S.n][T.n]) allsrc->getVals();
            for (int its=0; its<S.n; its++){
                for (int it=0; it<T.n; it++){
                    pall[0][its][it] = p[it]*par.mxx; // mxx component
                    pall[1][its][it] = p[it]*par.mzz; // mzz component
                    pall[2][its][it] = p[it]*par.mxz; // mxz component
                }
            }
            return allsrc;
        }
    }
    else
    {
        int ntr = src->getN123()/src->getHyper()->getAxis(1).n;
        if (par.nmodels==2 || (par.acoustic_elastic && par.acoustic_source))
        {
            std::string msg = std::to_string(par.ns) + " wavelets are expected, " + std::to_string(ntr) + " are provided\n";
            successCheck(ntr == par.ns,__FILE__,__LINE__,msg);
            if (verbose) fprintf(stderr,"Parameters fangle, mt, mxx, mzz, mxz are ignored\n");
        }
        else if (par.mt==false)
        {
            std::string msg = std::to_string(2*par.ns) + " wavelets are expected, " + std::to_string(ntr) + " are provided\n";
            successCheck(ntr == 2*par.ns,__FILE__,__LINE__,msg);
            if (verbose) fprintf(stderr,"Parameters mxx, mzz, mxz are ignored\n");
        }
        else
        {
            std::string msg = std::to_string(3*par.ns) + " wavelets are expected, " + std::to_string(ntr) + " are provided\n";
            successCheck(ntr == 3*par.ns,__FILE__,__LINE__,msg);
            if (verbose) fprintf(stderr,"Parameter fangle is ignored\n");
        }
        return src;
    }
}

void analyzeModel(const hypercube<data_t> &domain, std::shared_ptr<vecReg<data_t> > model, param &par)
{
    successCheck(model->getHyper()->getNdim()==3,__FILE__,__LINE__,"The model must have three axes\n");
    successCheck(par.nmodels==2 || par.nmodels==3 || par.nmodels==5,__FILE__,__LINE__,"The model must contain two (Vp, rho), three (Vp, Vs, rho) or five (Vp, Vs, rho, delta, eps) parameters\n");

    axis<data_t> Z = model->getHyper()->getAxis(1);
    axis<data_t> X = model->getHyper()->getAxis(2);
    int nx=X.n;
    int nz=Z.n;

    if (par.verbose>1) {
        fprintf(stderr,"\n==========================\n Subsurface model bounds\n==========================\n");
        if (par.nmodels==2) fprintf(stderr,"Model is assumed to contain Vp (km/s) and density Rho (g/cc) in that order\n");
        else if (par.nmodels==3) fprintf(stderr,"Model is assumed to contain Vp (km/s), Vs (km/s) and density Rho (g/cc) in that order\n");
        else if (par.nmodels==5) fprintf(stderr,"Model is assumed to contain Vp (km/s), Vs (km/s), density Rho (g/cc), Delta, and Epsilon in that order\n");
    }

    if (par.nmodels>=3)
    {
        data_t vpmin = model->min(0,nx*nz);
        data_t vpmax = model->max(0,nx*nz);
        data_t vsmin = model->min(nx*nz,2*nx*nz);
        data_t vsmax = model->max(nx*nz,2*nx*nz);
        data_t rhomin = model->min(2*nx*nz,3*nx*nz);
        data_t rhomax = model->max(2*nx*nz,3*nx*nz);
        data_t delmin, delmax, epsmin, epsmax;

        if (par.nmodels==5)
        {
            delmin = model->min(3*nx*nz,4*nx*nz);
            delmax = model->max(3*nx*nz,4*nx*nz);
            epsmin = model->min(4*nx*nz,5*nx*nz);
            epsmax = model->max(4*nx*nz,5*nx*nz);
        }

        if (par.verbose>1) {
            fprintf(stderr,"Vp bounds are %.2f - %.2f km/s\n",vpmin,vpmax);
            fprintf(stderr,"Vs bounds are %.2f - %.2f km/s\n",vsmin,vsmax);
            fprintf(stderr,"Rho bounds are %.2f - %.2f g/cc\n",rhomin,rhomax);
            if (par.nmodels==5){
                fprintf(stderr,"Delta bounds are %.2f - %.2f\n",delmin,delmax);
                fprintf(stderr,"Epsilon bounds are %.2f - %.2f\n",epsmin,epsmax);
            }
        }

        if (!par.soft_clip)
        {
            model->clip(par.vpmin,par.vpmax,0,nx*nz);
            model->clip(par.vsmin,par.vsmax,nx*nz,2*nx*nz);
            model->clip(par.rhomin,par.rhomax,2*nx*nz,3*nx*nz);

            if (par.nmodels==5)
            {
                model->clip(par.deltamin,par.deltamax,3*nx*nz,4*nx*nz);
                model->clip(par.epsilonmin,par.epsilonmax,4*nx*nz,5*nx*nz);
                delmin = model->min(3*nx*nz,4*nx*nz);
                delmax = model->max(3*nx*nz,4*nx*nz);
                epsmin = model->min(4*nx*nz,5*nx*nz);
                epsmax = model->max(4*nx*nz,5*nx*nz);
            }

            bool check=true;
            data_t (* pm) [nx][nz]= (data_t (*) [nx][nz]) model->getVals();
            for (int ix=0; ix<nx; ix++){
                for (int iz=0; iz<nz; iz++){
                    if (pm[0][ix][iz] < sqrt(2)*pm[1][ix][iz]){
                        check = false;
                        pm[1][ix][iz] = pm[0][ix][iz] / sqrt(2.0001);
                    }
                }
            }
            //successCheck(check,__FILE__,__LINE__,"Vp values must always be larger than sqrt(2)*Vs\n");
            if (par.verbose>1 && !check) fprintf(stderr,"WARNING: Vs exceeds Vp/sqrt(2) at some locations and will be clipped accordingly\n");

            vpmin = model->min(0,nx*nz);
            vpmax = model->max(0,nx*nz);
            vsmin = model->min(nx*nz,2*nx*nz);
            vsmax = model->max(nx*nz,2*nx*nz);
            rhomin = model->min(2*nx*nz,3*nx*nz);
            rhomax = model->max(2*nx*nz,3*nx*nz);

            if (par.verbose>1) {
                fprintf(stderr,"Vp bounds after hard clipping are %.2f - %.2f km/s\n",vpmin,vpmax);
                fprintf(stderr,"Vs bounds after hard clipping are %.2f - %.2f km/s\n",vsmin,vsmax);
                fprintf(stderr,"Rho bounds after hard clipping are %.2f - %.2f g/cc\n",rhomin,rhomax);
                if (par.nmodels==5){
                    fprintf(stderr,"Delta bounds after hard clipping are %.2f - %.2f\n",delmin,delmax);
                    fprintf(stderr,"Epsilon bounds after hard clipping are %.2f - %.2f\n",epsmin,epsmax);
                }
            }
        }

        par.vmax=vpmax;
        par.vmin=vsmin;

        // compute the maximum horizontal P-wave velocity
        if (par.nmodels==5) {
            std::shared_ptr<vecReg<data_t> > temp = std::make_shared<vecReg<data_t> >(hypercube<data_t>(Z,X));
            data_t * ptemp = temp->getVals();
            const data_t * pvp = model->getCVals();
            const data_t * pe = model->getCVals()+4*nx*nz;
            #pragma omp parallel for
            for (int i=0; i<nx*nz; i++) ptemp[i] = pvp[i]*sqrt(1+2*pe[i]);
            
            par.vmax = temp->max();
        }

        if (par.acoustic_elastic){
            par.vmax = std::max(par.vmax, par.water_velocity);
            par.vmin = std::min(par.vmin, par.water_velocity);
        }
    }

    else
    {
        data_t vpmin = model->min(0,nx*nz);
        data_t vpmax = model->max(0,nx*nz);
        data_t rhomin = model->min(nx*nz,2*nx*nz);
        data_t rhomax = model->max(nx*nz,2*nx*nz);

        if (par.verbose>1) {
            fprintf(stderr,"Vp bounds are %.2f - %.2f km/s\n",vpmin,vpmax);
            fprintf(stderr,"Rho bounds are %.2f - %.2f g/cc\n",rhomin,rhomax);
        }

        model->clip(par.vpmin,par.vpmax,0,nx*nz);
        model->clip(par.rhomin,par.rhomax,nx*nz,2*nx*nz);

        vpmin = model->min(0,nx*nz);
        vpmax = model->max(0,nx*nz);
        rhomin = model->min(nx*nz,2*nx*nz);
        rhomax = model->max(nx*nz,2*nx*nz);

        if (par.verbose>1) {
            fprintf(stderr,"Vp bounds after hard clipping are %.2f - %.2f km/s\n",vpmin,vpmax);
            fprintf(stderr,"Rho bounds after hard clipping are %.2f - %.2f g/cc\n",rhomin,rhomax);
        }

        par.vmax=vpmax;
        par.vmin=vpmin;
    }    

    if (par.verbose>2) {
        fprintf(stderr,"\n==========================\n Dispersion analysis\n==========================\n");
        fprintf(stderr,"Maximum frequency assumed by the user = %.2f Hz\n",par.fmax);
        fprintf(stderr,"Corresponding minimum number of grid points per wavelength according to minimum velocity = %.1f\n",par.vmin/(par.fmax*std::max(X.d,Z.d)));
    }
    
    if (par.verbose>2) fprintf(stderr,"\n==========================\n Time analysis\n==========================\n");
    axis<data_t> T = domain.getAxis(1);
    data_t tmax = (T.n-1)*T.d;
    if (par.verbose>2){
        fprintf(stderr,"Maximum duration of the input = %.3f sec, sampling = %.5f sec, number of samples = %d\n",tmax, T.d, T.n);
        fprintf(stderr,"Original Courant number = %.2f\n",par.courant);
        fprintf(stderr,"Corresponding propagation time step based on CFL condition = %.5f sec\n",par.courant*std::min(X.d,Z.d)/par.vmax);
    }

    if (par.dt==0) par.dt=par.courant*std::min(X.d,Z.d)/par.vmax;
    else if (par.dt > 0) {
        par.dt = std::min(par.dt,T.d);
        par.courant = par.dt * par.vmax / std::min(X.d,Z.d);
    }
    else {
        par.dt=par.courant*std::min(X.d,Z.d)/par.vmax;
        if (T.d > par.dt)
        {
            par.dt = T.d/ceil((T.d/par.dt));
        }
        else 
        {
            par.dt = T.d;
        }
        par.courant = par.dt * par.vmax / std::min(X.d,Z.d);
    }
    if (par.verbose>2) fprintf(stderr,"Updated Courant number = %.2f\n",par.courant);
    par.nt = round(tmax/par.dt) + 1;
    par.tmax=(par.nt-1)*par.dt;
    if (par.verbose>2){
        fprintf(stderr,"Maximum duration of the propagation = %.3f sec, sampling = %.5f sec, number of samples = %d\n", par.tmax, par.dt, par.nt);
        if (par.resampling!="linear") par.resampling="sinc";
        fprintf(stderr,"Time resampling method = %s interpolation\n", par.resampling.c_str());
    }

    if (par.sub < 0) par.sub = round(T.d/par.dt);
    if (par.sub > 0 && par.verbose>2) fprintf(stderr,"Full wavefield (gradient) will be saved (computed) every %d propagation time steps\n", par.sub);
}

void nl_we_op_e::convert_model(data_t * m, int n, bool forward) const
{
    data_t (* pm) [n] = (data_t (*) [n]) m;

    if (forward)
    {
        for (int i=0; i<n; i++){
            pm[0][i] = pm[2][i]*(pm[0][i]*pm[0][i]-2*pm[1][i]*pm[1][i]);
            pm[1][i] = pm[2][i]*pm[1][i]*pm[1][i];
        }
    }
    else
    {
        for (int i=0; i<n; i++){
            pm[0][i] = sqrt((pm[0][i]+2*pm[1][i])/pm[2][i]);
            pm[1][i] = sqrt(pm[1][i]/pm[2][i]);
        }
    }
}

// variable coefficients expressions for 2nd derivative SBP operators
static inline data_t expr1(const data_t ** par, int i){return par[0][i];} // e.g. = lambda
static inline data_t expr2(const data_t ** par, int i){return par[1][i];} // e.g. = mu
static inline data_t expr3(const data_t ** par, int i){return par[2][i];} // e.g. = rho
static inline data_t expr4a(const data_t ** par, int i){return par[0][i] + 2*par[1][i];} // e.g. = lambda + 2 mu
static inline data_t expr4b(const data_t ** par, int i){return 2*par[1][i];} // e.g. = 2 mu
static inline data_t expr5(const data_t ** par, int i){return 0.5*sqrt(par[1][i]*par[2][i]);} // e.g. = 1/2 sqrt(rho.mu)
static inline data_t expr6(const data_t ** par, int i){return 0.5*sqrt((par[0][i]+2*par[1][i])*par[2][i]);} // e.g. = 1/2 sqrt(rho.(lambda+2 mu))

void nl_we_op_e::compute_gradients(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const
{
    // Grad_lambda for wide stencil = div(adjoint).H.Ht.div(forward) = (adjointx_x + adjointz_z).H.Ht.(forwardx_x + forwardz_z)
    // Grad_mu for wide stencil = (adjointx_z+adjointz_x).H.Ht.(forwardx_z+forwardz_x) + 2.adjointx_x.H.Ht.forwardx_x + 2.adjointz_z.H.Ht.forwardz_z
    // Grad_rho = adjointx.H.Ht.forwardx_tt + adjointz.H.Ht.forwardz_tt
    // The H quadrature will be applied elsewhere to the final gradients (all shots included)

    int nxz = nx*nz;
    data_t (*pm) [nxz] = (data_t (*)[nxz]) model;
    data_t (*pfor) [2][nxz] = (data_t (*) [2][nxz]) u_full;
    data_t (*padj) [nxz] = (data_t (*)[nxz]) curr;
    data_t (*padj_x) [nxz] = (data_t (*)[nxz]) u_x;
    data_t (*padj_z) [nxz] = (data_t (*)[nxz]) u_z;
    data_t (*temp) [nxz] = (data_t (*)[nxz]) tmp;
    data_t (*g) [nxz] = (data_t (*)[nxz]) grad;

    Dx(false, pfor[it][0], temp[0], nx, nz, dx, 0, nx, 0, nz); // forwardx_x
    Dz(false, pfor[it][1], temp[1], nx, nz, dz, 0, nx, 0, nz); // forwardz_z
    Dx(false, pfor[it][1], temp[2], nx, nz, dx, 0, nx, 0, nz); // forwardz_x
    Dz(false, pfor[it][0], temp[3], nx, nz, dz, 0, nx, 0, nz); // forwardx_z

    // cosine square damping in a given time window, often around the source function (gradient close to the source)
    data_t t = it*dt;
    data_t w = 1;
    if (par.gdamp_tmax-t>0 && t-par.gdamp_tmin>=0) w = cos(par.gdamp_power*M_PI*(t-par.gdamp_tmin)/(par.gdamp_tmax-par.gdamp_tmin));
    w *= w;

    // different from first and last time sample
    if (it>0 && it<par.nt/par.sub){
        #pragma omp parallel for
        for (int i=0; i<nxz; i++) 
        {
            g[0][i] += w*dt*(padj_x[0][i] + padj_z[1][i])*(temp[0][i] + temp[1][i]); // lambda gradient
            g[1][i] += w*dt*((padj_z[0][i] + padj_x[1][i])*(temp[2][i] + temp[3][i]) + 2*padj_x[0][i]*temp[0][i] + 2*padj_z[1][i]*temp[1][i]); // mu gradient
            g[2][i] += w*1.0/dt*(padj[0][i]*(pfor[it+1][0][i]-2*pfor[it][0][i]+pfor[it-1][0][i]) + padj[1][i]*(pfor[it+1][1][i]-2*pfor[it][1][i]+pfor[it-1][1][i])); // rho gradient
        }
    }
    // time zero
    else if(it==0){
        #pragma omp parallel for
        for (int i=0; i<nxz; i++) 
        {
            g[0][i] += w*0.5*dt*(padj_x[0][i] + padj_z[1][i])*(temp[0][i] + temp[1][i]); // lambda gradient
            g[1][i] += w*0.5*dt*((padj_z[0][i] + padj_x[1][i])*(temp[2][i] + temp[3][i]) + 2*padj_x[0][i]*temp[0][i] + 2*padj_z[1][i]*temp[1][i]); // mu gradient
            g[2][i] += w*1.0/dt*(padj[0][i]*(pfor[it+1][0][i]-pfor[it][0][i]) + padj[1][i]*(pfor[it+1][1][i]-pfor[it][1][i])); // rho gradient
        }
    }
    // tmax
    else{
        #pragma omp parallel for
        for (int i=0; i<nxz; i++) 
        {
            g[0][i] += w*dt*(padj_x[0][i] + padj_z[1][i])*(temp[0][i] + temp[1][i]); // lambda gradient
            g[1][i] += w*dt*((padj_z[0][i] + padj_x[1][i])*(temp[2][i] + temp[3][i]) + 2*padj_x[0][i]*temp[0][i] + 2*padj_z[1][i]*temp[1][i]); // mu gradient
            g[2][i] += w*1.0/dt*(padj[0][i]*(-pfor[it][0][i]+pfor[it-1][0][i]) + padj[1][i]*(-pfor[it][1][i]+pfor[it-1][1][i])); // rho gradient
        }
    }
}

void nl_we_op_e::propagate(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const
{
    // wavefields allocations and pointers
    data_t * u  = new data_t[6*nx*nz];
    data_t * dux  = new data_t[2*nx*nz];
    data_t * duz  = new data_t[2*nx*nz];
    memset(u, 0, 6*nx*nz*sizeof(data_t));
    memset(dux, 0, 2*nx*nz*sizeof(data_t));
    memset(duz, 0, 2*nx*nz*sizeof(data_t));
    data_t (*prev) [nx*nz] = (data_t (*)[nx*nz]) u;
    data_t (*curr) [nx*nz] = (data_t (*)[nx*nz]) (u + 2*nx*nz);
    data_t (*next) [nx*nz] = (data_t (*)[nx*nz]) (u + 4*nx*nz);
    data_t (*bucket) [nx*nz];
    data_t (*u_x) [nx*nz] = (data_t (*)[nx*nz]) dux;
    data_t (*u_z) [nx*nz] = (data_t (*)[nx*nz]) duz;
    data_t (*u_full) [2][nx*nz];
    data_t * tmp;

    if (par.sub>0) u_full = (data_t (*) [2][nx*nz]) full_wfld;
    if (grad != nullptr) 
    {
        tmp = new data_t[4*nx*nz];
        memset(tmp, 0, 4*nx*nz*sizeof(data_t));
    }
    else if (par.version==2)
    {
        tmp = new data_t[nx*nz];
        memset(tmp, 0, nx*nz*sizeof(data_t));
    }
	const data_t* mod[3] = {model, model+nx*nz, model+2*nx*nz};

    // #############  PML stuff ##################
    int l = std::max(par.taper_top,par.taper_bottom);
    l = std::max(l, par.taper_left);
    l = std::max(l, par.taper_right);
    data_t * gx; data_t * gz; data_t * gxprime; data_t * gzprime;
    data_t * top_ac; data_t * top_an; data_t * top_s2;
    data_t * bottom_ac; data_t * bottom_an; data_t * bottom_s2;
    data_t * left_ac; data_t * left_an; data_t * left_s2; data_t * left_t2c; data_t * left_t2n; data_t * left_t3c; data_t * left_t3n; data_t * left_s5; data_t * left_s6;
    data_t * right_ac; data_t * right_an; data_t * right_s2; data_t * right_t2c; data_t * right_t2n; data_t * right_t3c; data_t * right_t3n; data_t * right_s5; data_t * right_s6;
    
    if (par.pml) 
    {
        pml_allocate(top_ac, top_an, top_s2, bottom_ac, bottom_an, bottom_s2,
        left_ac, left_an, left_s2, left_t2c, left_t2n, left_t3c, left_t3n, left_s5, left_s6,
        right_ac, right_an, right_s2, right_t2c, right_t2n, right_t3c, right_t3n, right_s5, right_s6,
        par.pml_T, par.pml_B, par.pml_L, par.pml_R, nx, nz, l);

        // stretching functions g = (p+1)/(2L)*vmax*log(1/R)(d/L)^p where L is the PML thickness, d=x or z, p is the polynomial power (e.g. =2)
        gx = new data_t [l];
        gz = new data_t [l];
        gxprime = new data_t [l];
        gzprime = new data_t [l];
        data_t Lx=l*dx; data_t Lz=dz*l;
        data_t gx0=(par.p+1)/(2*Lx)*par.vmax*log(1.0/par.R);
        data_t gz0=(par.p+1)/(2*Lz)*par.vmax*log(1.0/par.R);
        for (int i=0; i<l; i++)
        {
            gx[i]=gx0*pow((1.0*(l-i)/l),par.p); gxprime[i]=par.p/Lx*gx0*pow((1.0*(l-i)/l),par.p-1);
            gz[i]=gz0*pow((1.0*(l-i)/l),par.p); gzprime[i]=par.p/Lz*gz0*pow((1.0*(l-i)/l),par.p-1);
        }
    }
    // ###########################################
    
    // source and receiver components for injection/extraction
    int ns=par.ns;
    int nr=par.nr;
    if (adj)
    {
        ns = par.nr;
        nr = par.ns;
    }
    const data_t * srcx[2] = {allsrc, nullptr};
    const data_t * srcz[2] = {allsrc + ns*par.nt, nullptr};
    data_t * rcvx[2] = {allrcv, nullptr};
    data_t * rcvz[2] = {allrcv + nr*par.nt, nullptr};
    if (par.mt==true)
    {
        if (!adj)
        {
            srcx[1] = allsrc + 2*ns*par.nt;
            srcz[0] = allsrc + 2*ns*par.nt;
            srcz[1] = allsrc + ns*par.nt;
        }
        else
        {
            rcvx[1] = allrcv + 2*nr*par.nt;
            rcvz[0] = allrcv + 2*nr*par.nt;
            rcvz[1] = allrcv + nr*par.nt;
        }
    }

    // prev = 1/(2 * rho) * dt2 * src
    inj->inject(false, srcx, prev[0], nx, nz, par.nt, inj->_ntr, 0, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
    inj->inject(false, srcz, prev[1], nx, nz, par.nt, inj->_ntr, 0, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
    for (int i=0; i<nx*nz; i++)
    {
        prev[0][i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
        prev[1][i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
    }    

    int pct10 = round(par.nt/10);

    for (int it=0; it<par.nt-1; it++)
    {
        // copy the current wfld to the full wfld vector
        if ((par.sub>0) && (it%par.sub==0) && (grad==nullptr)) memcpy(u_full[it/par.sub], curr, 2*nx*nz*sizeof(data_t));

        // extract receivers
        if (grad == nullptr)
        {
            ext->extract(true, curr[0], rcvx, nx, nz, par.nt, ext->_ntr, it, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
            ext->extract(true, curr[1], rcvz, nx, nz, par.nt, ext->_ntr, it, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
        }

        // compute FWI gradients except for first and last time samples
        if ((grad != nullptr) && ((par.nt-1-it)%par.sub==0) && it!=0) compute_gradients(model, full_wfld, curr[0], u_x[0], u_z[0], tmp, grad, par, nx, nz, (par.nt-1-it)/par.sub, dx, dz, par.sub*par.dt);

        // apply spatial SBP operators
        if (par.version==1)
        {
            Dxx_var<expr4a>(false, curr[0], next[0], nx, nz, dx, 0, nx, 0, nz, mod, 1.0);
            Dzz_var<expr4a>(false, curr[1], next[1], nx, nz, dz, 0, nx, 0, nz, mod, 1.0);
        }
        else
        {
            Dxx_var<expr4b>(false, curr[0], next[0], nx, nz, dx, 0, nx, 0, nz, mod, 1.0);
            Dzz_var<expr4b>(false, curr[1], next[1], nx, nz, dz, 0, nx, 0, nz, mod, 1.0);
        }
        
        Dxx_var<expr2>(true, curr[1], next[1], nx, nz, dx, 0, nx, 0, nz, mod, 1.0);
        Dzz_var<expr2>(true, curr[0], next[0], nx, nz, dz, 0, nx, 0, nz, mod, 1.0);
        
        Dx(false, curr[0], u_x[0], nx, nz, dx, 0, nx, 0, nz);
        Dz(false, curr[1], u_z[1], nx, nz, dz, 0, nx, 0, nz);
        
        if (par.version==2)
        {
            #pragma omp parallel for
            for (int i=0; i<nx*nz; i++)
            {
                tmp[i] = u_x[0][i] + u_z[1][i];
            }
            mult_Dx(true, tmp, next[0], nx, nz, dx, 0, nx, 0, nz, mod[0], 1.0);
            mult_Dz(true, tmp, next[1], nx, nz, dz, 0, nx, 0, nz, mod[0], 1.0);
        }
        else
        {
            mult_Dx(true, u_z[1], next[0], nx, nz, dx, 0, nx, 0, nz, mod[0], 1.0);
            mult_Dz(true, u_x[0], next[1], nx, nz, dz, 0, nx, 0, nz, mod[0], 1.0);
        }

        Dx(false, curr[1], u_x[1], nx, nz, dx, 0, nx, 0, nz);
        Dz(false, curr[0], u_z[0], nx, nz, dz, 0, nx, 0, nz);
        
        mult_Dx(true, u_z[0], next[1], nx, nz, dx, 0, nx, 0, nz, mod[1], 1.0);
        mult_Dz(true, u_x[1], next[0], nx, nz, dz, 0, nx, 0, nz, mod[1], 1.0);

        // inject sources
        inj->inject(true, srcx, next[0], nx, nz, par.nt, inj->_ntr, it, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
        inj->inject(true, srcz, next[1], nx, nz, par.nt, inj->_ntr, it, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());

        // apply boundary conditions
        if (par.bc_top==1)
        {
            const data_t * in[2] = {curr[1],curr[0]};
            esat_neumann_top<expr2,expr2>(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            in[0] = curr[0]; in[1] = curr[1];
            if (par.version==1) esat_neumann_top<expr1,expr4a>(true, in, next[1], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            else
            {
                esat_neumann_top<expr1,expr4b>(true, in, next[1], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
                esat_Dz_top<expr1>(true, curr[1], next[1], nx, nz, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
        }
        else if (par.bc_top==2)
        {
            const data_t * in[3] = {curr[1],curr[0],prev[0]};
            esat_absorbing_top<expr2,expr2,expr5>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = prev[1];
            if (par.version==1) esat_absorbing_top<expr1,expr4a,expr6>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            else
            {
                esat_absorbing_top<expr1,expr4b,expr6>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
                esat_Dz_top<expr1>(true, curr[1], next[1], nx, nz, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
        }
        if (par.bc_bottom==1)
        {
            const data_t * in[2] = {curr[1],curr[0]};
            esat_neumann_bottom<expr2,expr2>(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            in[0] = curr[0]; in[1] = curr[1];
            if (par.version==1) esat_neumann_bottom<expr1,expr4a>(true, in, next[1], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            else
            {
                esat_neumann_bottom<expr1,expr4b>(true, in, next[1], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
                esat_Dz_bottom<expr1>(true, curr[1], next[1], nx, nz, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
        }
        else if (par.bc_bottom==2)
        {
            const data_t * in[3] = {curr[1],curr[0],prev[0]};
            esat_absorbing_bottom<expr2,expr2,expr5>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = prev[1];
            if (par.version==1) esat_absorbing_bottom<expr1,expr4a,expr6>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            else
            {
                esat_absorbing_bottom<expr1,expr4b,expr6>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
                esat_Dz_bottom<expr1>(true, curr[1], next[1], nx, nz, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
        }
        if (par.bc_left==1)
        {
            const data_t * in[2] = {curr[1],curr[0]};
            if (par.version==1) esat_neumann_left<expr1,expr4a>(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            else
            {
                esat_neumann_left<expr1,expr4b>(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
                esat_Dx_left<expr1>(true, curr[0], next[0], nx, nz, dx, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            }
            in[0] = curr[0]; in[1] = curr[1];
            esat_neumann_left<expr2,expr2>(true, in, next[1], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        else if (par.bc_left==2)
        {
            const data_t * in[3] = {curr[1],curr[0],prev[0]};
            if (par.version==1) esat_absorbing_left<expr1,expr4a,expr6>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            else
            {
                esat_absorbing_left<expr1,expr4b,expr6>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
                esat_Dx_left<expr1>(true, curr[0], next[0], nx, nz, dx, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            }
            in[0] = curr[0]; in[1] = curr[1]; in[2] = prev[1];
            esat_absorbing_left<expr2,expr2,expr5>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        if (par.bc_right==1)
        {
            const data_t * in[2] = {curr[1],curr[0]};
            if (par.version==1) esat_neumann_right<expr1,expr4a>(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            else
            {
                esat_neumann_right<expr1,expr4b>(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
                esat_Dx_right<expr1>(true, curr[0], next[0], nx, nz, dx, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            }
            in[0] = curr[0]; in[1] = curr[1];
            esat_neumann_right<expr2,expr2>(true, in, next[1], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        else if (par.bc_right==2)
        {
            const data_t * in[3] = {curr[1],curr[0],prev[0]};
            if (par.version==1) esat_absorbing_right<expr1,expr4a,expr6>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            else
            {
                esat_absorbing_right<expr1,expr4b,expr6>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
                esat_Dx_right<expr1>(true, curr[0], next[0], nx, nz, dx, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            }
            in[0] = curr[0]; in[1] = curr[1]; in[2] = prev[1];
            esat_absorbing_right<expr2,expr2,expr5>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }

        // update wfld with the 2 steps time recursion
        #pragma omp parallel for
        for (int i=0; i<nx*nz; i++)
        {
            next[0][i] = par.dt*par.dt*next[0][i]/mod[2][i] + 2*curr[0][i] - prev[0][i];
            next[1][i] = par.dt*par.dt*next[1][i]/mod[2][i] + 2*curr[1][i] - prev[1][i];
        }

        // scale boundaries when relevant (for locally absorbing BC only)
        data_t * in[2] = {next[0],next[1]};
        esat_scale_boundaries(in, nx, nz, dx, dz, 0, nx, 0, nz, mod, par.dt, par.bc_top==2, par.bc_bottom==2, par.bc_left==2, par.bc_right==2);
        
        // apply PML or taper when relevant
        if (par.pml==true)
        {
            if (par.pml_T) pml_top(curr[0], next[0], dux, duz, top_ac, top_an, top_s2, model, gz, gzprime, nx, nz, l, par.pml_L*l, nx-par.pml_R*l, dx, dz, par.dt,par.version);
            if (par.pml_B) pml_bottom(curr[0], next[0], dux, duz, bottom_ac, bottom_an, bottom_s2, model, gz, gzprime, nx, nz, l, par.pml_L*l, nx-par.pml_R*l, dx, dz, par.dt,par.version);
            if (par.pml_L) pml_left(curr[0], next[0], dux, duz, left_ac, left_an, left_s2, left_t2c, left_t2n, left_t3c, left_t3n, left_s5, left_s6, true, true, model, gx, gxprime, nx, nz, l, 0, nz, dx, dz, par.dt,par.version);
            if (par.pml_R) pml_right(curr[0], next[0], dux, duz, right_ac, right_an, right_s2, right_t2c, right_t2n, right_t3c, right_t3n, right_s5, right_s6, true, true, model, gx, gxprime, nx, nz, l, 0, nz, dx, dz, par.dt,par.version);
        }
        else
        {
            taperz(curr[0], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
            taperz(curr[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
            taperx(curr[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
            taperx(curr[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
            taperz(curr[1], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
            taperz(curr[1], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
            taperx(curr[1], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
            taperx(curr[1], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
            taperz(next[0], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
            taperz(next[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
            taperx(next[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
            taperx(next[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
            taperz(next[1], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
            taperz(next[1], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
            taperx(next[1], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
            taperx(next[1], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
        }

        bucket=prev;
        prev=curr;
        curr=next;
        next=bucket;

        if ((it+1) % pct10 == 0 && par.verbose>2) fprintf(stderr,"Propagation progress = %d\%\n",100*(it+1)/(par.nt-1));
        else if (it == par.nt-2 && par.verbose>2) fprintf(stderr,"Propagation progress = 100\%\n");
    }

    // copy the last wfld to the full wfld vector
    if ((par.sub>0) && ((par.nt-1)%par.sub==0) && (grad==nullptr)) memcpy(u_full[(par.nt-1)/par.sub], curr, 2*nx*nz*sizeof(data_t));

    // extract receivers last sample
    if (grad == nullptr)
    {
        ext->extract(true, curr[0], rcvx, nx, nz, par.nt, ext->_ntr, par.nt-1, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
        ext->extract(true, curr[1], rcvz, nx, nz, par.nt, ext->_ntr, par.nt-1, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
    }

    // last sample gradient
    if ( grad != nullptr ) compute_gradients(model, full_wfld, curr[0], u_x[0], u_z[0], tmp, grad, par, nx, nz, 0, dx, dz, par.sub*par.dt);
    
    delete [] u;
    delete [] dux;
    delete [] duz;
    if (grad != nullptr || par.version==2) delete [] tmp;
    if (par.pml==true)
    {
        pml_deallocate(top_ac, top_an, top_s2, bottom_ac, bottom_an, bottom_s2,
        left_ac, left_an, left_s2, left_t2c, left_t2n, left_t3c, left_t3n, left_s5, left_s6,
        right_ac, right_an, right_s2, right_t2c, right_t2n, right_t3c, right_t3n, right_s5, right_s6,
        par.pml_T, par.pml_B, par.pml_L, par.pml_R, nx, nz, l);

        delete [] gx;
        delete [] gz;
        delete [] gxprime;
        delete [] gzprime;
    }
}

void nl_we_op_e::apply_forward(bool add, const data_t * pmod, data_t * pdat)
{
//    if (_par.verbose>1) fprintf(stderr,"Start forward propagation\n");

    std::shared_ptr<vecReg<data_t> > model = std::make_shared<vecReg<data_t> >(_domain);
    memcpy(model->getVals(), pmod, _domain.getN123()*sizeof(data_t));
    analyzeModel(*_allsrc->getHyper(),model,_par);
    convert_model(model->getVals(), model->getN123()/_par.nmodels, true);
    const data_t * pm = model->getCVals();

    int nx = _domain.getAxis(2).n;
    int nz = _domain.getAxis(1).n;
    data_t dx = _domain.getAxis(2).d;
    data_t dz = _domain.getAxis(1).d;
    data_t ox = _domain.getAxis(2).o;
    data_t oz = _domain.getAxis(1).o;
    hypercube<data_t> domain = *_allsrc->getHyper();

    if (_par.sub>0){
        if (_par.nmodels>=3) _full_wfld=std::make_shared<vecReg<data_t> > (hypercube<data_t>(_domain.getAxis(1),_domain.getAxis(2),axis<data_t>(2,0,1), axis<data_t>(1+_par.nt/_par.sub,0,_par.dt*_par.sub), axis<data_t>(_par.ns,0,1)));
        else _full_wfld=std::make_shared<vecReg<data_t> > (hypercube<data_t>(_domain.getAxis(1),_domain.getAxis(2), axis<data_t>(1+_par.nt/_par.sub,0,_par.dt*_par.sub), axis<data_t>(_par.ns,0,1)));
        if (_par.acoustic_elastic && _par.acoustic_wavefield) {
            _full_wflda=std::make_shared<vecReg<data_t> > (hypercube<data_t>(axis<data_t>(_par.nza,0,dz),_domain.getAxis(2), axis<data_t>(1+_par.nt/_par.sub,0,_par.dt*_par.sub)));
        }
    }

    // setting up and apply the time resampling operator
    resampler * resamp;

    axis<data_t> T = domain.getAxis(1);
    axis<data_t> Xs = domain.getAxis(2);
    axis<data_t> Xr0 = _range.getAxis(2);
    axis<data_t> Cr0 = _range.getAxis(3);
    axis<data_t> Xr(2*Xr0.n,0,1);
    int ncomp=2;
    if (_par.nmodels==2) {Xr.n = Xr0.n; ncomp=1;}
    if (_par.acoustic_elastic) {Xr.n = 3*Xr0.n; ncomp=3;}
    Xs.n = domain.getN123()/T.n;
    data_t alpha = _par.dt/T.d; // ratio between sampling rates
    T.n = _par.nt;
    T.d = _par.dt;
    hypercube<data_t> hyper_s(T,Xs);
    hypercube<data_t> hyper_r(T,Xr);
    hypercube<data_t> hyper_r0(T,Xr0,Cr0);
    std::shared_ptr<vecReg<data_t> > allsrc = std::make_shared<vecReg<data_t> >(hyper_s);
    std::shared_ptr<vecReg<data_t> > allrcv = std::make_shared<vecReg<data_t> >(hyper_r);
    allsrc->zero();
    allrcv->zero();

    if (_par.resampling == "linear") resamp = new linear_resampler(domain, hyper_s);
    else resamp = new sinc_resampler(domain, hyper_s, _par.sinc_half_length);

    // resample in time (interpolation)
    //fprintf(stderr,"Resample in time all source traces\n");
    resamp->apply_forward(false,_allsrc->getCVals(),allsrc->getVals());

    int size=1, rank=0, rank0=0;
#ifdef ENABLE_MPI
    if (!_par.skip_mpi){
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    }
    MPI_Comm_rank(MPI_COMM_WORLD,&rank0);
#endif

    // loop over shots
    for (int s=rank; s<_par.ns; s += size)
    {
        if (_par.verbose>1) fprintf(stderr,"Start propagating shot %d by process %d\n",s,rank0);

        // cumulative number of receivers
        int nr=0;
        if (s>0) {for (int i=0; i<s; i++) nr += _par.rxz[i].size();}

        // setting up the injection and extraction operators
        injector * inj;
        if (_par.mt==true) inj = new ddelta_m3(_domain,{_par.sxz[s]});
        else inj = new delta_m3(_domain,{_par.sxz[s]});

        injector * ext;
        if (_par.gl<=0) ext = new delta_m3(_domain,_par.rxz[s]);
        else ext = new dipole_m3(_domain,_par.rxz[s],_par.gl);

        // perform the wave propagation
        data_t * full = nullptr;
        if (_par.sub>0) {
            if (_par.nmodels>=3) full = _full_wfld->getVals() + s*nx*nz*2*(1+_par.nt/_par.sub);
            else full = _full_wfld->getVals() + s*nx*nz*(1+_par.nt/_par.sub);
        }

#ifdef CUDA
        propagate_gpu(false, pm, allsrc->getCVals()+s*_par.nt, allrcv->getVals()+nr*_par.nt, inj, ext, full, nullptr, _par, nx, nz, dx, dz, ox, oz, s);
#else
        propagate(false, pm, allsrc->getCVals()+s*_par.nt, allrcv->getVals()+nr*_par.nt, inj, ext, full, nullptr, _par, nx, nz, dx, dz, ox, oz, s);
#endif

        if (_par.verbose>1) fprintf(stderr,"Finish propagating shot %d by process %d\n",s, rank0);

        delete inj;
        delete ext;
    }

#ifdef ENABLE_MPI
    if (!_par.skip_mpi && size>1){
        data_t * temp;
        if (rank==0) temp = new data_t[allrcv->getN123()];
        if (sizeof(data_t)==8) MPI_Reduce(allrcv->getCVals(), temp, allrcv->getN123(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        else MPI_Reduce(allrcv->getCVals(), temp, allrcv->getN123(), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank==0) {
            memcpy(allrcv->getVals(), temp, allrcv->getN123()*sizeof(data_t));
            delete [] temp;
        }

        if (_par.sub>0){
            data_t * temp;
            if (rank==0) temp = new data_t[_full_wfld->getN123()];
            if (sizeof(data_t)==8) MPI_Reduce(_full_wfld->getCVals(), temp, _full_wfld->getN123(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            else MPI_Reduce(_full_wfld->getCVals(), temp, _full_wfld->getN123(), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            if (rank==0) {
                memcpy(_full_wfld->getVals(), temp, _full_wfld->getN123()*sizeof(data_t));
                delete [] temp;
            }
        }
    }
#endif              

    // convert to DAS (strain) data when relevant
    if (_par.gl>0) 
    {
        //fprintf(stderr,"Convert 2-components dipole data to strain\n");
        if (!_par.acoustic_elastic) dipole_to_strain(false, allrcv->getVals(), _par.rdip.data(), Xr0.n, T.n, 0, Xr0.n, _par.gl);
        else dipole_to_strain(false, allrcv->getVals()+T.n*Xr0.n, _par.rdip.data(), Xr0.n, T.n, 0, Xr0.n, _par.gl);
    }

    // convert to particle velocity or strain rate or pressure when relevant: forward  mode
    if (_par.seismotype==1)
    {
        std::shared_ptr<vecReg<data_t> > allrcv_dt = std::make_shared<vecReg<data_t> >(hyper_r);
        allrcv_dt->zero();
        Dt(false, false, allrcv->getCVals(), allrcv_dt->getVals(), Xr.n, T.n, T.d, 0, Xr.n);
        allrcv = allrcv_dt;
        if (_par.acoustic_elastic) allrcv->scale(-1,0,Xr0.n*T.n);
    }

    // resample in time (decimation)
    //fprintf(stderr,"Resample in time all receiver traces\n");
    allrcv->scale(alpha);
    resamp->setDomainRange(_range, hyper_r0);
    resamp->apply_adjoint(add, pdat, allrcv->getCVals());
        
    delete resamp;
}

void nl_we_op_e::apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat)
{
//    if (_par.verbose>1) fprintf(stderr,"Start adjoint propagation\n");

    std::shared_ptr<vecReg<data_t> > model0 = std::make_shared<vecReg<data_t> >(_domain);
    memcpy(model0->getVals(), pmod0, _domain.getN123()*sizeof(data_t));
    int verbose = _par.verbose;
    _par.verbose=0;
    analyzeModel(*_allsrc->getHyper(),model0,_par);
    _par.verbose=verbose;
    convert_model(model0->getVals(), model0->getN123()/_par.nmodels, true);
    const data_t * pm0 = model0->getCVals();

    int nx = _domain.getAxis(2).n;
    int nz = _domain.getAxis(1).n;
    data_t dx = _domain.getAxis(2).d;
    data_t dz = _domain.getAxis(1).d;
    data_t ox = _domain.getAxis(2).o;
    data_t oz = _domain.getAxis(1).o;
    hypercube<data_t> domain = *_allsrc->getHyper();

    data_t * temp;
    int nm=std::min(3,_par.nmodels);
    if (!add) memset(pmod, 0, nm*nx*nz*sizeof(data_t));
    else { // copy pre-existant gradient to a temporary container
        temp = new data_t[nm*nx*nz];
        memcpy(temp, pmod, nm*nx*nz*sizeof(data_t));
    }
 
    // setting up and apply the time resampling operator
    resampler * resamp;

    axis<data_t> T = domain.getAxis(1);
    axis<data_t> Xs = domain.getAxis(2);
    axis<data_t> Xr0 = _range.getAxis(2);
    axis<data_t> Cr0 = _range.getAxis(3);
    axis<data_t> Xr(2*Xr0.n,0,1);
    if (_par.nmodels==2) Xr.n = Xr0.n;
    if (_par.acoustic_elastic) Xr.n = 3*Xr0.n;
    Xs.n = domain.getN123()/T.n;
    data_t alpha = _par.dt/T.d; // ratio between sampling rates
    T.n = _par.nt;
    T.d = _par.dt;
    hypercube<data_t> hyper_s(T,Xs);
    hypercube<data_t> hyper_r(T,Xr);
    hypercube<data_t> hyper_r0(T,Xr0,Cr0);
    std::shared_ptr<vecReg<data_t> > allsrc = std::make_shared<vecReg<data_t> >(hyper_s);
    std::shared_ptr<vecReg<data_t> > allrcv = std::make_shared<vecReg<data_t> >(hyper_r);
    allsrc->zero();
    allrcv->zero();

    if (_par.resampling == "linear") resamp = new linear_resampler(_range, hyper_r0);
    else resamp = new sinc_resampler(_range, hyper_r0, _par.sinc_half_length);

    // resample in time (interpolation)
    resamp->apply_forward(false, pdat, allrcv->getVals());

    // apply inverse time quadrature, revert in time and multiply by -1
    applyHt(true, false, allrcv->getCVals(), allrcv->getVals(), Xr.n, T.n, T.d, 0, Xr.n);
    allrcv->revert(1);
    allrcv->scale(-alpha);

    // convert from particle velocity or strain rate or pressure when relevant: adjoint mode
    if (_par.seismotype==1)
    {
        std::shared_ptr<vecReg<data_t> > allrcv_dt = std::make_shared<vecReg<data_t> >(hyper_r);
        allrcv_dt->zero();
        Dt(true, false, allrcv->getCVals(), allrcv_dt->getVals(), Xr.n, T.n, T.d, 0, Xr.n);
        allrcv = allrcv_dt;
        if (_par.acoustic_elastic) allrcv->scale(-1,0,Xr0.n*T.n);
    }

    // convert from DAS (strain) to dipole data when relevant
    if (_par.gl>0) 
    {
        if (!_par.acoustic_elastic) dipole_to_strain(true, allrcv->getVals(), _par.rdip.data(), Xr0.n, T.n, 0, Xr0.n, _par.gl);
        else dipole_to_strain(true, allrcv->getVals()+T.n*Xr0.n, _par.rdip.data(), Xr0.n, T.n, 0, Xr0.n, _par.gl);
    }

    int size=1, rank=0, rank0=0;
#ifdef ENABLE_MPI
    if (!_par.skip_mpi){
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    }
    MPI_Comm_rank(MPI_COMM_WORLD,&rank0);
#endif

    // loop over shots
    for (int s=rank; s<_par.ns; s += size)
    {
        if (_par.verbose>1) fprintf(stderr,"Start back-propagating shot %d by process %d\n",s,rank0);

        // cumulative number of receivers
        int nr=0;
        if (s>0) {for (int i=0; i<s; i++) nr += _par.rxz[i].size();}

        // setting up the injection and extraction operators
        injector * ext;
        if (_par.mt==true) ext = new ddelta_m3(_domain,{_par.sxz[s]});
        else ext = new delta_m3(_domain,{_par.sxz[s]});

        injector * inj;
        if (_par.gl<=0) inj = new delta_m3(_domain,_par.rxz[s]);
        else inj = new dipole_m3(_domain,_par.rxz[s],_par.gl);

        // perform the wave propagation
        data_t * full = nullptr;
        if (_par.sub>0) {
            if (_par.nmodels>=3) full = _full_wfld->getVals() + s*nx*nz*2*(1+_par.nt/_par.sub);
            else full = _full_wfld->getVals() + s*nx*nz*(1+_par.nt/_par.sub);
        }

#ifdef CUDA
        propagate_gpu(true, pm0, allrcv->getCVals()+nr*_par.nt, allsrc->getVals()+s*_par.nt, inj, ext, full, pmod, _par, nx, nz, dx, dz, ox, oz, s);
#else
        propagate(true, pm0, allrcv->getCVals()+nr*_par.nt, allsrc->getVals()+s*_par.nt, inj, ext, full, pmod, _par, nx, nz, dx, dz, ox, oz, s);
#endif

        if (_par.verbose>1) fprintf(stderr,"Finish back-propagating shot %d with gradient computation by process %d\n",s, rank0);

        delete inj;
        delete ext;
    }

#ifdef ENABLE_MPI
    if (!_par.skip_mpi && size>1){
        data_t * gtemp = new data_t[nm*nx*nz];
        if (sizeof(data_t)==8) MPI_Allreduce(pmod, gtemp, nm*nx*nz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        else MPI_Allreduce(pmod, gtemp, nm*nx*nz, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        memcpy(pmod, gtemp, nm*nx*nz*sizeof(data_t));
        delete [] gtemp;
    }
#endif

    // apply H quadrature to final gradient
    int nxz=nx*nz;
    applyHz(false, false, pmod, pmod, nx, nz, dz, 0, nx, 0, nz);
    applyHz(false, false, pmod+nxz, pmod+nxz, nx, nz, dz, 0, nx, 0, nz);
    if (nm==3) applyHz(false, false, pmod+2*nxz, pmod+2*nxz, nx, nz, dz, 0, nx, 0, nz);
    applyHx(false, false, pmod, pmod, nx, nz, dx, 0, nx, 0, nz);
    applyHx(false, false, pmod+nxz, pmod+nxz, nx, nz, dx, 0, nx, 0, nz);
    if (nm==3) applyHx(false, false, pmod+2*nxz, pmod+2*nxz, nx, nz, dx, 0, nx, 0, nz);

    // convert gradients from lambda, mu, rho to vp, vs, rho    or   from K-1, rho-1 to vp, rho
    // Grad_vp = 2.rho.vp.Grad_lambda = 2.sqrt(rho(lambda+2mu)).Grad_lambda
    // Grad_vs = 2.rho.vs.(Grad_mu - 2.Grad_lambda) = 2.sqrt(rho.mu).(Grad_mu - 2.Grad_lambda)
    // Grad_rho = (vp2 - 2.vs2).Grad_lambda + vs2.Grad_mu + Grad_rho = (lambda/rho).Grad_lambda + (mu/rho).Grad_mu + Grad_rho
    // or
    // Grad_vp = 2.rho.vp.Grad_K = 2.sqrt(K/rho-1).Grad_K
    // Grad_rho = vp^2.Grad_K - 1/rho^2.Grad_rho-1 = K.rho-1.Grad_K - (rho-1)^2.Grad_rho-1
    data_t (*g) [nxz] = (data_t (*)[nxz]) pmod;
    data_t (*pm) [nxz] = (data_t (*)[nxz]) pm0;
    if (nm>=3)
    {
        #pragma omp parallel for
        for(int i=0; i<nxz; i++)
        {
            g[2][i] = (pm[0][i]/pm[2][i])*g[0][i] + (pm[1][i]/pm[2][i])*g[1][i] + g[2][i]; // Rho Gradient
            g[1][i] = 2*sqrt(pm[1][i]*pm[2][i])*(g[1][i] - 2*g[0][i]); // Vs gradient
            g[0][i] = 2*sqrt(pm[2][i]*(pm[0][i]+2*pm[1][i]))*g[0][i]; // Vp gradient
        }
    }
    else
    {
        #pragma omp parallel for
        for(int i=0; i<nxz; i++)
        {
            g[1][i] = pm[1][i]*pm[0][i]*g[0][i] - pm[1][i]*pm[1][i]*g[1][i]; // rho gradient
            g[0][i] = 2*sqrt(pm[0][i]/pm[1][i])*g[0][i]; // Vp gradient
        }
    }

    if(add)
    {
        #pragma omp parallel for
        for(int i=0; i<nm*nxz; i++) pmod[i] += temp[i];

        delete [] temp;
    }

    delete resamp;
}

void l_we_op_e::apply_forward(bool add, const data_t * pmod, data_t * pdat)
{
    if (_par.verbose>1) fprintf(stderr,"Start forward propagation\n");

    int nx = _model->getHyper()->getAxis(2).n;
    int nz = _model->getHyper()->getAxis(1).n;
    data_t dx = _model->getHyper()->getAxis(2).d;
    data_t dz = _model->getHyper()->getAxis(1).d;
    data_t ox = _model->getHyper()->getAxis(2).o;
    data_t oz = _model->getHyper()->getAxis(1).o;

    // setting up and apply the time resampling operator
    resampler * resamp;

    axis<data_t> T = _domain.getAxis(1);
    axis<data_t> Xs = _domain.getAxis(2);
    axis<data_t> Xr0 = _range.getAxis(2);
    axis<data_t> Cr0 = _range.getAxis(3);
    axis<data_t> Xr(2*Xr0.n,0,1);
    if (_par.nmodels==2) Xr.n = Xr0.n;
    if (_par.acoustic_elastic) Xr.n = 3*Xr0.n;
    Xs.n = _domain.getN123()/T.n;
    data_t alpha = _par.dt/T.d; // ratio between sampling rates
    T.n = _par.nt;
    T.d = _par.dt;
    hypercube<data_t> hyper_s(T,Xs);
    hypercube<data_t> hyper_r(T,Xr);
    hypercube<data_t> hyper_r0(T,Xr0,Cr0);
    std::shared_ptr<vecReg<data_t> > allsrc = std::make_shared<vecReg<data_t> >(hyper_s);
    std::shared_ptr<vecReg<data_t> > allrcv = std::make_shared<vecReg<data_t> >(hyper_r);
    allsrc->zero();
    allrcv->zero();

    if (_par.resampling == "linear") resamp = new linear_resampler(_domain, hyper_s);
    else resamp = new sinc_resampler(_domain, hyper_s, _par.sinc_half_length);

    // resample in time (interpolation)
    resamp->apply_forward(false,pmod,allsrc->getVals());

    // loop over shots
    for (int s=0; s<_par.ns; s++)
    {
        if (_par.verbose>1) fprintf(stderr,"Start propagating shot %d\n",s);

        // cumulative number of receivers
        int nr=0;
        if (s>0) {for (int i=0; i<s; i++) nr += _par.rxz[i].size();}

        // setting up the injection and extraction operators
        injector * inj;
        if (_par.mt==true) inj = new ddelta_m3(*_model->getHyper(),{_par.sxz[s]});
        else inj = new delta_m3(*_model->getHyper(),{_par.sxz[s]});

        injector * ext;
        if (_par.gl<=0) ext = new delta_m3(*_model->getHyper(),_par.rxz[s]);
        else ext = new dipole_m3(*_model->getHyper(),_par.rxz[s],_par.gl);

        // perform the wave propagation
        data_t * full = nullptr;
        if (_par.sub>0) full = _full_wfld->getVals();
#ifdef CUDA
        propagate_gpu(false, _model->getCVals(), allsrc->getCVals()+s*_par.nt, allrcv->getVals()+nr*_par.nt, inj, ext, full, nullptr, _par, nx, nz, dx, dz, ox, oz, s);
#else
        propagate(false, _model->getCVals(), allsrc->getCVals()+s*_par.nt, allrcv->getVals()+nr*_par.nt, inj, ext, full, nullptr, _par, nx, nz, dx, dz, ox, oz, s);
#endif
        delete inj;
        delete ext;

        if (_par.verbose>1) fprintf(stderr,"Finish propagating shot %d\n",s);
    }

    // convert to DAS (strain) data when relevant
    if (_par.gl>0) 
    {
        if (!_par.acoustic_elastic) dipole_to_strain(false, allrcv->getVals(), _par.rdip.data(), Xr0.n, T.n, 0, Xr0.n, _par.gl);
        else dipole_to_strain(false, allrcv->getVals()+T.n*Xr0.n, _par.rdip.data(), Xr0.n, T.n, 0, Xr0.n, _par.gl);
    }

    // convert to particle velocity or strain rate or pressure when relevant: forward  mode
    if (_par.seismotype==1)
    {
        std::shared_ptr<vecReg<data_t> > allrcv_dt = std::make_shared<vecReg<data_t> >(hyper_r);
        allrcv_dt->zero();
        Dt(false, false, allrcv->getCVals(), allrcv_dt->getVals(), Xr.n, T.n, T.d, 0, Xr.n);
        allrcv = allrcv_dt;
        if (_par.acoustic_elastic) allrcv->scale(-1,0,Xr0.n*T.n);
    }

    // resample in time (decimation)
    allrcv->scale(alpha);
    resamp->setDomainRange(_range, hyper_r0);
    resamp->apply_adjoint(add, pdat, allrcv->getCVals());
    
    delete resamp;
}

void l_we_op_e::apply_adjoint(bool add, data_t * pmod, const data_t * pdat)
{
    if (_par.verbose>1) fprintf(stderr,"Start adjoint propagation\n");

    int nx = _model->getHyper()->getAxis(2).n;
    int nz = _model->getHyper()->getAxis(1).n;
    data_t dx = _model->getHyper()->getAxis(2).d;
    data_t dz = _model->getHyper()->getAxis(1).d;
    data_t ox = _model->getHyper()->getAxis(2).o;
    data_t oz = _model->getHyper()->getAxis(1).o;

    // setting up and apply the time resampling operator
    resampler * resamp;

    axis<data_t> T = _domain.getAxis(1);
    axis<data_t> Xs = _domain.getAxis(2);
    axis<data_t> Xr0 = _range.getAxis(2);
    axis<data_t> Cr0 = _range.getAxis(3);
    axis<data_t> Xr(2*Xr0.n,0,1);
    if (_par.nmodels==2) Xr.n = Xr0.n;
    if (_par.acoustic_elastic) Xr.n = 3*Xr0.n;
    Xs.n = _domain.getN123()/T.n;
    data_t alpha = _par.dt/T.d; // ratio between sampling rates
    T.n = _par.nt;
    T.d = _par.dt;
    hypercube<data_t> hyper_s(T,Xs);
    hypercube<data_t> hyper_r(T,Xr);
    hypercube<data_t> hyper_r0(T,Xr0,Cr0);
    std::shared_ptr<vecReg<data_t> > allsrc = std::make_shared<vecReg<data_t> >(hyper_s);
    std::shared_ptr<vecReg<data_t> > allrcv = std::make_shared<vecReg<data_t> >(hyper_r);
    allsrc->zero();
    allrcv->zero();

    if (_par.resampling == "linear") resamp = new linear_resampler(_range, hyper_r0);
    else resamp = new sinc_resampler(_range, hyper_r0, _par.sinc_half_length);

    // resample in time (interpolation)
    resamp->apply_forward(false, pdat, allrcv->getVals());

    // apply inverse time quadrature and revert in time
    applyHt(true, false, allrcv->getCVals(), allrcv->getVals(), Xr.n, T.n, T.d, 0, Xr.n);
    allrcv->revert(1);
    allrcv->scale(alpha);

    // convert from particle velocity or strain rate or pressure when relevant: adjoint mode
    if (_par.seismotype==1)
    {
        std::shared_ptr<vecReg<data_t> > allrcv_dt = std::make_shared<vecReg<data_t> >(hyper_r);
        allrcv_dt->zero();
        Dt(true, false, allrcv->getCVals(), allrcv_dt->getVals(), Xr.n, T.n, T.d, 0, Xr.n);
        allrcv = allrcv_dt;
        if (_par.acoustic_elastic) allrcv->scale(-1,0,Xr0.n*T.n);
    }

    // convert from DAS (strain) to dipole data when relevant
    if (_par.gl>0) 
    {
        if (!_par.acoustic_elastic) dipole_to_strain(true, allrcv->getVals(), _par.rdip.data(), Xr0.n, T.n, 0, Xr0.n, _par.gl);
        else dipole_to_strain(true, allrcv->getVals()+T.n*Xr0.n, _par.rdip.data(), Xr0.n, T.n, 0, Xr0.n, _par.gl);
    }

    // loop over shots
    for (int s=0; s<_par.ns; s++)
    {
        if (_par.verbose>1) fprintf(stderr,"Start back-propagating shot %d\n",s);

        // cumulative number of receivers
        int nr=0;
        if (s>0) {for (int i=0; i<s; i++) nr += _par.rxz[i].size();}
        
        // setting up the injection and extraction operators
        injector * ext;
        if (_par.mt==true) ext = new ddelta_m3(*_model->getHyper(),{_par.sxz[s]});
        else ext = new delta_m3(*_model->getHyper(),{_par.sxz[s]});

        injector * inj;
        if (_par.gl<=0) inj = new delta_m3(*_model->getHyper(),_par.rxz[s]);
        else inj = new dipole_m3(*_model->getHyper(),_par.rxz[s],_par.gl);

        // perform the wave propagation
        data_t * full = nullptr;
        if (_par.sub>0) full = _full_wfld->getVals();
#ifdef CUDA
        propagate_gpu(true, _model->getCVals(), allrcv->getCVals()+nr*_par.nt, allsrc->getVals()+s*_par.nt, inj, ext, full, nullptr, _par, nx, nz, dx, dz, ox, oz, s);
#else
        propagate(true, _model->getCVals(), allrcv->getCVals()+nr*_par.nt, allsrc->getVals()+s*_par.nt, inj, ext, full, nullptr, _par, nx, nz, dx, dz, ox, oz, s);
#endif
        delete inj;
        delete ext;

        if (_par.verbose>1) fprintf(stderr,"Finish back-propagating shot %d\n",s);
    }

    // apply time quadrature and revert back in time
    applyHt(false, false, allsrc->getCVals(), allsrc->getVals(), Xs.n, T.n, T.d, 0, Xs.n);
    allsrc->revert(1);

    // resample in time (decimation)
    resamp->setDomainRange(_domain, hyper_s);
    resamp->apply_adjoint(add, pmod, allsrc->getCVals());
    
    delete resamp;
}

void nl_we_op_vti::convert_model(data_t * m, int n, bool forward) const
{
    data_t (* pm) [n] = (data_t (*) [n]) m;
    data_t val=0;

    if (forward)
    {
        for (int i=0; i<n; i++){
            pm[0][i] = pm[2][i]*(pm[0][i]*pm[0][i]-2*pm[1][i]*pm[1][i]); // lambda
            pm[1][i] = pm[2][i]*pm[1][i]*pm[1][i]; // mu
            pm[3][i] = sqrt(2*(pm[0][i]+2*pm[1][i])*(pm[0][i]+pm[1][i])*pm[3][i] + (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) - pm[1][i]; // c13
        }
    }
    else
    {
        for (int i=0; i<n; i++){
            pm[3][i] = ((pm[3][i]+pm[1][i])*(pm[3][i]+pm[1][i]) - (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) / (2*(pm[0][i]+pm[1][i])*(pm[0][i]+2*pm[1][i])); // delta
            pm[0][i] = sqrt((pm[0][i]+2*pm[1][i])/pm[2][i]); // vp
            pm[1][i] = sqrt(pm[1][i]/pm[2][i]); // vs
        }
    }
}

// variable coefficients expressions for 2nd derivative SBP operators
static inline data_t vtiexpr1(const data_t ** par, int i){return par[0][i];} // e.g. = lambda
static inline data_t vtiexpr2(const data_t ** par, int i){return par[1][i];} // e.g. = mu
static inline data_t vtiexpr3(const data_t ** par, int i){return par[2][i];} // e.g. = rho
static inline data_t vtiexpr4(const data_t ** par, int i){return par[3][i];} // e.g. = c13
static inline data_t vtiexpr5(const data_t ** par, int i){return par[4][i];} // e.g. = eps
static inline data_t vtiexpr6a(const data_t ** par, int i){return par[0][i] + 2*par[1][i];} // e.g. = lambda + 2.mu
static inline data_t vtiexpr6b(const data_t ** par, int i){return 2*par[1][i];} // e.g. = 2 mu
static inline data_t vtiexpr7a(const data_t ** par, int i){return (1 + 2*par[4][i])*(par[0][i] + 2*par[1][i]);} // e.g. = c11 = (1 + 2.eps).(lambda + 2.mu)
static inline data_t vtiexpr7b(const data_t ** par, int i){return (1 + 2*par[4][i])*2*par[1][i];} // e.g. = (1 + 2.eps).2.mu
static inline data_t vtiexpr8(const data_t ** par, int i){return 0.5*sqrt(par[1][i]*par[2][i]);} // e.g. = 1/2 sqrt(rho.mu)
static inline data_t vtiexpr9(const data_t ** par, int i){return 0.5*sqrt((par[0][i]+2*par[1][i])*par[2][i]);} // e.g. = 1/2 sqrt(rho.(lambda+2.mu))
static inline data_t vtiexpr10(const data_t ** par, int i){return 0.5*sqrt((1+2*par[4][i])*(par[0][i]+2*par[1][i])*par[2][i]);} // e.g. = 1/2 sqrt(rho.(1+2.eps).(lambda+2.mu))

void nl_we_op_vti::compute_gradients(const data_t * model, const data_t * u_full, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const
{
    // Grad_lambda for wide stencil = (1+2.epsilon).(adjointx_x).H.Ht.(forwardx_x) + (adjointz_z).H.Ht.(forwardz_z) + d(C13)/d(lambda).(adjointx_x.H.Ht.forwardz_z + adjointz_z.H.Ht.forwardx_x)
    // Grad_mu for wide stencil = (adjointx_z+adjointz_x).H.Ht.(forwardx_z+forwardz_x) + 2.(1+2.epsilon).adjointx_x.H.Ht.forwardx_x + 2.adjointz_z.H.Ht.forwardz_z + d(C13)/d(mu).(adjointx_x.H.Ht.forwardz_z + adjointz_z.H.Ht.forwardx_x)
    // Grad_rho = adjointx.H.Ht.forwardx_tt + adjointz.H.Ht.forwardz_tt
    // d(C13)/d(lambda) = ((1+2.delta).lambda + (1+3.delta).mu)/sqrt(2.(lambda+2.mu).(lambda+mu).delta + (lambda+mu)^2)
    // d(C13)/d(mu) = ((1+3.delta).lambda + (1+4.delta).mu)/sqrt(2.(lambda+2.mu).(lambda+mu).delta + (lambda+mu)^2) - 1
    // The H quadrature will be applied elsewhere to the final gradients (all shots included)

    int nxz = nx*nz;
    data_t (*pm) [nxz] = (data_t (*)[nxz]) model;
    data_t (*pfor) [2][nxz] = (data_t (*) [2][nxz]) u_full;
    data_t (*padj) [nxz] = (data_t (*)[nxz]) curr;
    data_t (*padj_x) [nxz] = (data_t (*)[nxz]) u_x;
    data_t (*padj_z) [nxz] = (data_t (*)[nxz]) u_z;
    data_t (*temp) [nxz] = (data_t (*)[nxz]) tmp;
    data_t (*g) [nxz] = (data_t (*)[nxz]) grad;

    Dx(false, pfor[it][0], temp[0], nx, nz, dx, 0, nx, 0, nz); // forwardx_x
    Dz(false, pfor[it][1], temp[1], nx, nz, dz, 0, nx, 0, nz); // forwardz_z
    Dx(false, pfor[it][1], temp[2], nx, nz, dx, 0, nx, 0, nz); // forwardz_x
    Dz(false, pfor[it][0], temp[3], nx, nz, dz, 0, nx, 0, nz); // forwardx_z

    data_t val1=0, val2=0, val3=0, del=0;

    // cosine square damping in a given time window, often around the source function (gradient close to the source)
    data_t t = it*dt;
    data_t w = 1;
    if (par.gdamp_tmax-t>0 && t-par.gdamp_tmin>=0) w = cos(par.gdamp_power*M_PI*(t-par.gdamp_tmin)/(par.gdamp_tmax-par.gdamp_tmin));
    w *= w;

    // different from first and last time sample
    if (it>0 && it<par.nt/par.sub){
        #pragma omp parallel for private(val1,val2,val3,del)
        for (int i=0; i<nxz; i++) 
        {
            del = ((pm[3][i]+pm[1][i])*(pm[3][i]+pm[1][i]) - (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) / (2*(pm[0][i]+pm[1][i])*(pm[0][i]+2*pm[1][i]));((pm[3][i]+pm[1][i])*(pm[3][i]+pm[1][i]) - (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) / (2*(pm[0][i]+pm[1][i])*(pm[0][i]+2*pm[1][i]));
            val1 = sqrt(2*(pm[0][i]+2*pm[1][i])*(pm[0][i]+pm[1][i])*del + (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i]));
            val2 = ((1+2*del)*pm[0][i] + (1+3*del)*pm[1][i])/val1; // d(C13)/d(lambda)
            val3 = ((1+3*del)*pm[0][i] + (1+4*del)*pm[1][i])/val1 - 1; // d(C13)/d(mu)
            g[0][i] += w*dt*((1+2*pm[4][i])*padj_x[0][i]*temp[0][i] + padj_z[1][i]*temp[1][i] + val2*(padj_x[0][i]*temp[1][i] + padj_z[1][i]*temp[0][i])); // lambda gradient
            g[1][i] += w*dt*((padj_z[0][i] + padj_x[1][i])*(temp[2][i] + temp[3][i]) + 2*(1+2*pm[4][i])*padj_x[0][i]*temp[0][i] + 2*padj_z[1][i]*temp[1][i] + val3*(padj_x[0][i]*temp[1][i] + padj_z[1][i]*temp[0][i])); // mu gradient
            g[2][i] += w*1.0/dt*(padj[0][i]*(pfor[it+1][0][i]-2*pfor[it][0][i]+pfor[it-1][0][i]) + padj[1][i]*(pfor[it+1][1][i]-2*pfor[it][1][i]+pfor[it-1][1][i])); // rho gradient
            g[3][i] = 0;
            g[4][i] = 0;
        }
    }
    // time zero
    else if(it==0){
        #pragma omp parallel for private(val1,val2,val3,del)
        for (int i=0; i<nxz; i++) 
        {
            del = ((pm[3][i]+pm[1][i])*(pm[3][i]+pm[1][i]) - (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) / (2*(pm[0][i]+pm[1][i])*(pm[0][i]+2*pm[1][i]));((pm[3][i]+pm[1][i])*(pm[3][i]+pm[1][i]) - (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) / (2*(pm[0][i]+pm[1][i])*(pm[0][i]+2*pm[1][i]));
            val1 = sqrt(2*(pm[0][i]+2*pm[1][i])*(pm[0][i]+pm[1][i])*del + (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i]));
            val2 = ((1+2*del)*pm[0][i] + (1+3*del)*pm[1][i])/val1; // d(C13)/d(lambda)
            val3 = ((1+3*del)*pm[0][i] + (1+4*del)*pm[1][i])/val1 - 1; // d(C13)/d(mu)
            g[0][i] += w*0.5*dt*((1+2*pm[4][i])*padj_x[0][i]*temp[0][i] + padj_z[1][i]*temp[1][i] + val2*(padj_x[0][i]*temp[1][i] + padj_z[1][i]*temp[0][i])); // lambda gradient
            g[1][i] += w*0.5*dt*((padj_z[0][i] + padj_x[1][i])*(temp[2][i] + temp[3][i]) + 2*(1+2*pm[4][i])*padj_x[0][i]*temp[0][i] + 2*padj_z[1][i]*temp[1][i] + val3*(padj_x[0][i]*temp[1][i] + padj_z[1][i]*temp[0][i])); // mu gradient
            g[2][i] += w*1.0/dt*(padj[0][i]*(pfor[it+1][0][i]-pfor[it][0][i]) + padj[1][i]*(pfor[it+1][1][i]-pfor[it][1][i])); // rho gradient
            g[3][i] = 0;
            g[4][i] = 0;
        }
    }
    // tmax
    else{
        #pragma omp parallel for private(val1,val2,val3,del)
        for (int i=0; i<nxz; i++) 
        {
            del = ((pm[3][i]+pm[1][i])*(pm[3][i]+pm[1][i]) - (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) / (2*(pm[0][i]+pm[1][i])*(pm[0][i]+2*pm[1][i]));((pm[3][i]+pm[1][i])*(pm[3][i]+pm[1][i]) - (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i])) / (2*(pm[0][i]+pm[1][i])*(pm[0][i]+2*pm[1][i]));
            val1 = sqrt(2*(pm[0][i]+2*pm[1][i])*(pm[0][i]+pm[1][i])*del + (pm[0][i]+pm[1][i])*(pm[0][i]+pm[1][i]));
            val2 = ((1+2*del)*pm[0][i] + (1+3*del)*pm[1][i])/val1; // d(C13)/d(lambda)
            val3 = ((1+3*del)*pm[0][i] + (1+4*del)*pm[1][i])/val1 - 1; // d(C13)/d(mu)
            g[0][i] += w*dt*((1+2*pm[4][i])*padj_x[0][i]*temp[0][i] + padj_z[1][i]*temp[1][i] + val2*(padj_x[0][i]*temp[1][i] + padj_z[1][i]*temp[0][i])); // lambda gradient
            g[1][i] += w*dt*((padj_z[0][i] + padj_x[1][i])*(temp[2][i] + temp[3][i]) + 2*(1+2*pm[4][i])*padj_x[0][i]*temp[0][i] + 2*padj_z[1][i]*temp[1][i] + val3*(padj_x[0][i]*temp[1][i] + padj_z[1][i]*temp[0][i])); // mu gradient
            g[2][i] += w*1.0/dt*(padj[0][i]*(-pfor[it][0][i]+pfor[it-1][0][i]) + padj[1][i]*(-pfor[it][1][i]+pfor[it-1][1][i])); // rho gradient
            g[3][i] = 0;
            g[4][i] = 0;
        }
    }
}

void nl_we_op_vti::propagate(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const
{
    // wavefields allocations and pointers
    data_t * u  = new data_t[6*nx*nz];
    data_t * dux  = new data_t[2*nx*nz];
    data_t * duz  = new data_t[2*nx*nz];
    memset(u, 0, 6*nx*nz*sizeof(data_t));
    memset(dux, 0, 2*nx*nz*sizeof(data_t));
    memset(duz, 0, 2*nx*nz*sizeof(data_t));
    data_t (*prev) [nx*nz] = (data_t (*)[nx*nz]) u;
    data_t (*curr) [nx*nz] = (data_t (*)[nx*nz]) (u + 2*nx*nz);
    data_t (*next) [nx*nz] = (data_t (*)[nx*nz]) (u + 4*nx*nz);
    data_t (*bucket) [nx*nz];
    data_t (*u_x) [nx*nz] = (data_t (*)[nx*nz]) dux;
    data_t (*u_z) [nx*nz] = (data_t (*)[nx*nz]) duz;
    data_t (*u_full) [2][nx*nz];
    data_t * tmp;

    if (par.sub>0) u_full = (data_t (*) [2][nx*nz]) full_wfld;
    if (grad != nullptr) 
    {
        tmp = new data_t[4*nx*nz];
        memset(tmp, 0, 4*nx*nz*sizeof(data_t));
    }
    else if (par.version==2)
    {
        tmp = new data_t[nx*nz];
        memset(tmp, 0, nx*nz*sizeof(data_t));
    }
	const data_t* mod[5] = {model, model+nx*nz, model+2*nx*nz, model+3*nx*nz, model+4*nx*nz};
    data_t * mod_bis;
    if (par.version==2)
    {
        mod_bis = new data_t [nx*nz];
        for (int i=0; i<nx*nz; i++) mod_bis[i] = (1+2*mod[4][i])*mod[0][i]; // (1 + 2.eps).lambda
    }
    const data_t * mod6[1] = {mod_bis};

    // #############  PML stuff ##################
    int l = std::max(par.taper_top,par.taper_bottom);
    l = std::max(l, par.taper_left);
    l = std::max(l, par.taper_right);
    // ###########################################

    // source and receiver components for injection/extraction
    int ns=par.ns;
    int nr=par.nr;
    if (adj)
    {
        ns = par.nr;
        nr = par.ns;
    }
    const data_t * srcx[2] = {allsrc, nullptr};
    const data_t * srcz[2] = {allsrc + ns*par.nt, nullptr};
    data_t * rcvx[2] = {allrcv, nullptr};
    data_t * rcvz[2] = {allrcv + nr*par.nt, nullptr};
    if (par.mt==true)
    {
        if (!adj)
        {
            srcx[1] = allsrc + 2*ns*par.nt;
            srcz[0] = allsrc + 2*ns*par.nt;
            srcz[1] = allsrc + ns*par.nt;
        }
        else
        {
            rcvx[1] = allrcv + 2*nr*par.nt;
            rcvz[0] = allrcv + 2*nr*par.nt;
            rcvz[1] = allrcv + nr*par.nt;
        }
    }

    // prev = 1/(2 * rho) * dt2 * src
    inj->inject(false, srcx, prev[0], nx, nz, par.nt, inj->_ntr, 0, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
    inj->inject(false, srcz, prev[1], nx, nz, par.nt, inj->_ntr, 0, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
    for (int i=0; i<nx*nz; i++)
    {
        prev[0][i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
        prev[1][i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
    }    

    int pct10 = round(par.nt/10);

    for (int it=0; it<par.nt-1; it++)
    {
        // copy the current wfld to the full wfld vector
        if ((par.sub>0) && (it%par.sub==0) && (grad==nullptr)) memcpy(u_full[it/par.sub], curr, 2*nx*nz*sizeof(data_t));

        // extract receivers
        if (grad == nullptr)
        {
            ext->extract(true, curr[0], rcvx, nx, nz, par.nt, ext->_ntr, it, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
            ext->extract(true, curr[1], rcvz, nx, nz, par.nt, ext->_ntr, it, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
        }

        // compute FWI gradients except for first and last time samples
        if ((grad != nullptr) && ((par.nt-1-it)%par.sub==0) && it!=0) compute_gradients(model, full_wfld, curr[0], u_x[0], u_z[0], tmp, grad, par, nx, nz, (par.nt-1-it)/par.sub, dx, dz, par.sub*par.dt);
        
        // apply spatial SBP operators
        if (par.version==1)
        {
            Dxx_var<vtiexpr7a>(false, curr[0], next[0], nx, nz, dx, 0, nx, 0, nz, mod, 1.0);
            Dzz_var<vtiexpr6a>(false, curr[1], next[1], nx, nz, dz, 0, nx, 0, nz, mod, 1.0);
        }
        else
        {
            Dxx_var<vtiexpr7b>(false, curr[0], next[0], nx, nz, dx, 0, nx, 0, nz, mod, 1.0);
            Dzz_var<vtiexpr6b>(false, curr[1], next[1], nx, nz, dz, 0, nx, 0, nz, mod, 1.0);
        }
        
        Dxx_var<vtiexpr2>(true, curr[1], next[1], nx, nz, dx, 0, nx, 0, nz, mod, 1.0);
        Dzz_var<vtiexpr2>(true, curr[0], next[0], nx, nz, dz, 0, nx, 0, nz, mod, 1.0);
        
        Dx(false, curr[0], u_x[0], nx, nz, dx, 0, nx, 0, nz);
        Dz(false, curr[1], u_z[1], nx, nz, dz, 0, nx, 0, nz);
        
        if (par.version==2)
        {
            #pragma omp parallel for
            for (int i=0; i<nx*nz; i++)
            {
                tmp[i] = mod_bis[i]*u_x[0][i] + mod[3][i]*u_z[1][i]; // (1+2.eps).lambda*ux_x + c13*uz_z
            }
            Dx(true, tmp, next[0], nx, nz, dx, 0, nx, 0, nz);
            #pragma omp parallel for
            for (int i=0; i<nx*nz; i++)
            {
                tmp[i] = mod[0][i]*u_z[1][i] + mod[3][i]*u_x[0][i]; // lambda*uz_z + c13*ux_x
            }
            Dz(true, tmp, next[1], nx, nz, dz, 0, nx, 0, nz);
        }
        else
        {
            mult_Dx(true, u_z[1], next[0], nx, nz, dx, 0, nx, 0, nz, mod[3], 1.0);
            mult_Dz(true, u_x[0], next[1], nx, nz, dz, 0, nx, 0, nz, mod[3], 1.0);
        }

        Dx(false, curr[1], u_x[1], nx, nz, dx, 0, nx, 0, nz);
        Dz(false, curr[0], u_z[0], nx, nz, dz, 0, nx, 0, nz);
        
        mult_Dx(true, u_z[0], next[1], nx, nz, dx, 0, nx, 0, nz, mod[1], 1.0);
        mult_Dz(true, u_x[1], next[0], nx, nz, dz, 0, nx, 0, nz, mod[1], 1.0);

        // inject sources
        inj->inject(true, srcx, next[0], nx, nz, par.nt, inj->_ntr, it, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
        inj->inject(true, srcz, next[1], nx, nz, par.nt, inj->_ntr, it, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
    
        // apply boundary conditions
        if (par.bc_top==1)
        {
            const data_t * in[2] = {curr[1],curr[0]};
            esat_neumann_top<vtiexpr2,vtiexpr2>(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            in[0] = curr[0]; in[1] = curr[1];
            if (par.version==1) esat_neumann_top<vtiexpr4,vtiexpr6a>(true, in, next[1], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            else
            {
                esat_neumann_top<vtiexpr4,vtiexpr6b>(true, in, next[1], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
                esat_Dz_top<vtiexpr1>(true, curr[1], next[1], nx, nz, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
        }
        else if (par.bc_top==2)
        {
            const data_t * in[3] = {curr[1],curr[0],prev[0]};
            esat_absorbing_top<vtiexpr2,vtiexpr2,vtiexpr8>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = prev[1];
            if (par.version==1) esat_absorbing_top<vtiexpr4,vtiexpr6a,vtiexpr9>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            else
            {
                esat_absorbing_top<vtiexpr4,vtiexpr6b,vtiexpr9>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
                esat_Dz_top<vtiexpr1>(true, curr[1], next[1], nx, nz, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
        }
        if (par.bc_bottom==1)
        {
            const data_t * in[2] = {curr[1],curr[0]};
            esat_neumann_bottom<vtiexpr2,vtiexpr2>(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            in[0] = curr[0]; in[1] = curr[1];
            if (par.version==1) esat_neumann_bottom<vtiexpr4,vtiexpr6a>(true, in, next[1], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            else
            {
                esat_neumann_bottom<vtiexpr4,vtiexpr6b>(true, in, next[1], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
                esat_Dz_bottom<vtiexpr1>(true, curr[1], next[1], nx, nz, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
        }
        else if (par.bc_bottom==2)
        {
            const data_t * in[3] = {curr[1],curr[0],prev[0]};
            esat_absorbing_bottom<vtiexpr2,vtiexpr2,vtiexpr8>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = prev[1];
            if (par.version==1) esat_absorbing_bottom<vtiexpr4,vtiexpr6a,vtiexpr9>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            else
            {
                esat_absorbing_bottom<vtiexpr4,vtiexpr6b,vtiexpr9>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
                esat_Dz_bottom<vtiexpr1>(true, curr[1], next[1], nx, nz, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
        }
        if (par.bc_left==1)
        {
            const data_t * in[2] = {curr[1],curr[0]};
            if (par.version==1) esat_neumann_left<vtiexpr4,vtiexpr7a>(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            else
            {
                esat_neumann_left<vtiexpr4,vtiexpr7b>(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
                esat_Dx_left<vtiexpr1>(true, curr[0], next[0], nx, nz, dx, par.pml_T*l, nz-par.pml_B*l, mod6, 1.0);
            }
            in[0] = curr[0]; in[1] = curr[1];
            esat_neumann_left<vtiexpr2,vtiexpr2>(true, in, next[1], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        else if (par.bc_left==2)
        {
            const data_t * in[3] = {curr[1],curr[0],prev[0]};
            if (par.version==1) esat_absorbing_left<vtiexpr4,vtiexpr7a,vtiexpr10>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            else
            {
                esat_absorbing_left<vtiexpr4,vtiexpr7b,vtiexpr10>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
                esat_Dx_left<vtiexpr1>(true, curr[0], next[0], nx, nz, dx, par.pml_T*l, nz-par.pml_B*l, mod6, 1.0);
            }
            in[0] = curr[0]; in[1] = curr[1]; in[2] = prev[1];
            esat_absorbing_left<vtiexpr2,vtiexpr2,vtiexpr8>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        if (par.bc_right==1)
        {
            const data_t * in[2] = {curr[1],curr[0]};
            if (par.version==1) esat_neumann_right<vtiexpr4,vtiexpr7a>(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            else
            {
                esat_neumann_right<vtiexpr4,vtiexpr7b>(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
                esat_Dx_right<vtiexpr1>(true, curr[0], next[0], nx, nz, dx, par.pml_T*l, nz-par.pml_B*l, mod6, 1.0);
            }
            in[0] = curr[0]; in[1] = curr[1];
            esat_neumann_right<vtiexpr2,vtiexpr2>(true, in, next[1], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        else if (par.bc_right==2)
        {
            const data_t * in[3] = {curr[1],curr[0],prev[0]};
            if (par.version==1) esat_absorbing_right<vtiexpr4,vtiexpr7a,vtiexpr10>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            else
            {
                esat_absorbing_right<vtiexpr4,vtiexpr7b,vtiexpr10>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
                esat_Dx_right<vtiexpr1>(true, curr[0], next[0], nx, nz, dx, par.pml_T*l, nz-par.pml_B*l, mod6, 1.0);
            }
            in[0] = curr[0]; in[1] = curr[1]; in[2] = prev[1];
            esat_absorbing_right<vtiexpr2,vtiexpr2,vtiexpr8>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }

        // update wfld with the 2 steps time recursion
        #pragma omp parallel for
        for (int i=0; i<nx*nz; i++)
        {
            next[0][i] = par.dt*par.dt*next[0][i]/mod[2][i] + 2*curr[0][i] - prev[0][i];
            next[1][i] = par.dt*par.dt*next[1][i]/mod[2][i] + 2*curr[1][i] - prev[1][i];
        }

        // scale boundaries when relevant (for locally absorbing BC only)
        data_t * in[2] = {next[0],next[1]};
        vtisat_scale_boundaries(in, nx, nz, dx, dz, 0, nx, 0, nz, mod, par.dt, par.bc_top==2, par.bc_bottom==2, par.bc_left==2, par.bc_right==2);
        
        // apply taper when relevant
        taperz(curr[0], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
        taperz(curr[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        taperx(curr[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        taperx(curr[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
        taperz(curr[1], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
        taperz(curr[1], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        taperx(curr[1], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        taperx(curr[1], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
        taperz(next[0], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
        taperz(next[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        taperx(next[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        taperx(next[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
        taperz(next[1], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
        taperz(next[1], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        taperx(next[1], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        taperx(next[1], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);

        bucket=prev;
        prev=curr;
        curr=next;
        next=bucket;

        if ((it+1) % pct10 == 0 && par.verbose>2) fprintf(stderr,"Propagation progress = %d\%\n",100*(it+1)/(par.nt-1));
        else if (it == par.nt-2 && par.verbose>2) fprintf(stderr,"Propagation progress = 100\%\n");
    }

    // copy the last wfld to the full wfld vector
    if ((par.sub>0) && ((par.nt-1)%par.sub==0) && (grad==nullptr)) memcpy(u_full[(par.nt-1)/par.sub], curr, 2*nx*nz*sizeof(data_t));

    // extract receivers last sample
    if (grad == nullptr)
    {
        ext->extract(true, curr[0], rcvx, nx, nz, par.nt, ext->_ntr, par.nt-1, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
        ext->extract(true, curr[1], rcvz, nx, nz, par.nt, ext->_ntr, par.nt-1, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
    }

    // last sample gradient
    if (grad != nullptr ) compute_gradients(model, full_wfld, curr[0], u_x[0], u_z[0], tmp, grad, par, nx, nz, 0, dx, dz, par.sub*par.dt);
    
    delete [] u;
    delete [] dux;
    delete [] duz;
    if (par.version==2) delete [] mod_bis;
    if (grad != nullptr || par.version==2) delete [] tmp;
}

void nl_we_op_a::convert_model(data_t * m, int n, bool forward) const
{
    data_t (* pm) [n] = (data_t (*) [n]) m;

    if (forward)
    {
        for (int i=0; i<n; i++){
            pm[0][i] = pm[1][i]*pm[0][i]*pm[0][i];
            pm[1][i] = 1.0/pm[1][i];
        }
    }
    else
    {
        for (int i=0; i<n; i++){
            pm[0][i] = sqrt(pm[0][i]*pm[1][i]);
            pm[1][i] = 1.0/pm[1][i];
        }
    }
}

void nl_we_op_a::compute_gradients(const data_t * model, const data_t * u_full, const data_t * curr, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const
{
    // Grad_K-1 = adjoint.H.Ht.forward_tt
    // Grad_K = -(1/K^2).Grad_K-1
    // Grad_rho-1 = adjoint_x.H.Ht.forward_x + adjoint_z.H.Ht.forward_z
    // The H quadrature will be applied elsewhere to the final gradients (all shots included)

    int nxz = nx*nz;
    data_t (*pm) [nxz] = (data_t (*)[nxz]) model;
    data_t (*pfor) [nxz] = (data_t (*) [nxz]) u_full;
    const data_t *padj = curr;
    data_t (*temp) [nxz] = (data_t (*)[nxz]) tmp;
    data_t (*g) [nxz] = (data_t (*)[nxz]) grad;

    Dx(false, pfor[it], temp[0], nx, nz, dx, 0, nx, 0, nz); // forward_x
    Dz(false, pfor[it], temp[1], nx, nz, dz, 0, nx, 0, nz); // forward_z
    Dx(false, curr, temp[2], nx, nz, dx, 0, nx, 0, nz); // adjoint_x
    Dz(false, curr, temp[3], nx, nz, dz, 0, nx, 0, nz); // adjoint_z

    // cosine square damping in a given time window, often around the source function (gradient close to the source)
    data_t t = it*dt;
    data_t w = 1;
    if (par.gdamp_tmax-t>0 && t-par.gdamp_tmin>=0) w = cos(par.gdamp_power*M_PI*(t-par.gdamp_tmin)/(par.gdamp_tmax-par.gdamp_tmin));
    w *= w;

    // different from first and last time sample
    if (it>0 && it<par.nt/par.sub){
        #pragma omp parallel for
        for (int i=0; i<nxz; i++) 
        {
            g[0][i] -= w*1.0/(dt*pm[0][i]*pm[0][i])*padj[i]*(pfor[it+1][i]-2*pfor[it][i]+pfor[it-1][i]); // K gradient
            g[1][i] += w*dt*(temp[2][i]*temp[0][i] + temp[3][i]*temp[1][i]); // rho-1 gradient
        }
    }
    // time zero
    else if(it==0){
        #pragma omp parallel for
        for (int i=0; i<nxz; i++) 
        {
            g[0][i] -= w*1.0/(dt*pm[0][i]*pm[0][i])*padj[i]*(pfor[it+1][i]-pfor[it][i]); // K gradient
            g[1][i] += w*0.5*dt*(temp[2][i]*temp[0][i] + temp[3][i]*temp[1][i]); // rho-1 gradient
        }
    }
    // tmax
    else{
        #pragma omp parallel for
        for (int i=0; i<nxz; i++) 
        {
            g[0][i] -= w*1.0/(dt*pm[0][i]*pm[0][i])*padj[i]*(-pfor[it][i]+pfor[it-1][i]); // K gradient
            g[1][i] += w*dt*(temp[2][i]*temp[0][i] + temp[3][i]*temp[1][i]); // rho-1 gradient
        }
    }
}

void nl_we_op_a::propagate(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const
{
    // wavefields allocations and pointers
    int nxz=nx*nz;
    data_t * u  = new data_t[3*nxz];
    memset(u, 0, 3*nxz*sizeof(data_t));
    data_t (*prev) [nxz] = (data_t (*)[nxz]) u;
    data_t (*curr) [nxz] = (data_t (*)[nxz]) (u + nxz);
    data_t (*next) [nxz] = (data_t (*)[nxz]) (u + 2*nxz);
    data_t (*bucket) [nxz];
    data_t (*u_full) [nxz];
    data_t * tmp;

    if (par.sub>0) u_full = (data_t (*) [nxz]) full_wfld;
    if (grad != nullptr) 
    {
        tmp = new data_t[4*nxz];
        memset(tmp, 0, 4*nxz*sizeof(data_t));
    }
	const data_t* mod[2] = {model, model+nxz};

    // free surface stiffness
    data_t alpha = 0.2508560249 / par.free_surface_stiffness;

    // #############  PML stuff ##################
    int l = std::max(par.taper_top,par.taper_bottom);
    l = std::max(l, par.taper_left);
    l = std::max(l, par.taper_right);
    // ###########################################
    
    // source and receiver components for injection/extraction
    int ns=par.ns;
    int nr=par.nr;
    if (adj)
    {
        ns = par.nr;
        nr = par.ns;
    }
    const data_t * src[1] = {allsrc};
    data_t * rcv[1] = {allrcv};

    // prev = K/2 * dt2 * src
    inj->inject(false, src, prev[0], nx, nz, par.nt, inj->_ntr, 0, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
    for (int i=0; i<nxz; i++) prev[0][i]  *= 0.5 * mod[0][i]*par.dt*par.dt;  

    int pct10 = round(par.nt/10);

    for (int it=0; it<par.nt-1; it++)
    {
        // copy the current wfld to the full wfld vector
        if ((par.sub>0) && (it%par.sub==0) && (grad==nullptr)) memcpy(u_full[it/par.sub], curr, nxz*sizeof(data_t));

        // extract receivers
        if (grad == nullptr)
        {
            ext->extract(true, curr[0], rcv, nx, nz, par.nt, ext->_ntr, it, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
        }

        // compute FWI gradients except for first and last time samples
        if ((grad != nullptr) && ((par.nt-1-it)%par.sub==0) && it!=0) compute_gradients(model, full_wfld, curr[0], tmp, grad, par, nx, nz, (par.nt-1-it)/par.sub, dx, dz, par.sub*par.dt);

        // apply spatial SBP operators
       
        Dxx_var<expr2>(false, curr[0], next[0], nx, nz, dx, 0, nx, 0, nz, mod, 1.0);
        Dzz_var<expr2>(true, curr[0], next[0], nx, nz, dz, 0, nx, 0, nz, mod, 1.0);

        // inject sources
        inj->inject(true, src, next[0], nx, nz, par.nt, inj->_ntr, it, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());

        // apply boundary conditions
        if (par.bc_top==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_top(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, alpha);
        }
        else if (par.bc_top==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_absorbing_top(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
        }
        else if (par.bc_top==3)
        {
            const data_t * in[1] = {curr[0]};
            asat_neumann_top(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
        }
        else if (par.bc_top==4)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_neumann_absorbing_top(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, _transmission->getCVals());
        }
        else if (par.bc_top==5)
        {
            const data_t * in[1] = {curr[0]};
            asat_neumann_dirichlet_top(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, alpha, _transmission->getCVals());
        }
        if (par.bc_bottom==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_bottom(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, alpha);
        }
        else if (par.bc_bottom==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_absorbing_bottom(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
        }
        else if (par.bc_bottom==3)
        {
            const data_t * in[1] = {curr[0]};
            asat_neumann_bottom(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
        }
        if (par.bc_left==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_left(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, alpha);
        }
        else if (par.bc_left==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_absorbing_left(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        else if (par.bc_left==3)
        {
            const data_t * in[1] = {curr[0]};
            asat_neumann_left(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        if (par.bc_right==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_right(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, alpha);
        }
        else if (par.bc_right==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_absorbing_right(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        else if (par.bc_right==3)
        {
            const data_t * in[1] = {curr[0]};
            asat_neumann_right(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }

        // update wfld with the 2 steps time recursion
        #pragma omp parallel for
        for (int i=0; i<nxz; i++)
        {
            next[0][i] = par.dt*par.dt*next[0][i]*mod[0][i] + 2*curr[0][i] - prev[0][i];
            //next[0][i] = par.dt*par.dt*next[0][i]*mod[0][i] + curr[0][i] - prev[0][i];
        }

        //taperz(next[0], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
        //taperz(next[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        //taperx(next[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        //taperx(next[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);

        //#pragma omp parallel for
        //for (int i=0; i<nxz; i++)
        //{
        //    next[0][i] += curr[0][i];
        //}

        // scale boundaries when relevant (for locally absorbing BC only)
        data_t * in[1] = {next[0]};
        if (par.bc_top!=4) asat_scale_boundaries(in, nx, nz, dx, dz, 0, nx, 0, nz, mod, par.dt, par.bc_top==2, par.bc_bottom==2, par.bc_left==2, par.bc_right==2);
        else asat_scale_boundaries_bis(in, nx, nz, dx, dz, 0, nx, 0, nz, mod, _transmission->getCVals(), par.dt, par.bc_top==4, par.bc_bottom==2, par.bc_left==2, par.bc_right==2);
        
        taperz(curr[0], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
        taperz(curr[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        taperx(curr[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        taperx(curr[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
        taperz(next[0], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
        taperz(next[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        taperx(next[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        taperx(next[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);

        bucket=prev;
        prev=curr;
        curr=next;
        next=bucket;

        if ((it+1) % pct10 == 0 && par.verbose>2) fprintf(stderr,"Propagation progress = %d\%\n",100*(it+1)/(par.nt-1));
        else if (it == par.nt-2 && par.verbose>2) fprintf(stderr,"Propagation progress = 100\%\n");
    }

    // copy the last wfld to the full wfld vector
    if ((par.sub>0) && ((par.nt-1)%par.sub==0) && (grad==nullptr)) memcpy(u_full[(par.nt-1)/par.sub], curr, nxz*sizeof(data_t));

    // extract receivers last sample
    if (grad == nullptr)
    {
        ext->extract(true, curr[0], rcv, nx, nz, par.nt, ext->_ntr, par.nt-1, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
    }

    // last sample gradient
    if (grad != nullptr ) compute_gradients(model, full_wfld, curr[0], tmp, grad, par, nx, nz, 0, dx, dz, par.sub*par.dt);
    
    delete [] u;
    if (grad != nullptr) delete [] tmp;
}


void nl_we_op_ae::propagate(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const
{
    // wavefields allocations and pointers
    int nxz=nx*nz;
    int nza = par.nza;
    int nxza=nx*nza; // for water layer
    data_t * u  = new data_t[6*nxz];
    data_t * dux  = new data_t[2*nxz];
    data_t * duz  = new data_t[2*nxz];
    data_t * phi  = new data_t[3*nxza]; // momentum potential in the water layer
    data_t * modela  = new data_t[2*nxza]; // acoustic model
    memset(u, 0, 6*nxz*sizeof(data_t));
    memset(dux, 0, 2*nxz*sizeof(data_t));
    memset(duz, 0, 2*nxz*sizeof(data_t));
    memset(phi, 0, 3*nxza*sizeof(data_t));
    data_t (*prev) [nxz] = (data_t (*)[nxz]) u;
    data_t (*curr) [nxz] = (data_t (*)[nxz]) (u + 2*nxz);
    data_t (*next) [nxz] = (data_t (*)[nxz]) (u + 4*nxz);
    data_t (*bucket) [nxz];
    data_t (*preva) [nxza] = (data_t (*)[nxza]) phi;
    data_t (*curra) [nxza] = (data_t (*)[nxza]) (phi + nxza);
    data_t (*nexta) [nxza] = (data_t (*)[nxza]) (phi + 2*nxza);
    data_t (*bucketa) [nxza];
    data_t (*u_x) [nxz] = (data_t (*)[nxz]) dux;
    data_t (*u_z) [nxz] = (data_t (*)[nxz]) duz;
    data_t (*u_full) [2][nxz];
    data_t (*phi_full) [nxza];
    data_t * tmp;

    data_t K = par.water_density*par.water_velocity*par.water_velocity; // bulk modulus
    data_t rho = par.water_density; // density
    for (int i=0; i<nxza; i++) {
        modela[i] = K;
        modela[nxza+i] = 1.0/rho;
    }

    if (par.sub>0) u_full = (data_t (*) [2][nxz]) full_wfld;
    if (par.sub>0 && par.acoustic_wavefield) phi_full = (data_t (*) [nxza]) _full_wflda->getVals();
    if (grad != nullptr) 
    {
        tmp = new data_t[4*nxz];
        memset(tmp, 0, 4*nxz*sizeof(data_t));
    }
	const data_t* mod[3] = {model, model+nxz, model+2*nxz};
	const data_t* moda[2] = {modela, modela+nxza};

    // free surface stiffness
    data_t alpha = 0.2508560249 / par.free_surface_stiffness;

    // #############  PML stuff ##################
    int l = std::max(par.taper_top,par.taper_bottom);
    l = std::max(l, par.taper_left);
    l = std::max(l, par.taper_right);
    // ###########################################
    
    // source and receiver components for injection/extraction
    int ns=par.ns;
    int nr=par.nr;
    if (adj)
    {
        ns = par.nr;
        nr = par.ns;
    }

    // in the elastic medium
    const data_t * srcx[2] = {nullptr, nullptr};
    const data_t * srcz[2] = {nullptr, nullptr};
    data_t * rcvx[2] = {nullptr, nullptr};
    data_t * rcvz[2] = {nullptr, nullptr};

    if (!adj){
        rcvx[0] = allrcv + nr*par.nt;
        rcvz[0] = allrcv + 2*nr*par.nt;
        if (!par.acoustic_source){
            srcx[0] = allsrc;
            srcz[0] = allsrc + ns*par.nt;
        }
    }
    else{
        srcx[0] = allsrc + ns*par.nt;
        srcz[0] = allsrc + 2*ns*par.nt;
        if (!par.acoustic_source){
            rcvx[0] = allrcv;
            rcvz[0] = allrcv + nr*par.nt;
        }
    }
    
    if (par.mt==true)
    {
        if (!adj)
        {
            srcx[1] = allsrc + 2*ns*par.nt;
            srcz[0] = allsrc + 2*ns*par.nt;
            srcz[1] = allsrc + ns*par.nt;
        }
        else
        {
            rcvx[1] = allrcv + 2*nr*par.nt;
            rcvz[0] = allrcv + 2*nr*par.nt;
            rcvz[1] = allrcv + nr*par.nt;
        }
    }

    // in the water layer
    const data_t * src[1] = {allsrc};
    data_t * rcv[1] = {allrcv};

    // injector/extractor in the water layer
    hypercube<data_t> domaina (axis<data_t>(par.nza,0,dz),axis<data_t>(nx,ox,dx));
    injector *inja=nullptr, *exta=nullptr;
    if (!adj){
        if (par.acoustic_source) inja = new delta_m3(domaina,{par.sxz[s]});
        exta = new delta_m3(domaina,par.rxza[s]);
    }
    else{
        if (par.acoustic_source) exta = new delta_m3(domaina,{par.sxz[s]});
        inja = new delta_m3(domaina,par.rxza[s]);
    }
    
    // prev = 1/(2 * rho) * dt2 * src
    if (!par.acoustic_source) {
        inj->inject(false, srcx, prev[0], nx, nz, par.nt, inj->_ntr, 0, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
        inj->inject(false, srcz, prev[1], nx, nz, par.nt, inj->_ntr, 0, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
        for (int i=0; i<nxz; i++)
        {
            prev[0][i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
            prev[1][i]  *= 0.5 * par.dt*par.dt/ mod[2][i];
        }
    }
    else { // preva = K/2 * dt2 * src
        inja->inject(false, src, preva[0], nx, nza, par.nt, inja->_ntr, 0, 0, inja->_ntr, inja->_xind.data(), inja->_zind.data(), inja->_xw.data(), inja->_zw.data());
        for (int i=0; i<nxza; i++) preva[0][i]  *= 0.5 * moda[0][i]*par.dt*par.dt;
    }   

    int pct10 = round(par.nt/10);

    for (int it=0; it<par.nt-1; it++)
    {
        // copy the current wfld to the full wfld vector
        if ((par.sub>0) && (it%par.sub==0) && (grad==nullptr)) memcpy(u_full[it/par.sub], curr, 2*nxz*sizeof(data_t));
        if ((par.sub>0) && par.acoustic_wavefield && (it%par.sub==0) && (grad==nullptr)) memcpy(phi_full[it/par.sub], curra, nxza*sizeof(data_t));

        // extract receivers
        if (grad == nullptr)
        {
            if (!adj){
                ext->extract(true, curr[0], rcvx, nx, nz, par.nt, ext->_ntr, it, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
                ext->extract(true, curr[1], rcvz, nx, nz, par.nt, ext->_ntr, it, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
                exta->extract(true, curra[0], rcv, nx, nza, par.nt, exta->_ntr, it, 0, exta->_ntr, exta->_xind.data(), exta->_zind.data(), exta->_xw.data(), exta->_zw.data());
            }
            else{
                if (!par.acoustic_source) {
                    ext->extract(true, curr[0], rcvx, nx, nz, par.nt, ext->_ntr, it, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
                    ext->extract(true, curr[1], rcvz, nx, nz, par.nt, ext->_ntr, it, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
                }
                else exta->extract(true, curra[0], rcv, nx, nza, par.nt, exta->_ntr, it, 0, exta->_ntr, exta->_xind.data(), exta->_zind.data(), exta->_xw.data(), exta->_zw.data());
            }
        }

        // compute FWI gradients except for first and last time samples
        if ((grad != nullptr) && ((par.nt-1-it)%par.sub==0) && it!=0) compute_gradients(model, full_wfld, curr[0], u_x[0], u_z[0], tmp, grad, par, nx, nz, (par.nt-1-it)/par.sub, dx, dz, par.sub*par.dt);

        // apply spatial SBP operators
        Dxx_var<expr2>(false, curra[0], nexta[0], nx, nza, dx, 0, nx, 0, nza, moda, 1.0);
        Dzz_var<expr2>(true, curra[0], nexta[0], nx, nza, dz, 0, nx, 0, nza, moda, 1.0);

        if (par.version==1)
        {
            Dxx_var<expr4a>(false, curr[0], next[0], nx, nz, dx, 0, nx, 0, nz, mod, 1.0);
            Dzz_var<expr4a>(false, curr[1], next[1], nx, nz, dz, 0, nx, 0, nz, mod, 1.0);
        }
        else
        {
            Dxx_var<expr4b>(false, curr[0], next[0], nx, nz, dx, 0, nx, 0, nz, mod, 1.0);
            Dzz_var<expr4b>(false, curr[1], next[1], nx, nz, dz, 0, nx, 0, nz, mod, 1.0);
        }
        
        Dxx_var<expr2>(true, curr[1], next[1], nx, nz, dx, 0, nx, 0, nz, mod, 1.0);
        Dzz_var<expr2>(true, curr[0], next[0], nx, nz, dz, 0, nx, 0, nz, mod, 1.0);
        
        Dx(false, curr[0], u_x[0], nx, nz, dx, 0, nx, 0, nz);
        Dz(false, curr[1], u_z[1], nx, nz, dz, 0, nx, 0, nz);
        
        if (par.version==2)
        {
            mult_Dx(true, u_x[0], next[0], nx, nz, dx, 0, nx, 0, nz, mod[0], 1.0);
            mult_Dz(true, u_z[1], next[1], nx, nz, dz, 0, nx, 0, nz, mod[0], 1.0);
        }
        
        mult_Dx(true, u_z[1], next[0], nx, nz, dx, 0, nx, 0, nz, mod[0], 1.0);
        mult_Dz(true, u_x[0], next[1], nx, nz, dz, 0, nx, 0, nz, mod[0], 1.0);

        Dx(false, curr[1], u_x[1], nx, nz, dx, 0, nx, 0, nz);
        Dz(false, curr[0], u_z[0], nx, nz, dz, 0, nx, 0, nz);
        
        mult_Dx(true, u_z[0], next[1], nx, nz, dx, 0, nx, 0, nz, mod[1], 1.0);
        mult_Dz(true, u_x[1], next[0], nx, nz, dz, 0, nx, 0, nz, mod[1], 1.0);

        // inject sources
        if (adj){
            inj->inject(true, srcx, next[0], nx, nz, par.nt, inj->_ntr, it, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
            inj->inject(true, srcz, next[1], nx, nz, par.nt, inj->_ntr, it, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
            inja->inject(true, src, nexta[0], nx, nza, par.nt, inja->_ntr, it, 0, inja->_ntr, inja->_xind.data(), inja->_zind.data(), inja->_xw.data(), inja->_zw.data());
        }
        else{
            if (!par.acoustic_source){
                inj->inject(true, srcx, next[0], nx, nz, par.nt, inj->_ntr, it, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
                inj->inject(true, srcz, next[1], nx, nz, par.nt, inj->_ntr, it, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
            }
            else{
                inja->inject(true, src, nexta[0], nx, nza, par.nt, inja->_ntr, it, 0, inja->_ntr, inja->_xind.data(), inja->_zind.data(), inja->_xw.data(), inja->_zw.data());
            }
        }

        // apply boundary conditions
        if (par.bc_top==1)
        {
            const data_t * in[1] = {curra[0]};
            asat_dirichlet_top(true, in, nexta[0], nx, nza, dx, dz, par.pml_L*l, nx-par.pml_R*l, moda, alpha);
        }
        else if (par.bc_top==2)
        {
            const data_t * in[2] = {curra[0],preva[0]};
            asat_absorbing_top(true, in, nexta[0], nx, nza, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, moda, 1.0);
        }
        else if (par.bc_top==3)
        {
            const data_t * in[1] = {curra[0]};
            asat_neumann_top(true, in, nexta[0], nx, nza, dx, dz, par.pml_L*l, nx-par.pml_R*l, moda, 1.0);
        }
        if (par.bc_bottom==1)
        {
            const data_t * in[2] = {curr[1],curr[0]};
            esat_neumann_bottom<expr2,expr2>(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            in[0] = curr[0]; in[1] = curr[1];
            if (par.version==1) esat_neumann_bottom<expr1,expr4a>(true, in, next[1], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            else
            {
                esat_neumann_bottom<expr1,expr4b>(true, in, next[1], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
                esat_Dz_bottom<expr1>(true, curr[1], next[1], nx, nz, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
        }
        else if (par.bc_bottom==2)
        {
            const data_t * in[3] = {curr[1],curr[0],prev[0]};
            esat_absorbing_bottom<expr2,expr2,expr5>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            in[0] = curr[0]; in[1] = curr[1]; in[2] = prev[1];
            if (par.version==1) esat_absorbing_bottom<expr1,expr4a,expr6>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            else
            {
                esat_absorbing_bottom<expr1,expr4b,expr6>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
                esat_Dz_bottom<expr1>(true, curr[1], next[1], nx, nz, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
        }
        if (par.bc_left==1)
        {
            const data_t * ina[1] = {curra[0]};
            asat_dirichlet_left(true, ina, nexta[0], nx, nza, dx, dz, par.pml_T*l, nza-par.pml_B*l, moda, alpha);

            const data_t * in[2] = {curr[1],curr[0]};
            if (par.version==1) esat_neumann_left<expr1,expr4a>(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            else
            {
                esat_neumann_left<expr1,expr4b>(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
                esat_Dx_left<expr1>(true, curr[0], next[0], nx, nz, dx, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            }
            in[0] = curr[0]; in[1] = curr[1];
            esat_neumann_left<expr2,expr2>(true, in, next[1], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        else if (par.bc_left==2)
        {
            const data_t * ina[2] = {curra[0],preva[0]};
            asat_absorbing_left(true, ina, nexta[0], nx, nza, dx, dz, par.dt, par.pml_T*l, nza-par.pml_B*l, moda, 1.0);

            const data_t * in[3] = {curr[1],curr[0],prev[0]};
            if (par.version==1) esat_absorbing_left<expr1,expr4a,expr6>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            else
            {
                esat_absorbing_left<expr1,expr4b,expr6>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
                esat_Dx_left<expr1>(true, curr[0], next[0], nx, nz, dx, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            }
            in[0] = curr[0]; in[1] = curr[1]; in[2] = prev[1];
            esat_absorbing_left<expr2,expr2,expr5>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        if (par.bc_right==1)
        {
            const data_t * ina[1] = {curra[0]};
            asat_dirichlet_right(true, ina, nexta[0], nx, nza, dx, dz, par.pml_T*l, nza-par.pml_B*l, moda, alpha);

            const data_t * in[2] = {curr[1],curr[0]};
            if (par.version==1) esat_neumann_right<expr1,expr4a>(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            else
            {
                esat_neumann_right<expr1,expr4b>(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
                esat_Dx_right<expr1>(true, curr[0], next[0], nx, nz, dx, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            }
            in[0] = curr[0]; in[1] = curr[1];
            esat_neumann_right<expr2,expr2>(true, in, next[1], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        else if (par.bc_right==2)
        {
            const data_t * ina[2] = {curra[0],preva[0]};
            asat_absorbing_right(true, ina, nexta[0], nx, nza, dx, dz, par.dt, par.pml_T*l, nza-par.pml_B*l, moda, 1.0);

            const data_t * in[3] = {curr[1],curr[0],prev[0]};
            if (par.version==1) esat_absorbing_right<expr1,expr4a,expr6>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            else
            {
                esat_absorbing_right<expr1,expr4b,expr6>(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
                esat_Dx_right<expr1>(true, curr[0], next[0], nx, nz, dx, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            }
            in[0] = curr[0]; in[1] = curr[1]; in[2] = prev[1];
            esat_absorbing_right<expr2,expr2,expr5>(true, in, next[1], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }

        // apply boundary conditions at the interface without coupling 
        {
            const data_t * ina[1] = {curra[0]};
            asat_neumann_bottom(true, ina, nexta[0], nx, nza, dx, dz, par.pml_L*l, nx-par.pml_R*l, moda, 1.0);

            const data_t * in[2] = {curr[1],curr[0]};
            esat_neumann_top<expr2,expr2>(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            in[0] = curr[0]; in[1] = curr[1];
            if (par.version==1) esat_neumann_top<expr1,expr4a>(true, in, next[1], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            else
            {
                esat_neumann_top<expr1,expr4b>(true, in, next[1], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
                esat_Dz_top<expr1>(true, curr[1], next[1], nx, nz, dz, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
        }

        // update wfld with the 2 steps time recursion
        #pragma omp parallel for
        for (int i=0; i<nxza; i++)
        {
            nexta[0][i] = par.dt*par.dt*nexta[0][i]*moda[0][i] + 2*curra[0][i] - preva[0][i];
        }
        #pragma omp parallel for
        for (int i=0; i<nxz; i++)
        {
            next[0][i] = par.dt*par.dt*next[0][i]/mod[2][i] + 2*curr[0][i] - prev[0][i];
            next[1][i] = par.dt*par.dt*next[1][i]/mod[2][i] + 2*curr[1][i] - prev[1][i];
        }

        // scale boundaries when relevant (for locally absorbing BC only)
        data_t * ina[1] = {nexta[0]};
        data_t * in[2] = {next[0],next[1]};
        asat_scale_boundaries(ina, nx, nza, dx, dz, 0, nx, 0, nza, moda, par.dt, par.bc_top==2, false, par.bc_left==2, par.bc_right==2);
        esat_scale_boundaries(in, nx, nz, dx, dz, 0, nx, 0, nz, mod, par.dt, false, par.bc_bottom==2, par.bc_left==2, par.bc_right==2);

        // apply the coupling at the interface
        // next_phi = next_phi +- K.dt/2.1/Hz.(next_uz - prev_uz)
        // next_uz = next_uz -+ dt/(2.rho).1/Hz.(next_phi - prev_phi)
        int sign=(1-2*adj); // = 1 for forward modeling and = -1 for adjoint modeling
        data_t h = dz*(17.0/48);
        data_t val0 = par.dt/(2*h);
        if (it==-1){
            #pragma omp parallel for
            for (int ix=0; ix<nx; ix++)
            {
                data_t val1 =  K*2*val0*(next[1][ix*nz]-curr[1][ix*nz] + sign*2*val0/mod[2][ix*nz]*curra[0][ix*nza+nza-1]); // K.dt.1/Hz.(next_uz - curr_uz +- dt/rho.1/Hz.curr_phi)
                data_t val2 = 1 + 4*val0*val0*K/mod[2][ix*nz]; // 1 + K.dt2/rho.1/Hz.1/Hz
                nexta[0][ix*nza+nza-1] = (nexta[0][ix*nza+nza-1] + sign*val1)/val2;
                next[1][ix*nz] = next[1][ix*nz] - sign*2*val0/mod[2][ix*nz]*(nexta[0][ix*nza+nza-1] - curra[0][ix*nza+nza-1]);
            }
        }
        else{
            #pragma omp parallel for
            for (int ix=0; ix<nx; ix++)
            {
                data_t val1 =  K*val0*(next[1][ix*nz]-prev[1][ix*nz] + sign*val0/mod[2][ix*nz]*preva[0][ix*nza+nza-1]); // K.dt/2.1/Hz.(next_uz - prev_uz +- dt/(2.rho).1/Hz.prev_phi)
                data_t val2 = 1 + val0*val0*K/mod[2][ix*nz]; // 1 + K.dt2/(4.rho).1/Hz.1/Hz
                nexta[0][ix*nza+nza-1] = (nexta[0][ix*nza+nza-1] + sign*val1)/val2;
                next[1][ix*nz] = next[1][ix*nz] - sign*val0/mod[2][ix*nz]*(nexta[0][ix*nza+nza-1] - preva[0][ix*nza+nza-1]);
            }
        }

        taperz(curra[0], nx, nza, 0, nx, par.taper_top, 0, par.taper_strength);
        taperx(curra[0], nx, nza, 0, nza, par.taper_left, 0, par.taper_strength);
        taperx(curra[0], nx, nza, 0, nza, nx-par.taper_right, nx, par.taper_strength);
        taperz(nexta[0], nx, nza, 0, nx, par.taper_top, 0, par.taper_strength);
        taperx(nexta[0], nx, nza, 0, nza, par.taper_left, 0, par.taper_strength);
        taperx(nexta[0], nx, nza, 0, nza, nx-par.taper_right, nx, par.taper_strength);

        taperz(curr[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        taperx(curr[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        taperx(curr[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
        taperz(curr[1], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        taperx(curr[1], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        taperx(curr[1], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
        taperz(next[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        taperx(next[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        taperx(next[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
        taperz(next[1], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        taperx(next[1], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        taperx(next[1], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
        
        bucketa=preva;
        preva=curra;
        curra=nexta;
        nexta=bucketa;

        bucket=prev;
        prev=curr;
        curr=next;
        next=bucket;

        if ((it+1) % pct10 == 0 && par.verbose>2) fprintf(stderr,"Propagation progress = %d\%\n",100*(it+1)/(par.nt-1));
        else if (it == par.nt-2 && par.verbose>2) fprintf(stderr,"Propagation progress = 100\%\n");
    }

    // copy the last wfld to the full wfld vector
    if ((par.sub>0) && ((par.nt-1)%par.sub==0) && (grad==nullptr)) memcpy(u_full[(par.nt-1)/par.sub], curr, 2*nxz*sizeof(data_t));
    if ((par.sub>0) && par.acoustic_wavefield && ((par.nt-1)%par.sub==0) && (grad==nullptr)) memcpy(phi_full[(par.nt-1)/par.sub], curra, nxza*sizeof(data_t));

    // extract receivers last sample
    if (grad == nullptr)
    {
        if (!adj){
            ext->extract(true, curr[0], rcvx, nx, nz, par.nt, ext->_ntr, par.nt-1, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
            ext->extract(true, curr[1], rcvz, nx, nz, par.nt, ext->_ntr, par.nt-1, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
            exta->extract(true, curra[0], rcv, nx, nza, par.nt, exta->_ntr, par.nt-1, 0, exta->_ntr, exta->_xind.data(), exta->_zind.data(), exta->_xw.data(), exta->_zw.data());
        }
        else{
            if (!par.acoustic_source) {
                ext->extract(true, curr[0], rcvx, nx, nz, par.nt, ext->_ntr, par.nt-1, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
                ext->extract(true, curr[1], rcvz, nx, nz, par.nt, ext->_ntr, par.nt-1, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
            }
            else exta->extract(true, curra[0], rcv, nx, nza, par.nt, exta->_ntr, par.nt-1, 0, exta->_ntr, exta->_xind.data(), exta->_zind.data(), exta->_xw.data(), exta->_zw.data());
        }
    }
    
    delete [] u;
    delete [] dux;
    delete [] duz;
    delete [] phi;
    delete [] modela;
    delete inja;
    delete exta;
    if (grad != nullptr) delete [] tmp;
}



void born_op_a::compute_gradients(const data_t * model, const data_t * u_full, const data_t * curr, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const
{
    // Grad_s2 = adjoint.H.Ht.forward_tt
    // The H quadrature will be applied elsewhere to the final gradients (all shots included)

    int nxz = nx*nz;
    data_t (*pfor) [nxz] = (data_t (*) [nxz]) u_full;
    const data_t *padj = curr;
    data_t (*g) [nxz] = (data_t (*)[nxz]) grad;

    // different from first and last time sample
    if (it>0 && it<par.nt/par.sub){
        #pragma omp parallel for
        for (int i=0; i<nxz; i++) 
        {
            g[0][i] += 1.0/dt*padj[i]*(pfor[it+1][i] - 2*pfor[it][i] + pfor[it-1][i]);
            // g[0][i] += dt*padj[i]*pfor[it][i];
        }
    }
    // time zero
    else if(it==0){
        #pragma omp parallel for
        for (int i=0; i<nxz; i++) 
        {
            g[0][i] += 1.0/dt*padj[i]*(pfor[it+1][i] - pfor[it][i]);
            // g[0][i] += 0.5*dt*padj[i]*pfor[it][i];
        }
    }
    // tmax
    else{
        #pragma omp parallel for
        for (int i=0; i<nxz; i++) 
        {
            g[0][i] += 1.0/dt*padj[i]*(- pfor[it][i] + pfor[it-1][i] );
            // g[0][i] += 0.5*dt*padj[i]*pfor[it][i];
        }
    }
}

void born_op_a::propagate(bool adj, const data_t * model, const data_t * allsrc, data_t * allrcv, const injector * inj, const injector * ext, data_t * full_wfld, data_t * grad, const param &par, int nx, int nz, data_t dx, data_t dz, data_t ox, data_t oz, int s) const
{
    // wavefields allocations and pointers
    int nxz=nx*nz;
    data_t * u  = new data_t[3*nxz]; // background wavefield
    memset(u, 0, 3*nxz*sizeof(data_t));
    data_t * du  = nullptr; // scattered wavefield
    if (!adj) {
        du  = new data_t[3*nxz];
        memset(du, 0, 3*nxz*sizeof(data_t));
    }
    data_t (*prev) [nxz] = (data_t (*)[nxz]) u;
    data_t (*curr) [nxz] = (data_t (*)[nxz]) (u + nxz);
    data_t (*next) [nxz] = (data_t (*)[nxz]) (u + 2*nxz);
    data_t (*dprev) [nxz] = (data_t (*)[nxz]) du;
    data_t (*dcurr) [nxz] = (data_t (*)[nxz]) (du + nxz);
    data_t (*dnext) [nxz] = (data_t (*)[nxz]) (du + 2*nxz);
    data_t (*bucket) [nxz];
    data_t (*u_full) [nxz];

    if (par.sub>0) u_full = (data_t (*) [nxz]) full_wfld;

    // build K,1/rho model with rho=1
    data_t * tmp=new data_t[nxz];
    data_t * tmp2=new data_t[nxz];
    for (int i=0; i<nxz; i++) {
        tmp[i] = model[i]*model[i];
        tmp2[i] = 1.0;
    }
    const data_t* mod[2] = {tmp, tmp2};


    // free surface stiffness
    data_t alpha = 0.2508560249 / par.free_surface_stiffness;

    // #############  PML stuff ##################
    int l = std::max(par.taper_top,par.taper_bottom);
    l = std::max(l, par.taper_left);
    l = std::max(l, par.taper_right);
    // ###########################################
    
    // source and receiver components for injection/extraction
    int ns=par.ns;
    int nr=par.nr;
    if (adj)
    {
        ns = par.nr;
        nr = par.ns;
    }
    const data_t * src[1] = {allsrc};
    data_t * rcv[1] = {allrcv};

    // prev = K/2 * dt2 * src
    inj->inject(false, src, prev[0], nx, nz, par.nt, inj->_ntr, 0, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());
    for (int i=0; i<nxz; i++) prev[0][i] *= 0.5 * mod[0][i]*par.dt*par.dt;
    if (!adj) for (int i=0; i<nxz; i++) dprev[0][i] = -mod[0][i]*prev[0][i]*grad[i];


    int pct10 = round(par.nt/10);

    for (int it=0; it<par.nt-1; it++)
    {
        // copy the current wfld to the full wfld vector
        if ((par.sub>0) && (it%par.sub==0) && (!adj)) memcpy(u_full[it/par.sub], curr, nxz*sizeof(data_t));

        // compute FWI gradients except for first and last time samples
        if ((adj) && ((par.nt-1-it)%par.sub==0) && it!=0) compute_gradients(model, full_wfld, curr[0], tmp, grad, par, nx, nz, (par.nt-1-it)/par.sub, dx, dz, par.sub*par.dt);

        // apply spatial SBP operators
        Dxx(false, curr[0], next[0], nx, nz, dx, 0, nx, 0, nz);
        Dzz(true, curr[0], next[0], nx, nz, dz, 0, nx, 0, nz);

        // inject sources
        inj->inject(true, src, next[0], nx, nz, par.nt, inj->_ntr, it, 0, inj->_ntr, inj->_xind.data(), inj->_zind.data(), inj->_xw.data(), inj->_zw.data());

        // apply boundary conditions
        if (par.bc_top==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_top(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, alpha);
        }
        else if (par.bc_top==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_absorbing_top(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
        }
        if (par.bc_bottom==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_bottom(true, in, next[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, alpha);
        }
        else if (par.bc_bottom==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_absorbing_bottom(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
        }
        if (par.bc_left==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_left(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, alpha);
        }
        else if (par.bc_left==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_absorbing_left(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }
        if (par.bc_right==1)
        {
            const data_t * in[1] = {curr[0]};
            asat_dirichlet_right(true, in, next[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, alpha);
        }
        else if (par.bc_right==2)
        {
            const data_t * in[2] = {curr[0],prev[0]};
            asat_absorbing_right(true, in, next[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
        }

        // update wfld with the 2 steps time recursion
        #pragma omp parallel for
        for (int i=0; i<nxz; i++) next[0][i] = par.dt*par.dt*mod[0][i]*next[0][i] + 2*curr[0][i] - prev[0][i];

        // scale boundaries when relevant (for locally absorbing BC only)
        data_t * in[1] = {next[0]};
        asat_scale_boundaries(in, nx, nz, dx, dz, 0, nx, 0, nz, mod, par.dt, par.bc_top==2, par.bc_bottom==2, par.bc_left==2, par.bc_right==2);

        taperz(curr[0], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
        taperz(curr[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        taperx(curr[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        taperx(curr[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
        taperz(next[0], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
        taperz(next[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
        taperx(next[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
        taperx(next[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);

        // copy the second time derivative of the current wfld to the full wfld vector
        // if ((par.sub>0) && (it%par.sub==0) && (!adj)) {
        //     #pragma omp parallel for
        //     for (int i=0; i<nxz; i++) u_full[it/par.sub][i] = (next[0][i] - 2*curr[0][i] + prev[0][i])/(par.dt*par.dt);
        // }

        // propagate scattered wavefield
        if (!adj)
        {
            // extract receivers
            ext->extract(true, dcurr[0], rcv, nx, nz, par.nt, ext->_ntr, it, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());

            // apply spatial SBP operators
            Dxx(false, dcurr[0], dnext[0], nx, nz, dx, 0, nx, 0, nz);
            Dzz(true, dcurr[0], dnext[0], nx, nz, dz, 0, nx, 0, nz);

            // add source for the scattered wfld
            #pragma omp parallel for
            for (int i=0; i<nxz; i++) dnext[0][i] -= (next[0][i] - 2*curr[0][i] + prev[0][i])/(par.dt*par.dt)*grad[i];

            // apply boundary conditions to the scattered wfld
            if (par.bc_top==1)
            {
                    const data_t * din[1] = {dcurr[0]};
                    asat_dirichlet_top(true, din, dnext[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, alpha);
            }
            else if (par.bc_top==2)
            {
                    const data_t * din[2] = {dcurr[0],dprev[0]};
                    asat_absorbing_top(true, din, dnext[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
            if (par.bc_bottom==1)
            {
                    const data_t * din[1] = {dcurr[0]};
                    asat_dirichlet_bottom(true, din, dnext[0], nx, nz, dx, dz, par.pml_L*l, nx-par.pml_R*l, mod, alpha);
            }
            else if (par.bc_bottom==2)
            {
                    const data_t * din[2] = {dcurr[0],dprev[0]};
                    asat_absorbing_bottom(true, din, dnext[0], nx, nz, dx, dz, par.dt, par.pml_L*l, nx-par.pml_R*l, mod, 1.0);
            }
            if (par.bc_left==1)
            {
                    const data_t * din[1] = {dcurr[0]};
                    asat_dirichlet_left(true, din, dnext[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, alpha);
            }
            else if (par.bc_left==2)
            {
                    const data_t * din[2] = {dcurr[0],dprev[0]};
                    asat_absorbing_left(true, din, dnext[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            }
            if (par.bc_right==1)
            {
                    const data_t * din[1] = {dcurr[0]};
                    asat_dirichlet_right(true, din, dnext[0], nx, nz, dx, dz, par.pml_T*l, nz-par.pml_B*l, mod, alpha);
            }
            else if (par.bc_right==2)
            {
                    const data_t * din[2] = {dcurr[0],dprev[0]};
                    asat_absorbing_right(true, din, dnext[0], nx, nz, dx, dz, par.dt, par.pml_T*l, nz-par.pml_B*l, mod, 1.0);
            }

            // update scattered wfld with the 2 steps time recursion
            #pragma omp parallel for
            for (int i=0; i<nxz; i++) dnext[0][i] = par.dt*par.dt*dnext[0][i]*mod[0][i] + 2*dcurr[0][i] - dprev[0][i];

            // scale boundaries when relevant (for locally absorbing BC only)
            data_t * din[1] = {dnext[0]};
            asat_scale_boundaries(din, nx, nz, dx, dz, 0, nx, 0, nz, mod, par.dt, par.bc_top==2, par.bc_bottom==2, par.bc_left==2, par.bc_right==2);

            taperz(dcurr[0], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
            taperz(dcurr[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
            taperx(dcurr[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
            taperx(dcurr[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);
            taperz(dnext[0], nx, nz, 0, nx, par.taper_top, 0, par.taper_strength);
            taperz(dnext[0], nx, nz, 0, nx, nz-par.taper_bottom, nz, par.taper_strength);
            taperx(dnext[0], nx, nz, 0, nz, par.taper_left, 0, par.taper_strength);
            taperx(dnext[0], nx, nz, 0, nz, nx-par.taper_right, nx, par.taper_strength);

            bucket=dprev;
            dprev=dcurr;
            dcurr=dnext;
            dnext=bucket;
        }

        bucket=prev;
        prev=curr;
        curr=next;
        next=bucket;

        if ((it+1) % pct10 == 0 && par.verbose>2) fprintf(stderr,"Propagation progress = %d\%\n",100*(it+1)/(par.nt-1));
        else if (it == par.nt-2 && par.verbose>2) fprintf(stderr,"Propagation progress = 100\%\n");
    }

    // copy the last wfld (or second time derivative) to the full wfld vector (not needed)
    if ((par.sub>0) && ((par.nt-1)%par.sub==0) && (!adj)) memcpy(u_full[(par.nt-1)/par.sub], curr, nxz*sizeof(data_t));
    // if ((par.sub>0) && ((par.nt-1)%par.sub==0) && (!adj)) {
    //     #pragma omp parallel for
    //     //for (int i=0; i<nxz; i++) u_full[(par.nt-1)/par.sub][i] = ( - 2*curr[0][i] + prev[0][i])/(par.dt*par.dt);
    //     for (int i=0; i<nxz; i++) u_full[(par.nt-1)/par.sub][i] = 0;
    // }

    // extract receivers last sample
    if (!adj)
    {
        ext->extract(true, dcurr[0], rcv, nx, nz, par.nt, ext->_ntr, par.nt-1, 0, ext->_ntr, ext->_xind.data(), ext->_zind.data(), ext->_xw.data(), ext->_zw.data());
    }

    // last sample gradient
    if (adj) compute_gradients(model, full_wfld, curr[0], tmp, grad, par, nx, nz, 0, dx, dz, par.sub*par.dt);
    
    delete [] u;
    if (!adj) delete [] du;
    delete [] tmp;
    delete [] tmp2;
}

void born_op_a::apply_forward(bool add, const data_t * pmod, data_t * pdat)
{
    const data_t * pm = _model->getCVals();

    int nx = _domain.getAxis(2).n;
    int nz = _domain.getAxis(1).n;
    data_t dx = _domain.getAxis(2).d;
    data_t dz = _domain.getAxis(1).d;
    data_t ox = _domain.getAxis(2).o;
    data_t oz = _domain.getAxis(1).o;
    hypercube<data_t> domain = *_allsrc->getHyper();

    if (_par.sub>0){
        _background_wfld=std::make_shared<vecReg<data_t> > (hypercube<data_t>(_domain.getAxis(1),_domain.getAxis(2), axis<data_t>(1+_par.nt/_par.sub,0,_par.dt*_par.sub), axis<data_t>(_par.ns,0,1)));
    }

    // setting up and apply the time resampling operator
    resampler * resamp;

    axis<data_t> T = domain.getAxis(1);
    axis<data_t> Xs = domain.getAxis(2);
    axis<data_t> Xr0 = _range.getAxis(2);
    axis<data_t> Cr0 = _range.getAxis(3);
    axis<data_t> Xr=Xr0;
    int ncomp=1;
    Xs.n = domain.getN123()/T.n;
    data_t alpha = _par.dt/T.d; // ratio between sampling rates
    T.n = _par.nt;
    T.d = _par.dt;
    hypercube<data_t> hyper_s(T,Xs);
    hypercube<data_t> hyper_r(T,Xr);
    hypercube<data_t> hyper_r0(T,Xr0,Cr0);
    std::shared_ptr<vecReg<data_t> > allsrc = std::make_shared<vecReg<data_t> >(hyper_s);
    std::shared_ptr<vecReg<data_t> > allrcv = std::make_shared<vecReg<data_t> >(hyper_r);
    allsrc->zero();
    allrcv->zero();

    if (_par.resampling == "linear") resamp = new linear_resampler(domain, hyper_s);
    else resamp = new sinc_resampler(domain, hyper_s, _par.sinc_half_length);

    // resample in time (interpolation)
    resamp->apply_forward(false,_allsrc->getCVals(),allsrc->getVals());

    int size=1, rank=0, rank0=0;
#ifdef ENABLE_MPI
    if (!_par.skip_mpi){
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    }
    MPI_Comm_rank(MPI_COMM_WORLD,&rank0);
#endif

    // loop over shots
    for (int s=rank; s<_par.ns; s += size)
    {
        if (_par.verbose>1) fprintf(stderr,"Start propagating shot %d by process %d\n",s,rank0);

        // cumulative number of receivers
        int nr=0;
        if (s>0) {for (int i=0; i<s; i++) nr += _par.rxz[i].size();}

        // setting up the injection and extraction operators
        injector * inj;
        if (_par.mt==true) inj = new ddelta_m3(_domain,{_par.sxz[s]});
        else inj = new delta_m3(_domain,{_par.sxz[s]});

        injector * ext;
        if (_par.gl<=0) ext = new delta_m3(_domain,_par.rxz[s]);
        else ext = new dipole_m3(_domain,_par.rxz[s],_par.gl);

        // perform the wave propagation
        data_t * full = nullptr;
        if (_par.sub>0) {
            full = _background_wfld->getVals() + s*nx*nz*(1+_par.nt/_par.sub);
        }

        data_t * ds2 = const_cast<data_t *> (pmod);
        propagate(false, pm, allsrc->getCVals()+s*_par.nt, allrcv->getVals()+nr*_par.nt, inj, ext, full, ds2, _par, nx, nz, dx, dz, ox, oz, s);

        if (_par.verbose>1) fprintf(stderr,"Finish propagating shot %d by process %d\n",s, rank0);

        delete inj;
        delete ext;
    }

#ifdef ENABLE_MPI
    if (!_par.skip_mpi && size>1){
        data_t * temp;
        if (rank==0) temp = new data_t[allrcv->getN123()];
        if (sizeof(data_t)==8) MPI_Reduce(allrcv->getCVals(), temp, allrcv->getN123(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        else MPI_Reduce(allrcv->getCVals(), temp, allrcv->getN123(), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank==0) {
            memcpy(allrcv->getVals(), temp, allrcv->getN123()*sizeof(data_t));
            delete [] temp;
        }

        if (_par.sub>0){
            data_t * temp;
            if (rank==0) temp = new data_t[_background_wfld->getN123()];
            if (sizeof(data_t)==8) MPI_Reduce(_background_wfld->getCVals(), temp, _background_wfld->getN123(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            else MPI_Reduce(_background_wfld->getCVals(), temp, _background_wfld->getN123(), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            if (rank==0) {
                memcpy(_background_wfld->getVals(), temp, _background_wfld->getN123()*sizeof(data_t));
                delete [] temp;
            }
        }
    }
#endif              

    // resample in time (decimation)
    //fprintf(stderr,"Resample in time all receiver traces\n");
    allrcv->scale(alpha);
    resamp->setDomainRange(_range, hyper_r0);
    resamp->apply_adjoint(add, pdat, allrcv->getCVals());
        
    delete resamp;
}

void born_op_a::apply_adjoint(bool add, data_t * pmod, const data_t * pdat)
{
    const data_t * pm0 = _model->getCVals();

    int nx = _domain.getAxis(2).n;
    int nz = _domain.getAxis(1).n;
    data_t dx = _domain.getAxis(2).d;
    data_t dz = _domain.getAxis(1).d;
    data_t ox = _domain.getAxis(2).o;
    data_t oz = _domain.getAxis(1).o;
    hypercube<data_t> domain = *_allsrc->getHyper();

    data_t * temp;
    int nm=1;
    if (!add) memset(pmod, 0, nm*nx*nz*sizeof(data_t));
    else { // copy pre-existant gradient to a temporary container
        temp = new data_t[nm*nx*nz];
        memcpy(temp, pmod, nm*nx*nz*sizeof(data_t));
    }
 
    // setting up and apply the time resampling operator
    resampler * resamp;

    axis<data_t> T = domain.getAxis(1);
    axis<data_t> Xs = domain.getAxis(2);
    axis<data_t> Xr0 = _range.getAxis(2);
    axis<data_t> Cr0 = _range.getAxis(3);
    axis<data_t> Xr=Xr0;
    Xs.n = domain.getN123()/T.n;
    data_t alpha = _par.dt/T.d; // ratio between sampling rates
    T.n = _par.nt;
    T.d = _par.dt;
    hypercube<data_t> hyper_s(T,Xs);
    hypercube<data_t> hyper_r(T,Xr);
    hypercube<data_t> hyper_r0(T,Xr0,Cr0);
    std::shared_ptr<vecReg<data_t> > allsrc = std::make_shared<vecReg<data_t> >(hyper_s);
    std::shared_ptr<vecReg<data_t> > allrcv = std::make_shared<vecReg<data_t> >(hyper_r);
    allsrc->zero();
    allrcv->zero();

    if (_par.resampling == "linear") resamp = new linear_resampler(_range, hyper_r0);
    else resamp = new sinc_resampler(_range, hyper_r0, _par.sinc_half_length);

    // resample in time (interpolation)
    resamp->apply_forward(false, pdat, allrcv->getVals());

    // apply inverse time quadrature, revert in time and multiply by -1
    applyHt(true, false, allrcv->getCVals(), allrcv->getVals(), Xr.n, T.n, T.d, 0, Xr.n);
    allrcv->revert(1);
    allrcv->scale(-alpha);

    int size=1, rank=0, rank0=0;
#ifdef ENABLE_MPI
    if (!_par.skip_mpi){
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    }
    MPI_Comm_rank(MPI_COMM_WORLD,&rank0);
#endif

    // loop over shots
    for (int s=rank; s<_par.ns; s += size)
    {
        if (_par.verbose>1) fprintf(stderr,"Start back-propagating shot %d by process %d\n",s,rank0);

        // cumulative number of receivers
        int nr=0;
        if (s>0) {for (int i=0; i<s; i++) nr += _par.rxz[i].size();}

        // setting up the injection and extraction operators
        injector * ext;
        if (_par.mt==true) ext = new ddelta_m3(_domain,{_par.sxz[s]});
        else ext = new delta_m3(_domain,{_par.sxz[s]});

        injector * inj;
        if (_par.gl<=0) inj = new delta_m3(_domain,_par.rxz[s]);
        else inj = new dipole_m3(_domain,_par.rxz[s],_par.gl);

        // perform the wave propagation
        data_t * full = nullptr;
        if (_par.sub>0) {
            full = _background_wfld->getVals() + s*nx*nz*(1+_par.nt/_par.sub);
        }

        propagate(true, pm0, allrcv->getCVals()+nr*_par.nt, allsrc->getVals()+s*_par.nt, inj, ext, full, pmod, _par, nx, nz, dx, dz, ox, oz, s);
        
        if (_par.verbose>1) fprintf(stderr,"Finish back-propagating shot %d with gradient computation by process %d\n",s, rank0);
        
        delete inj;
        delete ext;
    }

#ifdef ENABLE_MPI
    if (!_par.skip_mpi && size>1){
        data_t * gtemp = new data_t[nm*nx*nz];
        if (sizeof(data_t)==8) MPI_Allreduce(pmod, gtemp, nm*nx*nz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        else MPI_Allreduce(pmod, gtemp, nm*nx*nz, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        memcpy(pmod, gtemp, nm*nx*nz*sizeof(data_t));
        delete [] gtemp;
    }
#endif

    // apply H quadrature to final gradient
    int nxz=nx*nz;
    applyHz(false, false, pmod, pmod, nx, nz, dz, 0, nx, 0, nz);
    applyHx(false, false, pmod, pmod, nx, nz, dx, 0, nx, 0, nz);

    if(add)
    {
        #pragma omp parallel for
        for(int i=0; i<nm*nxz; i++) pmod[i] += temp[i];

        delete [] temp;
    }

    delete resamp;
}