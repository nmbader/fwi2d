#pragma once

#include <fstream>
#include <sstream>
#include "vecReg.hpp"

// structure containing all the necessary parameters for wave propagation and inversion
struct param{
    // time parameters
    int nt=0;
    int sub=0;
    data_t courant=0.6, dt=-1, tmax=0, fmax=10;
    std::string resampling = "sinc";
    int sinc_half_length=11;

    // sources and receivers geometry
    int ns=1, nr=1;
    data_t sx0=0, sz0=0, rx0=0, rz0=0, rdip0=0, sxinc=0, szinc=0, rxinc=0, rzinc=0, rdipinc=0;
    std::string srcoord = "none";
    bool srcoord_from_file = false;
    std::vector<std::vector<data_t> > sxz;
    std::vector<std::vector<std::vector<data_t> > > rxz;
    std::vector<std::vector<std::vector<data_t> > > rxza;
    std::vector<data_t> rdip;
    data_t gl=0;

    // source mechanism and receiver type
    bool mt=false;
    data_t fangle=0;
    data_t mxx=1, mzz=1, mxz=1;
    int seismotype=0;

    // boundary parameters
    int bc_top=1, bc_bottom=1, bc_left=1, bc_right=1, taper_top=0, taper_bottom=0, taper_left=0, taper_right=0;
    data_t free_surface_stiffness=1.05, taper_strength=0.05;
    bool pml=false, pml_T=false, pml_B=false, pml_L=false, pml_R=false;
    data_t R = 1e-3; int p=2; 

    // model bounds
    data_t vpmin=0.2, vpmax=8, vsmin=0.1, vsmax=5, rhomin=0, rhomax=8, deltamin=-0.5, deltamax=1, epsilonmin=0, epsilonmax=1;
    data_t vmax=8, vmin=0.2;

    // for acoustic-elastic coupling
    bool acoustic_elastic=false, acoustic_source=true, acoustic_wavefield=false;
    int nza=0;
    data_t water_velocity=1.5, water_density=1;

    // model parameterization
    int nmodels=3;
    bool bsplines=false, soft_clip=false, inversion1d=false;
    int bs_nx=3, bs_nz=3;
    std::vector<int> bs_mx, bs_mz;
    std::vector<data_t> bs_controlx, bs_controlz;
    std::string horizon_file = "none";

    // inversion parameters
    std::string nlsolver="lbfgs", lsearch="regular_wolfe", mask_file="none", weights_file="none", inverse_diagonal_hessian_file="none";
    std::string prior_file = "none";
    data_t threshold=0, lambda=0, reg_xweight=1, reg_zweight=1, scale_source_log_clip=1; // threshold to stop the solver, lambda to control total regularization weight, reg_x-zweight to control regularization directional weights, scale_... to clip source scalers 
    int lbfgs_m=5, niter=0, max_trial=10, isave=10, envelop=0, regularization=-1, scale_source_times=0; // envelop in {0,1,2}, regul in {-1,0,1,2}, scale_... specifies for how many trials the source scalers will be computed
    bool normalize=0, integrate=0, double_difference=0, interferometry=0; // trace-by-trace normalization, time integration, double difference between consecutive traces, interferometry by deconvolution of consecutive traces, over residuals
    data_t wmin=0, wmax=1, epsilon=0.1; int smth_half_length=5; // parameters needed for interferometric objective function
    
    // Dev parameters with limited implementation
    data_t gdamp_tmin=0, gdamp_tmax=0, gdamp_power=1; // gradient damping during time integration ; cosine square damping between time gdamp_tmin and gdamp_tmax with power gdamp_power;

    // advanced line search parameters (check nlsolver.hpp)
    data_t ls_a0=1, ls_a1=0, ls_c1=1e-4, ls_c2=0.9, ls_max_step=1e6;
    bool ls_version=0;

    // miscallenous
    int version=2; // spatial operators version (1 or 2) for elastic WE
    int device=0; // first gpu device to use
    int nthreads=0; // maximum number of threads to use in openMP
    int verbose=1; // logging level {0,1,2,3}
    bool skip_mpi=false; // parameter to be used when running inversion with MPI
    bool format=0; // data format for read/write : SEP (0) or binary (1)
    std::string datapath="none"; // datapath for binary outputs
};

template <typename T> T convert_to (const std::string &str)
{
    std::istringstream ss(str);
    T num;
    ss >> num;
    return num;
}

// read parameter from command line arguments and store its value into a variable
// the argument must be of the form "param=value"
template<typename T>
void readParam(int argc, char **argv, std::string par, T &value){
    
    bool answer = false;
	int i = 1;
	par = par + "=";
	std::string str;

	while ( (answer == false) && (i < argc) ) {
		str = argv[i];
		if (str.substr(0, par.length()) == par){
			answer = true;
			str.erase(0, par.length());
            value = convert_to<T>(str);
		}
		i++;
	}
}

// read parameter list from command line arguments and store their values into a vector
// the argument must be of the form "param=val1,val2,val3..."
template<typename T>
void readParam(int argc, char **argv, std::string par, std::vector<T> &value){
    
    bool answer = false;
	int i = 1;
	par = par + "=";
	std::string str;

	while ( (answer == false) && (i < argc) ) {
		str = argv[i];
		if (str.substr(0, par.length()) == par){
			answer = true;
            value = {};
            str.erase(0, par.length());
			std::size_t current, previous = 0;
            current = str.find(',');
            while (current != std::string::npos) {
                value.push_back(convert_to<T>(str.substr(previous, current - previous)));
                previous = current + 1;
                current = str.find(',', previous);
            }
            value.push_back(convert_to<T>(str.substr(previous, current - previous)));
        }
		i++;
	}
}

void readCoord(std::string srcoord, param &par);
void readParam(int argc, char **argv, param &par);
void readParam(std::string parfile, param &par);
void readParameters(int argc, char **argv, param &par);

void analyzeNLInversion(param &par);
void analyzeBsplines(const hypercube<data_t> &domain, param &par);
