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
    std::vector<data_t> rdip;
    data_t gl=0;

    // source mechanism and receiver type
    bool mt=false;
    data_t fangle=0;
    data_t mxx=1, mzz=1, mxz=1;
    int seismotype=0;

    // boundary parameters
    int bc_top=1, bc_bottom=1, bc_left=1, bc_right=1, taper_top=0, taper_bottom=0, taper_left=0, taper_right=0;
    data_t alpha=0.2508560249/1.05, taper_strength=0.05;
    bool pml=false, pml_T=false, pml_B=false, pml_L=false, pml_R=false;
    data_t R = 1e-3; int p=2; 

    // model bounds
    data_t vpmin=0.2, vpmax=7, vsmin=0.2, vsmax=7, rhomin=0, rhomax=7, deltamin=-0.5, deltamax=1, epsilonmin=0, epsilonmax=1;
    data_t vmax=7, vmin=0.3;

    // model parameterization
    int nmodels=3;
    bool bsplines=false;
    int bs_nx=3, bs_nz=3;
    std::vector<int> bs_mx, bs_mz;
    std::vector<data_t> bs_controlx, bs_controlz;

    // inversion parameters
    std::string nlsolver="lbfgs", lsearch="regular_wolfe", mask_file="none", weights_file="none";
    data_t threshold=0;
    int niter=0, max_trial=10, isave=10, envelop=0;
    bool solver_verbose=true, normalize=0;

    // miscallenous
    int version=2;
    std::vector<int> device;
    bool verbose=false;
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
