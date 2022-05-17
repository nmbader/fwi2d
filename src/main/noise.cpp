#include <unistd.h>
#include <ctime>
#include <random>

#include "param.hpp"
#include "IO.hpp"
#include "seplib.h"

#include <iostream>
#include <fstream>

typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;

void printdoc(){
    std::string doc = "\nDescription:\n"
    "   Generate noise and add it (or replace) to the input.\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' or '< file.H' and output as 'output=file.H' or '> file.H'. '<' and '>' are valid for SEPlib format only.\n"
    "\nParameters:\n"
    "   type - string - ['uniform']:\n\t\toptions: 'uniform', 'normal'.\n"
    "   mean - float - [0]:\n\t\tmean value for the normal distribution.\n"
    "   sigma - float - [1]:\n\t\tstandard deviation for the normal distribution.\n"
    "   min - float - [-1]:\n\t\tminimum value for the uniform distribution.\n"
    "   max - float - [1]:\n\t\tmaximum value for the uniform distribution.\n"
    "   seed - float - [date]:\n\t\tseed for random numbers generation.\n"
    "   replace - int - [0]:\n\t\t0 = add noise to input, 1 = replace input.\n"
    "   format - bool - [0]:\n\t\tdata format for IO. 0 for SEPlib, 1 for binary with description file.\n"
    "   datapath - string - ['none']:\n\t\tpath for output binaries.\n"
    "\nExample:\n"
    "   NOISE.x < infile.H type=uniform min=-10 max=20 > oufile.H.\n"
    "   NOISE.x < infile.H type=normal mean=3 sigma=5 replace=1 > oufile.H\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);
    omp_set_num_threads(1);

    std::string input_file="in", output_file="out", type="uniform", datapath="none";
    data_t mean=0, sigma=1, min=-1, max=1, seed=0;
    bool replace=0, format=0;
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "type", type);
	readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<data_t>(argc, argv, "mean", mean);
    readParam<data_t>(argc, argv, "sigma", sigma);
    readParam<data_t>(argc, argv, "min", min);
    readParam<data_t>(argc, argv, "max", max);
    readParam<data_t>(argc, argv, "seed", seed);
    readParam<bool>(argc, argv, "replace", replace);
    readParam<bool>(argc, argv, "format", format);

    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");
    
    std::shared_ptr<vec> input = read<data_t>(input_file, format);

    if (seed == 0) seed = (data_t) time(0);

    std::default_random_engine generator;
    generator.seed(seed);
    std::uniform_real_distribution<data_t> uniform(min,max);
    std::normal_distribution<data_t> normal(mean,sigma);

    data_t * pin = input->getVals();
    int n = input->getN123();
    data_t a = 1-replace;

    if (type=="uniform") {
        for (int i=0; i<n; i++) pin[i] = a*pin[i] + uniform(generator);
    }
    else if (type=="normal"){
        for (int i=0; i<n; i++) pin[i] = a*pin[i] + normal(generator);
    }
    else{
        successCheck(false, __FILE__,__LINE__,"Type is not implemented\n");
    }
 
    if (output_file!="none") write<data_t>(input, output_file, format, datapath);

    return 0;
}