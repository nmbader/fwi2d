#ifdef DOUBLE_PRECISION
    typedef double data_t;
#else
    typedef float data_t;
#endif

#include <unistd.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "vecReg.hpp"
#include "IO.hpp"

typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;

void printdoc(){
    std::string doc = "\nDescription:\n"
    "   Convert SEPlib data to/from binary with a description file similar to SEPlib .H files.\n"
    "\nInput/Output:\n"
    "   Provide input as 'input=file.H' and output as 'output=file.H' where 'file.H' is a description file.\n"
    "\nParameters:\n"
    "   mode - int - [0]:\n\t\t0: convert sep file to binary. otherwise: convert binary file to sep format.\n"
    "   datapath - string - ['none']:\n\t\tdatapath to save binary data in mode 0. If 'none' it will default to DATAPATH or PWD environment vavriables.\n"
    "\nExample:\n"
    "   SEP2BIN.x input=seplib_data.H output=description_file.H mode=0 datapath=/path/to/binaries\n"
    "   SEP2BIN.x input=description_file.H output=seplib_data.H mode=1\n"
    "\n";
    fprintf(stderr,doc.c_str());
}

int main(int argc, char **argv){

    if (argc == 1 && isatty(STDIN_FILENO)==true) {printdoc(); return 0;}

	initpar(argc,argv);

    int mode=0;
    std::string input_file="none", output_file="none", datapath="none";
    
    readParam<std::string>(argc, argv, "input", input_file);
	readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "datapath", datapath);
    readParam<int>(argc, argv, "mode", mode);
    
    successCheck(input_file!="none",__FILE__,__LINE__,"Input file is not provided\n");

    if (mode==0)
    {
        std::shared_ptr<vec> v = sepRead<data_t>(input_file);
        if (output_file!="none") binWrite<data_t>(v, output_file, datapath);
    }

    else
    {
        std::shared_ptr<vec> v = binRead<data_t>(input_file);
        if (output_file!="none") sepWrite<data_t>(v, output_file);
    }

    return 0;
}