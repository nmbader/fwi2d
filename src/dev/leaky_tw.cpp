#ifdef DOUBLE_PRECISION
    typedef double data_t;
#else
    typedef float data_t;
#endif

#include <string.h>
#include "we_op.hpp"
#include "IO.hpp"
#include "seplib.h"


typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;

// Executable to model 2D seismic data and optionally save the full wavefield

int main(int argc, char **argv){

    int rank=0, size=0;
#ifdef ENABLE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    fprintf (stderr,"\n====================\nSize of MPI communicator = %d ; current rank = %d\n====================\n",size,rank);
#endif

    initpar(argc,argv);

    // Read parameters for wave propagation
    param par;
    readParameters(argc, argv, par);
    if (rank>0) par.verbose=0;
    par.device+=rank;

    // Read inputs/outputs files
    std::string source_file="none", model_file="none", transmission_file="none", wavefield_file="none", output_file="none";
    readParam<std::string>(argc, argv, "source", source_file);
    readParam<std::string>(argc, argv, "model", model_file);
    readParam<std::string>(argc, argv, "transmission", transmission_file);
    readParam<std::string>(argc, argv, "wavefield", wavefield_file);
    readParam<std::string>(argc, argv, "output", output_file);

    successCheck(source_file!="none",__FILE__,__LINE__,"Source wavelet is not provided\n");
    successCheck(model_file!="none",__FILE__,__LINE__,"Earth model is not provided\n");

    std::shared_ptr<vec> src = read<data_t>(source_file, par.format);
    std::shared_ptr<vec> model = read<data_t>(model_file, par.format);
    
    successCheck(model->getHyper()->getNdim()==3,__FILE__,__LINE__,"The model must have 3 axes\n");
    successCheck(model->getHyper()->getAxis(3).n==2,__FILE__,__LINE__,"The model must constain exactly 2 parameters: Vp and density\n");

    std::shared_ptr<vec> transmission;
    if (transmission_file=="none") {
        if (par.verbose>0 && par.bc_top==4) fprintf(stderr,"\nThe top boundary transmission coefficients are assumed to be = 0 (rigid-wall BC)\n");
        if (par.verbose>0 && par.bc_top==5) fprintf(stderr,"\nThe top boundary mask is assumed to be = 0 (rigid-wall BC)\n");
        transmission=std::make_shared<vec> (hyper(model->getHyper()->getAxis(2)));
        transmission->zero();
    }
    else {
        if (par.verbose>0) fprintf(stderr,"\nThe top boundary transmission coefficients (or mask) are read from provided file\n");
        transmission = read<data_t>(transmission_file, par.format);
        successCheck(transmission->min()>=0,__FILE__,__LINE__,"The transmission coefficients (or mask) of the top boundary must be non-negative\n");
        //successCheck(transmission->max()<=1,__FILE__,__LINE__,"The transmission coefficients (or mask) must be at most = 1\n");
    }
    
    // Analyze the input source time function and duplicate if necessary, analyze geometry
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0); 

    // Build the appropriate wave equation operator
    nl_we_op_a * op = new nl_we_op_a(*model->getHyper(),allsrc,par);
    op->_transmission = transmission;
    
    // Run the forward modeling
    std::shared_ptr<vec> allrcv = std::make_shared<vec> (*op->getRange());
    op->forward(false,model,allrcv);

    if ((rank==0) && (wavefield_file!="none") && (op->_par.sub)>0) write<data_t>(op->_full_wfld, wavefield_file, par.format, par.datapath);
    if ((rank==0) && (output_file!="none")) write<data_t>(allrcv, output_file, par.format, par.datapath);

    delete op;

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif

return 0;
}