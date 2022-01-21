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
    std::string source_file="none", model_file="none", wavefield_file="none", output_file="none";
    readParam<std::string>(argc, argv, "source", source_file);
    readParam<std::string>(argc, argv, "model", model_file);
    readParam<std::string>(argc, argv, "wavefield", wavefield_file);
    readParam<std::string>(argc, argv, "output", output_file);

    successCheck(source_file!="none",__FILE__,__LINE__,"Source wavelet is not provided\n");
    successCheck(model_file!="none",__FILE__,__LINE__,"Earth model is not provided\n");

    std::shared_ptr<vec> src = read<data_t>(source_file, par.format);
    std::shared_ptr<vec> model = read<data_t>(model_file, par.format);
    
    // Analyze the input source time function and duplicate if necessary, analyze geometry
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0);
    //analyzeModel(*allsrc->getHyper(),model,par);
    //l_we_op_ae lop(*allsrc->getHyper(),model,par);
    //lop.dotProduct();
    
    // If more than one shot is modeled, don't save the wavefield
    // if (par.ns>1) par.sub=0;

    // Build the appropriate wave equation operator
    nl_we_op_e * op;
    if (par.nmodels==2) op=new nl_we_op_a(*model->getHyper(),allsrc,par);
    else if (par.nmodels==3 && !par.acoustic_elastic) op=new nl_we_op_e(*model->getHyper(),allsrc,par);
    else if (par.nmodels==3 && par.acoustic_elastic) op=new nl_we_op_ae(*model->getHyper(),allsrc,par);
    else if (par.nmodels==5) op=new nl_we_op_vti(*model->getHyper(),allsrc,par);

    // Run the forward modeling
    std::shared_ptr<vec> allrcv = std::make_shared<vec> (*op->getRange());
    op->forward(false,model,allrcv);

    if ((rank==0) && (wavefield_file!="none") && (op->_par.sub)>0) write<data_t>(op->_full_wfld, wavefield_file, par.format, par.datapath);
    if ((rank==0) && (wavefield_file!="none") && (op->_par.sub)>0 && par.acoustic_elastic && par.acoustic_wavefield) write<data_t>(op->_full_wflda, wavefield_file+"a", par.format, par.datapath);
    if ((rank==0) && (output_file!="none")) write<data_t>(allrcv, output_file, par.format, par.datapath);

    delete op;

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif

return 0;
}