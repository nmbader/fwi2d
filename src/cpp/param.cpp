#include "param.hpp"

void readCoord(std::string srcoord, param &par){

    std::ifstream ifs (srcoord, std::ifstream::in);
    if (ifs.is_open())
    {
        par.sxz.clear();
        par.rxz.clear();
        par.srcoord_from_file=true;
        std::string line;
        int counter=0;
        int shot_counter=0;
        int curr=0, prev=0;
        while ( getline (ifs, line) )
        {
            if (line[0] != '#')
            {
                std::istringstream ss(line);

                ss >> curr; // read shot id
                if (counter==0) // first line
                {
                    // read shot x z coordinates
                    std::vector<data_t> sxz(2);
                    ss >> sxz[0];
                    ss >> sxz[1];

                    // read receiver x z dip coordinates
                    std::vector<data_t> rxz(3);
                    ss >> rxz[0];
                    ss >> rxz[1];
                    ss >> rxz[2];

                    par.sxz.push_back(sxz);
                    std::vector<std::vector<data_t> > vec = {rxz};
                    par.rxz.push_back(vec);

                    counter++;
                    prev=curr;
                }
                else // subsequent lines
                {
                    // read shot x z coordinates
                    std::vector<data_t> sxz(2);
                    ss >> sxz[0];
                    ss >> sxz[1];

                    // read receiver x z dip coordinates
                    std::vector<data_t> rxz(3);
                    ss >> rxz[0];
                    ss >> rxz[1];
                    ss >> rxz[2];

                    if (curr==prev)
                    {
                        par.rxz[shot_counter].push_back(rxz);
                    }
                    else
                    {   
                        par.sxz.push_back(sxz);
                        std::vector<std::vector<data_t> > vec = {rxz};
                        par.rxz.push_back(vec);
                        shot_counter++;
                    }
                    counter++;
                    prev=curr;

                }
            }
        }
        ifs.close();
    }
}

void readParam(int argc, char **argv, param &par){
    readParam<std::string>(argc, argv, "resampling", par.resampling);
    readParam<std::string>(argc, argv, "srcoord", par.srcoord);
    readParam<std::string>(argc, argv, "nlsolver", par.nlsolver);
    readParam<std::string>(argc, argv, "lsearch", par.lsearch);
    readParam<data_t>(argc, argv, "courant", par.courant);
    readParam<data_t>(argc, argv, "dt", par.dt);
    readParam<data_t>(argc, argv, "sx0", par.sx0);
    readParam<data_t>(argc, argv, "sz0", par.sz0);
    readParam<data_t>(argc, argv, "rx0", par.rx0);
    readParam<data_t>(argc, argv, "rz0", par.rz0);
    readParam<data_t>(argc, argv, "rdip0", par.rdip0);
    readParam<data_t>(argc, argv, "sxinc", par.sxinc);
    readParam<data_t>(argc, argv, "szinc", par.szinc);
    readParam<data_t>(argc, argv, "rxinc", par.rxinc);
    readParam<data_t>(argc, argv, "rzinc", par.rzinc);
    readParam<data_t>(argc, argv, "rdipinc", par.rdipinc);
    readParam<data_t>(argc, argv, "alpha", par.alpha);
    readParam<data_t>(argc, argv, "taper_strength", par.taper_strength);
    readParam<data_t>(argc, argv, "R", par.R);
    readParam<data_t>(argc, argv, "gl", par.gl);
    readParam<data_t>(argc, argv, "mxx", par.mxx);
    readParam<data_t>(argc, argv, "mzz", par.mzz);
    readParam<data_t>(argc, argv, "mxz", par.mxz);
    readParam<data_t>(argc, argv, "fangle", par.fangle);
    readParam<data_t>(argc, argv, "vpmin", par.vpmin);
    readParam<data_t>(argc, argv, "vpmax", par.vpmax);
    readParam<data_t>(argc, argv, "vsmin", par.vsmin);
    readParam<data_t>(argc, argv, "vsmax", par.vsmax);
    readParam<data_t>(argc, argv, "rhomin", par.rhomin);
    readParam<data_t>(argc, argv, "rhomax", par.rhomax);
    readParam<data_t>(argc, argv, "threshold", par.threshold);
    readParam<int>(argc, argv, "ns", par.ns);
    readParam<int>(argc, argv, "nr", par.nr);
    readParam<int>(argc, argv, "seismotype", par.seismotype);
    readParam<int>(argc, argv, "bc_top", par.bc_top);
    readParam<int>(argc, argv, "bc_bottom", par.bc_bottom);
    readParam<int>(argc, argv, "bc_left", par.bc_left);
    readParam<int>(argc, argv, "bc_right", par.bc_right);
    readParam<int>(argc, argv, "taper_top", par.taper_top);
    readParam<int>(argc, argv, "taper_bottom", par.taper_bottom);
    readParam<int>(argc, argv, "taper_left", par.taper_left);
    readParam<int>(argc, argv, "taper_right", par.taper_right);
    readParam<int>(argc, argv, "p", par.p);
    readParam<int>(argc, argv, "sinc_half_length", par.sinc_half_length);
    readParam<int>(argc, argv, "sub", par.sub);
    readParam<int>(argc, argv, "niter", par.niter);
    readParam<int>(argc, argv, "max_trial", par.max_trial);
    readParam<int>(argc, argv, "isave", par.isave);
    readParam<int>(argc, argv, "version", par.version);
    readParam<bool>(argc, argv, "mt", par.mt);
    readParam<bool>(argc, argv, "pml", par.pml);
    readParam<bool>(argc, argv, "verbose", par.verbose);
    readParam<bool>(argc, argv, "solver_verbose", par.solver_verbose);
    readParam<int>(argc, argv, "device", par.device);
}

void readParam(std::string parfile, param &par){
    
    std::ifstream ifs (parfile, std::ifstream::in);
    if (ifs.is_open())
    {
        std::string line;
        std::vector<std::string> lines = {"empty_line"};
        while ( getline (ifs, line) )
        {
            line = line + " ";
            if (line [0]!= '#')
            {
                lines.push_back(line);
            }
        }
        int argc=lines.size();
        char * argv[argc];
        for (int i=0; i<argc; i++){
            argv[i] = (char *) (lines[i].c_str());
        }
        readParam(argc, argv, par);
        ifs.close();
    }
}

void readParameters(int argc, char **argv, param &par){
    std::string parfile="none";
    readParam<std::string>(argc, argv, "parfile", parfile);
    readParam(parfile, par);
    readParam(argc, argv, par);
    readCoord(par.srcoord, par);
    if (par.srcoord_from_file==true) {
        par.ns=par.sxz.size();
    }
    else {
        par.sxz.clear();
        par.rxz.clear();
        for (int s=0; s<par.ns; s++){
            // compute shot x z coordinates
            std::vector<data_t> sxz = {par.sx0+s*par.sxinc, par.sz0+s*par.szinc};

            // compute receiver x z dip coordinates
            std::vector<std::vector<data_t> > vec;
            for (int r=0; r<par.nr; r++){
                std::vector<data_t> rxz = {par.rx0+r*par.rxinc, par.rz0+r*par.rzinc, par.rdip0+r*par.rdipinc};
                vec.push_back(rxz);
            }
            par.sxz.push_back(sxz);
            par.rxz.push_back(vec);
        }
    }
}

void analyzeNLInversion(param &par)
{
    fprintf(stderr,"\n==========================\n Inversion parameters\n==========================\n");
    par.isave=std::max(1,par.isave);
    fprintf(stderr,"Non-linear solver = %s\n",par.nlsolver.c_str());
    fprintf(stderr,"Line search = %s\n",par.lsearch.c_str());
    fprintf(stderr,"Number of iterations = %d\n",par.niter);
    fprintf(stderr,"Maximum number of trials per iteration = %d\n",par.max_trial);
    fprintf(stderr,"Threshold to stop the inversion = %f\n",par.threshold);
    fprintf(stderr,"If \"ioutput\" is provided, the model, gradient, and residual will be saved every %d iterations\n",par.isave);
}