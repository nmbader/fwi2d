#pragma once

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "vecReg.hpp"
#include "misc.hpp"
#include "param.hpp"
#include "seplib.h"

// write to SEP format (.H file)
template <typename T>
void sepWrite(const std::shared_ptr<vecReg<T> > vec, std::string output){
    
    long n123 = vec->getHyper()->getN123();
    std::vector<axis<T> > axes = vec->getHyper()->getAxes();
    std::string ni, oi, di, labeli, esize="esize";
    float of, df;
    const char* out=output.c_str();
    int size=sizeof(T);
    for (int i=0; i<axes.size(); i++){
        ni = "n" + std::to_string(i+1);
        oi = "o" + std::to_string(i+1);
        di = "d" + std::to_string(i+1);
        labeli = "label" + std::to_string(i+1);
        of = (float)axes[i].o;
        df = (float)axes[i].d;

        if (output == "out"){
            putch((char*)esize.c_str(),(char*)"d",(int*)&size);
            putch((char*)ni.c_str(),(char*)"d",(int*)&axes[i].n);
            putch((char*)oi.c_str(),(char*)"f",(void*)&of);
            putch((char*)di.c_str(),(char*)"f",(void*)&df);
            putch((char*)labeli.c_str(),(char*)"s",(void*)axes[i].label.c_str());
        }
        else{
            auxputch(esize.c_str(),"d",&size,out);
            auxputch(ni.c_str(),"d",&axes[i].n,out);
            auxputch(oi.c_str(),"f",&of,out);
            auxputch(di.c_str(),"f",&df,out);
            auxputch(labeli.c_str(),"s",axes[i].label.c_str(),out);
        }
    }
    successCheck(size*n123 == srite_big(out,vec->getVals(),size*n123),__FILE__,__LINE__, "Cannot write data\n");
}


// read from SEP format (.H file)
template <typename T>
std::shared_ptr<vecReg<T> > sepRead(std::string input, int ndim0=1){
    
    int n1=1, n2=1, n3=1, n4=1, n5=1, n6=1, size=4;
	float o1=0, o2=0, o3=0, o4=0, o5=0, o6=0, d1=1, d2=1, d3=1, d4=1, d5=1, d6=1;
	char label1[20]={""}, label2[20]={""}, label3[20]={""}, label4[20]={""}, label5[20]={""}, label6[20]={""};
    const char* data = input.c_str();

    if (input == "in"){
        hetch((char*)"esize",(char*)"d",&size);
        hetch((char*)"n1",(char*)"d",&n1);
        hetch((char*)"n2",(char*)"d",&n2);
        hetch((char*)"n3",(char*)"d",&n3);
        hetch((char*)"n4",(char*)"d",&n4);
        hetch((char*)"n5",(char*)"d",&n5);
        hetch((char*)"n6",(char*)"d",&n6);
        hetch((char*)"o1",(char*)"f",&o1);
        hetch((char*)"o2",(char*)"f",&o2);
        hetch((char*)"o3",(char*)"f",&o3);
        hetch((char*)"o4",(char*)"f",&o4);
        hetch((char*)"o5",(char*)"f",&o5);
        hetch((char*)"o6",(char*)"f",&o6);
        hetch((char*)"d1",(char*)"f",&d1);
        hetch((char*)"d2",(char*)"f",&d2);
        hetch((char*)"d3",(char*)"f",&d3);
        hetch((char*)"d4",(char*)"f",&d4);
        hetch((char*)"d5",(char*)"f",&d5);
        hetch((char*)"d6",(char*)"f",&d6);
        hetch((char*)"label1",(char*)"s",&label1);
        hetch((char*)"label2",(char*)"s",&label2);
        hetch((char*)"label3",(char*)"s",&label3);
        hetch((char*)"label4",(char*)"s",&label4);
        hetch((char*)"label5",(char*)"s",&label5);
        hetch((char*)"label6",(char*)"s",&label6);
    }
    else{
        auxpar((char*)"esize",(char*)"d",&size,data);
        auxpar((char*)"n1",(char*)"d",&n1,data);
        auxpar((char*)"n2",(char*)"d",&n2,data);
        auxpar((char*)"n3",(char*)"d",&n3,data);
        auxpar((char*)"n4",(char*)"d",&n4,data);
        auxpar((char*)"n5",(char*)"d",&n5,data);
        auxpar((char*)"n6",(char*)"d",&n6,data);
        auxpar((char*)"o1",(char*)"f",&o1,data);
        auxpar((char*)"o2",(char*)"f",&o2,data);
        auxpar((char*)"o3",(char*)"f",&o3,data);
        auxpar((char*)"o4",(char*)"f",&o4,data);
        auxpar((char*)"o5",(char*)"f",&o5,data);
        auxpar((char*)"o6",(char*)"f",&o6,data);
        auxpar((char*)"d1",(char*)"f",&d1,data);
        auxpar((char*)"d2",(char*)"f",&d2,data);
        auxpar((char*)"d3",(char*)"f",&d3,data);
        auxpar((char*)"d4",(char*)"f",&d4,data);
        auxpar((char*)"d5",(char*)"f",&d5,data);
        auxpar((char*)"d6",(char*)"f",&d6,data);
        auxpar((char*)"label1",(char*)"s",&label1,data);
        auxpar((char*)"label2",(char*)"s",&label2,data);
        auxpar((char*)"label3",(char*)"s",&label3,data);
        auxpar((char*)"label4",(char*)"s",&label4,data);
        auxpar((char*)"label5",(char*)"s",&label5,data);
        auxpar((char*)"label6",(char*)"s",&label6,data);
    }
	
	std::string lab1(&label1[0]);
	std::string lab2(&label2[0]);
    std::string lab3(&label3[0]);
	std::string lab4(&label4[0]);
	std::string lab5(&label5[0]);
	std::string lab6(&label6[0]);

	axis<T> ax1(n1, (T)o1, (T)d1, lab1);
	axis<T> ax2(n2, (T)o2, (T)d2, lab2);
    axis<T> ax3(n3, (T)o3, (T)d3, lab3);
    axis<T> ax4(n4, (T)o4, (T)d4, lab4);
    axis<T> ax5(n5, (T)o5, (T)d5, lab5);
    axis<T> ax6(n6, (T)o6, (T)d6, lab6);

	long n123 = (long)n1*(long)n2*(long)n3*(long)n4*(long)n5*(long)n6;

    std::vector<axis<T> > axes = {ax1, ax2, ax3, ax4, ax5, ax6};
    int ndim = 6;
    if ((n6==1) && (ndim0!=6)){
        ndim -=1;
        axes.pop_back();
        if ((n5==1) && (ndim0!=5)){
            ndim -=1;
            axes.pop_back();
            if ((n4==1) && (ndim0!=4)){
                ndim -= 1;
                axes.pop_back();
                if ((n3==1) && (ndim0!=3)){
                    ndim -= 1;
                    axes.pop_back();
                    if ((n2==1) && (ndim0!=2)){
                        ndim -= 1;
                        axes.pop_back();
                    }
                }
            }
        }
    }

    hypercube<T> hyper(axes);
    std::shared_ptr<vecReg<T> > vec = std::make_shared<vecReg<T> > (hyper);

    if (sizeof(T)==size){
        successCheck(size*n123 == sreed_big(data,vec->getVals(),size*n123),__FILE__,__LINE__, "Cannot read data\n");
    }
    else {
        if (size==4){
            float* val = new float[n123];
            successCheck(size*n123 == sreed_big(data,val,size*n123),__FILE__,__LINE__, "Cannot read data\n");
            for (long i=0; i<n123; i++) vec->getVals()[i]=val[i];
            delete [] val;
        }
        else if (size==8){
            double* val = new double[n123];
            successCheck(size*n123 == sreed_big(data,val,size*n123),__FILE__,__LINE__, "Cannot read data\n");
            for (long i=0; i<n123; i++) vec->getVals()[i]=val[i];
            delete [] val;
        }
        else successCheck(false,__FILE__,__LINE__,"Cannot read into the current template vector.\n");
    }
    return vec;
}



// write to binary with a text file following SEP convention
template <typename T>
void binWrite(const std::shared_ptr<vecReg<T> > vec, std::string output, std::string datapath="none"){
    
    long n123 = vec->getHyper()->getN123();
    std::vector<axis<T> > axes = vec->getHyper()->getAxes();
    std::string ni, oi, di, labeli;

    std::ofstream output_file;
    output_file.open(output);
    std::string msg = "The output file "+output+" could not be created. Check for the path\n";
    successCheck(output_file.is_open(),__FILE__,__LINE__,msg.c_str());
    for (int i=0; i<axes.size(); i++){
        ni = "n" + std::to_string(i+1);
        oi = "o" + std::to_string(i+1);
        di = "d" + std::to_string(i+1);
        labeli = "label" + std::to_string(i+1);

        output_file << ni << "=" << axes[i].n << "\n";
        output_file << oi << "=" << axes[i].o << "\n";
        output_file << di << "=" << axes[i].d << "\n";
    }
    output_file << "esize=" << sizeof(T) <<"\n";

    char* pPath;
    std::string bin_output;
    if (datapath=="none"){
        pPath = getenv ("DATAPATH");
        if (pPath==NULL) pPath = getenv("PWD");
        if (pPath==NULL) pPath = ".";
        bin_output = std::string(pPath);
    }
    else bin_output = datapath;
    bin_output=bin_output+"/"+output+"@";
    output_file << "in=" << bin_output <<"\n";
    output_file.close();

    FILE * pFile;
    T * pdata = vec->getVals();
    pFile = fopen (bin_output.c_str(), "wb");
    msg = "The output file "+bin_output+" could not be created. Check the path\n";
    successCheck(pFile!=NULL,__FILE__,__LINE__,msg.c_str());
    fwrite (pdata , sizeof(T), n123, pFile);
    fclose (pFile);
}


// read from a binary with a text file following SEP convention
template <typename T>
std::shared_ptr<vecReg<T> > binRead(std::string input, int ndim0=1){
    
    int n1=1, n2=1, n3=1, n4=1, n5=1, n6=1, esize=4;
	T o1=0, o2=0, o3=0, o4=0, o5=0, o6=0, d1=1, d2=1, d3=1, d4=1, d5=1, d6=1;
	std::string label1="", label2="", label3="", label4="", label5="", label6="";
    std::string in="";

    std::ifstream ifs (input, std::ifstream::in);
    std::string msg = "The input file "+input+" could not be opened. Check the path or the file name\n";
    successCheck(ifs.is_open(),__FILE__,__LINE__,msg.c_str());

    std::string line;
    while ( getline (ifs, line) )
    {
        if (line[0] != '#')
        {
            std::istringstream ss(line);
            std::string word;
            while (ss >> word)
            {
                if (word.substr(0, 3) == "n1=") {
                    word.erase(0, 3);
                    n1=convert_to<int>(word);
                }
                else if (word.substr(0, 3) == "n2=") {
                    word.erase(0, 3);
                    n2=convert_to<int>(word);
                }
                else if (word.substr(0, 3) == "n3=") {
                    word.erase(0, 3);
                    n3=convert_to<int>(word);
                }
                else if (word.substr(0, 3) == "n4=") {
                    word.erase(0, 3);
                    n4=convert_to<int>(word);
                }
                else if (word.substr(0, 3) == "n5=") {
                    word.erase(0, 3);
                    n5=convert_to<int>(word);
                }
                else if (word.substr(0, 3) == "n6=") {
                    word.erase(0, 3);
                    n6=convert_to<int>(word);
                }
                else if (word.substr(0, 3) == "o1=") {
                    word.erase(0, 3);
                    o1=convert_to<T>(word);
                }
                else if (word.substr(0, 3) == "o2=") {
                    word.erase(0, 3);
                    o2=convert_to<T>(word);
                }
                else if (word.substr(0, 3) == "o3=") {
                    word.erase(0, 3);
                    o3=convert_to<T>(word);
                }
                else if (word.substr(0, 3) == "o4=") {
                    word.erase(0, 3);
                    o4=convert_to<T>(word);
                }
                else if (word.substr(0, 3) == "o5=") {
                    word.erase(0, 3);
                    o5=convert_to<T>(word);
                }
                else if (word.substr(0, 3) == "o6=") {
                    word.erase(0, 3);
                    o6=convert_to<T>(word);
                }
                else if (word.substr(0, 3) == "d1=") {
                    word.erase(0, 3);
                    d1=convert_to<T>(word);
                }
                else if (word.substr(0, 3) == "d2=") {
                    word.erase(0, 3);
                    d2=convert_to<T>(word);
                }
                else if (word.substr(0, 3) == "d3=") {
                    word.erase(0, 3);
                    d3=convert_to<T>(word);
                }
                else if (word.substr(0, 3) == "d4=") {
                    word.erase(0, 3);
                    d4=convert_to<T>(word);
                }
                else if (word.substr(0, 3) == "d5=") {
                    word.erase(0, 3);
                    d5=convert_to<T>(word);
                }
                else if (word.substr(0, 3) == "d6=") {
                    word.erase(0, 3);
                    d6=convert_to<T>(word);
                }
                else if (word.substr(0, 7) == "label1=") {
                    word.erase(0, 7);
                    label1=word;
                }
                else if (word.substr(0, 7) == "label2=") {
                    word.erase(0, 7);
                    label2=word;
                }
                else if (word.substr(0, 7) == "label3=") {
                    word.erase(0, 7);
                    label3=word;
                }
                else if (word.substr(0, 7) == "label4=") {
                    word.erase(0, 7);
                    label4=word;
                }
                else if (word.substr(0, 7) == "label5=") {
                    word.erase(0, 7);
                    label5=word;
                }
                else if (word.substr(0, 7) == "label6=") {
                    word.erase(0, 7);
                    label6=word;
                }
                else if (word.substr(0, 6) == "esize=") {
                    word.erase(0, 6);
                    esize=convert_to<int>(word);
                }
                else if (word.substr(0, 3) == "in=") {
                    word.erase(0, 3);
                    if (word.at(0)=='"') {word.erase(0,1); word.pop_back();}
                    in=word;
                }
            }
        }
    }
    ifs.close();

    if (esize != sizeof(T)) fprintf(stderr,"\nWARNING: 'esize' is inconsistent with the data type.\n");

  	axis<T> ax1(n1, o1, d1, label1);
	axis<T> ax2(n2, o2, d2, label2);
    axis<T> ax3(n3, o3, d3, label3);
    axis<T> ax4(n4, o4, d4, label4);
    axis<T> ax5(n5, o5, d5, label5);
    axis<T> ax6(n6, o6, d6, label6);

	long n123 = (long)n1*(long)n2*(long)n3*(long)n4*(long)n5*(long)n6;
    std::vector<axis<T> > axes = {ax1, ax2, ax3, ax4, ax5, ax6};
    int ndim = 6;
    if ((n6==1) && (ndim0!=6)){
        ndim -=1;
        axes.pop_back();
        if ((n5==1) && (ndim0!=5)){
            ndim -=1;
            axes.pop_back();
            if ((n4==1) && (ndim0!=4)){
                ndim -= 1;
                axes.pop_back();
                if ((n3==1) && (ndim0!=3)){
                    ndim -= 1;
                    axes.pop_back();
                    if ((n2==1) && (ndim0!=2)){
                        ndim -= 1;
                        axes.pop_back();
                    }
                }
            }
        }
    }

    hypercube<T> hyper(axes);
    std::shared_ptr<vecReg<T> > vec = std::make_shared<vecReg<T> > (hyper);
    T * pdata = vec->getVals();

    FILE * pFile;
    pFile = fopen(in.c_str(), "rb");
    msg = "The input file "+in+" could not be read. Check the path or the file name\n";
    successCheck(pFile!=NULL,__FILE__,__LINE__,msg.c_str());

    // obtain file size in bytes:
    fseek (pFile , 0 , SEEK_END);
    long lSize = ftell (pFile);
    rewind (pFile);

    msg = "The binary file contains "+std::to_string(lSize)+" bytes whereas the description file suggests "+std::to_string(n123*esize)+" bytes\n";
    successCheck(lSize==n123*esize,__FILE__,__LINE__,msg.c_str());

    long result;
    if (sizeof(T)==esize) result = fread(pdata, esize, n123, pFile);
    else{
        if (esize==4){
            float* val = new float[n123];
            result = fread(val, esize, n123, pFile);
            for (long i=0; i<n123; i++) pdata[i]=val[i];
            delete [] val;
        }
        else if (esize==8){
            double* val = new double[n123];
            result = fread(val, esize, n123, pFile);
            for (long i=0; i<n123; i++) pdata[i]=val[i];
            delete [] val;
        }
        else successCheck(false,__FILE__,__LINE__,"Cannot read into the current template vector.\n");
    }

    msg = "The binary file was not read properly. It may be that the requested data type and esize are inconsistent\n";
    successCheck(result==n123,__FILE__,__LINE__,msg.c_str());

    return vec;
}

// wrappers write and read functions to choose among sepWrite/sepRead and binWrite/binRead
template <typename T>
void write(const std::shared_ptr<vecReg<T> > vec, std::string output, bool format=false, std::string datapath="none"){
    if (format==false) sepWrite<T>(vec, output);
    else binWrite<T>(vec,output,datapath); 
}

template <typename T>
std::shared_ptr<vecReg<T> > read(std::string input, bool format=false, int ndim0=1){
    if (format==0) return sepRead<T>(input, ndim0);
    else return binRead<T>(input, ndim0);
}