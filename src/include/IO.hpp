#pragma once

#include <string>
#include <sstream>
#include "vecReg.hpp"
#include "misc.hpp"
#include "seplib.h"

// write to SEP format (.H file)
template <typename T>
void sepWrite(const std::shared_ptr<vecReg<T> > vec, std::string output){
    
    long long n123 = vec->getHyper()->getN123();
    std::vector<axis<T> > axes = vec->getHyper()->getAxes();
    std::string ni, oi, di, labeli;
    float of, df;
    const char* out=output.c_str();
    for (int i=0; i<axes.size(); i++){
        ni = "n" + std::to_string(i+1);
        oi = "o" + std::to_string(i+1);
        di = "d" + std::to_string(i+1);
        labeli = "label" + std::to_string(i+1);
        of = (float)axes[i].o;
        df = (float)axes[i].d;

        if (output == "out"){
            putch((char*)ni.c_str(),(char*)"d",(int*)&axes[i].n);
            putch((char*)oi.c_str(),(char*)"f",(void*)&of);
            putch((char*)di.c_str(),(char*)"f",(void*)&df);
            putch((char*)labeli.c_str(),(char*)"s",(void*)axes[i].label.c_str());
        }
        else{
            auxputch(ni.c_str(),"d",&axes[i].n,out);
            auxputch(oi.c_str(),"f",&of,out);
            auxputch(di.c_str(),"f",&df,out);
            auxputch(labeli.c_str(),"s",axes[i].label.c_str(),out);
        }
    }
    
    if (sizeof(T)==4){
        successCheck(4*n123 == srite_big(out,vec->getVals(),4*n123), __FILE__,__LINE__,"Cannot write data\n");
    }
    else if (sizeof(T)==8){
        float* val = new float[n123];
        for (int i=0; i<n123; i++) val[i] = vec->getVals()[i];
        successCheck(4*n123 == srite_big(out,val,4*n123), __FILE__,__LINE__,"Cannot write data\n");
        delete val;
    }
    else{
        successCheck(false, __FILE__,__LINE__,"Cannot write the current template vector.\n");
    }
}


// read from SEP format (.H file)
template <typename T>
std::shared_ptr<vecReg<T> > sepRead(std::string input, int ndim0=1){
    
    int n1=1, n2=1, n3=1, n4=1;
	float o1=0, o2=0, o3=0, o4=0, d1=1, d2=1, d3=1, d4=1;
	char label1[20]={""}, label2[20]={""}, label3[20]={""}, label4[20]={""};
    const char* data = input.c_str();

    if (input == "in"){
        hetch((char*)"n1",(char*)"d",&n1);
        hetch((char*)"n2",(char*)"d",&n2);
        hetch((char*)"n3",(char*)"d",&n3);
        hetch((char*)"n4",(char*)"d",&n4);
        hetch((char*)"o1",(char*)"f",&o1);
        hetch((char*)"o2",(char*)"f",&o2);
        hetch((char*)"o3",(char*)"f",&o3);
        hetch((char*)"o4",(char*)"f",&o4);
        hetch((char*)"d1",(char*)"f",&d1);
        hetch((char*)"d2",(char*)"f",&d2);
        hetch((char*)"d3",(char*)"f",&d3);
        hetch((char*)"d4",(char*)"f",&d4);
        hetch((char*)"label1",(char*)"s",&label1);
        hetch((char*)"label2",(char*)"s",&label2);
        hetch((char*)"label3",(char*)"s",&label3);
        hetch((char*)"label4",(char*)"s",&label4);
    }
    else{
        auxpar((char*)"n1",(char*)"d",&n1,data);
        auxpar((char*)"n2",(char*)"d",&n2,data);
        auxpar((char*)"n3",(char*)"d",&n3,data);
        auxpar((char*)"n4",(char*)"d",&n4,data);
        auxpar((char*)"o1",(char*)"f",&o1,data);
        auxpar((char*)"o2",(char*)"f",&o2,data);
        auxpar((char*)"o3",(char*)"f",&o3,data);
        auxpar((char*)"o4",(char*)"f",&o4,data);
        auxpar((char*)"d1",(char*)"f",&d1,data);
        auxpar((char*)"d2",(char*)"f",&d2,data);
        auxpar((char*)"d3",(char*)"f",&d3,data);
        auxpar((char*)"d4",(char*)"f",&d4,data);
        auxpar((char*)"label1",(char*)"s",&label1,data);
        auxpar((char*)"label2",(char*)"s",&label2,data);
        auxpar((char*)"label3",(char*)"s",&label3,data);
        auxpar((char*)"label4",(char*)"s",&label4,data);
    }
	
	std::string lab1(&label1[0]);
	std::string lab2(&label2[0]);
    std::string lab3(&label3[0]);
	std::string lab4(&label4[0]);

	axis<T> ax1(n1, (T)o1, (T)d1, lab1);
	axis<T> ax2(n2, (T)o2, (T)d2, lab2);
    axis<T> ax3(n3, (T)o3, (T)d3, lab3);
    axis<T> ax4(n4, (T)o4, (T)d4, lab4);

	long long n123 = n1*n2*n3*n4;

    std::vector<axis<T> > axes = {ax1, ax2, ax3, ax4};
    int ndim = 4;
    if ((n4==1) && (ndim0!=4)){
        ndim -=1;
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

    hypercube<T> hyper(axes);
    std::shared_ptr<vecReg<T> > vec = std::make_shared<vecReg<T> > (hyper);

    if (sizeof(T)==4){
        successCheck(4*n123 == sreed(data,vec->getVals(),4*n123), __FILE__,__LINE__,"Cannot read data\n");
    }
    else if (sizeof(T)==8){
        float* val = new float[n123];
        successCheck(4*n123 == sreed(data,val,4*n123), __FILE__,__LINE__,"Cannot read data\n");
        for (int i=0; i<n123; i++) vec->getVals()[i]=val[i];
        delete val;
    }
	else{
        successCheck(false, __FILE__,__LINE__,"Cannot read into the current template vector.\n");
    }

    return vec;
}