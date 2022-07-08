FROM ubuntu:18.04 as builder
MAINTAINER nmbader@sep.stanford.edu
RUN apt-get -y update
RUN apt-get -y install build-essential
RUN apt-get -y install libelf-dev libffi-dev
RUN apt-get -y install pkg-config
RUN apt-get -y install wget git gcc g++ gfortran make cmake vim lsof
RUN apt-get -y install flex libxaw7-dev
RUN apt-get -y install libfftw3-3 libfftw3-dev libssl-dev

RUN apt-get -y  install python3-pip
RUN python3 -m pip install --no-cache-dir --upgrade pip

RUN python3 -m pip install --no-cache-dir numpy &&\
    python3 -m pip install --no-cache-dir jupyter &&\
    python3 -m pip install --no-cache-dir scipy &&\
    python3 -m pip install --no-cache-dir matplotlib

RUN mkdir -p /opt/ispc/bin
RUN mkdir -p /home
RUN mkdir -p /opt/fwi2d
RUN mkdir -p /home/fwi2d
WORKDIR /home

RUN wget https://github.com/ispc/ispc/releases/download/v1.17.0/ispc-v1.17.0-linux.tar.gz  &&\
    tar -xvf ispc-v1.17.0-linux.tar.gz &&\
    cp ispc-v1.17.0-linux/bin/ispc /opt/ispc/bin/ &&\
    rm -f ispc-v1.17.0-linux.tar.gz  &&\
    rm -rf ispc-v1.17.0-linux

ADD src /home/fwi2d/src
ADD external /home/fwi2d/external
ADD examples /home/fwi2d/examples
ADD CMakeLists.txt /home/fwi2d
ADD LICENSE /home/fwi2d
ADD README.md /home/fwi2d

RUN cd /home/fwi2d &&\
    mkdir -p build &&\
    cd external/SEP &&\
    bash ./buildit.sh &&\
    cd ../../build  &&\
    cmake -DCMAKE_INSTALL_PREFIX=/opt/fwi2d/ -DISPC=/opt/ispc/bin/ispc ../  &&\
    make -j12  &&\
    make install &&\
    cd ../ &&\
    rm -rf build

RUN apt-get -y clean

ENV HOME=/home 
ENV PATH="/opt/fwi2d/bin:${PATH}"
ENV DATAPATH="/tmp/"
RUN echo 'alias python=python3' >> ~/.bashrc