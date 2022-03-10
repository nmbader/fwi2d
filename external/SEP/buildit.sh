# /bin/bash

rm -rf build
rm -rf local
mkdir build
mkdir local
cd build
cmake -DCMAKE_INSTALL_PREFIX=../local/ \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
  -DCMAKE_C_COMPILER=/usr/bin/gcc ../src/
make && make install
rm -rf build
cd ..
