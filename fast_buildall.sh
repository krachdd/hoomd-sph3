#!/bin/bash
export GIT_SRC=$(pwd)
export CUDACXX=/usr/local/cuda-13.2/bin/nvcc
export PATH=/usr/local/cuda-13.2/bin:$PATH
./link_pgsd_module.sh
cd dependencies/pgsd-sph/pgsd/
rm -rf build
mkdir build 
cd build
CC=/usr/bin/mpicc CXX=/usr/bin/mpicxx cmake .. 
make 
cd $GIT_SRC
cd dependencies/gsd-sph/gsd/
rm -rf build
mkdir build
cd build 
cmake ..
make 

cd $GIT_SRC
cd hoomd-blue/
rm -rf build
cmake -B build -S . -GNinja -DENABLE_GPU=ON -DCUDA_ARCH_LIST=75
cd build 
ninja
