#!/bin/bash
export GIT_SRC=$(pwd)
./link_pgsd_module.sh
cd dependencies/pgsd-sph/pgsd-3.2.0/
rm -rf build
mkdir build 
cd build
CC=/usr/local.nfs/software/openmpi/4.1.4_gcc-11.3_cuda-11.7/bin/mpicc CXX=/usr/local.nfs/software/openmpi/4.1.4_gcc-11.3_cuda-11.7/bin/mpicxx cmake .. 
make 
cd $GIT_SRC
cd dependencies/gsd-sph/gsd-3.2.0/
rm -rf build
mkdir build
cd build 
cmake ..
make 

cd $GIT_SRC
cd hoomd-blue/
rm -rf build
mkdir build
cd build 
cmake ..
make -j4 
