#!/bin/bash
export GIT_SRC=$(pwd)
./link_pgsd_module.sh
cd dependencies/pgsd-sph/pgsd-3.2.0/
rm -rf build
mkdir build 
cd build
CC=/usr/bin/mpicc CXX=/usr/bin/mpicxx cmake .. 
make 
cd $GIT_SRC
cd dependencies/gsd-sph/gsd-3.4.1/
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
