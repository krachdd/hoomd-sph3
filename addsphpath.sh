#!/bin/bash
export GIT_SRC=$(pwd)
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/hoomd-blue/build
# export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/dependencies/gsd-2.5.3/build
# export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/dependencies/gsd-2.7.0/build
# export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/dependencies/gsd-2.8.0/build
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/dependencies/gsd-3.2.0/build
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/dependencies/pgsd-3.2.0/build
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/helper_modules/gsd2vtu
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/helper_modules/delete_solid_sphparticles
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/helper_modules/read_input_from_txt
