#!/bin/bash
export GIT_SRC=$(pwd)
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/hoomd-blue/build
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/dependencies/gsd-sph/gsd-3.2.0/build
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/dependencies/pgsd-sph/pgsd-3.2.0/build
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/helper_modules/gsd2vtu
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/helper_modules/
