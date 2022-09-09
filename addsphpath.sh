#!/bin/bash
export GIT_SRC=$(pwd)
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/hoomd-blue/build
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/dependencies/gsd-2.5.3/build
export PYTHONPATH=$PYTHONPATH:${GIT_SRC}/dependencies/gsd2vtu

