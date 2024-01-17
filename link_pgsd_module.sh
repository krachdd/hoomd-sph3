#!/bin/bash
export GIT_SRC=$(pwd)
ln -s ${GIT_SRC}/dependencies/pgsd-sph/pgsd-3.2.0/pgsd/pgsd.c ${GIT_SRC}/hoomd-blue/hoomd/extern/pgsd.c
ln -s ${GIT_SRC}/dependencies/pgsd-sph/pgsd-3.2.0/pgsd/pgsd.h ${GIT_SRC}/hoomd-blue/hoomd/extern/pgsd.h