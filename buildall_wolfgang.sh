#!/bin/bash
# GPU build of hoomd-sph3 for the wolfgang cluster (wolfgang-gpu0[1-4]: 4x A100, sm_80).
#
# Environment philosophy follows the existing CPU build on wolfgang: the conda
# env `sph3` is self-contained (Python, Open MPI 4.1.6 via mpicc/mpicxx), all
# job scripts do `module purge`.  Only CUDA (and a matching host gcc) come from
# modules.  The conda Open MPI is NOT CUDA-aware — fine for single-GPU and for
# multi-GPU with host-staged MPI buffers (HOOMD's default); switch to one of
# the `openmpi/*_cuda-*` modules only if you rebuild mpi4py against it too.
#
# Usage (login node, from the repo root):
#     conda activate sph3
#     ./buildall_wolfgang.sh            # gsd dependency + hoomd-blue
#     ./buildall_wolfgang.sh --hoomd    # hoomd-blue only (gsd already built)
#
# The build is a plain CPU compile job for nvcc; no GPU is needed on the login
# node (the "No NVidia GPU found" warning from `module load cuda` is harmless).
set -euo pipefail

export GIT_SRC=$(pwd)

# ── Environment ────────────────────────────────────────────────────────────
module purge
# CUDA 12.3 supports host gcc <= 12.2; pair them explicitly so nvcc does not
# pick up an unsupported system compiler.
module load gcc/12.2.0
module load cuda/12.3
# HOOMD sets CUDA_STANDARD 17, which needs CMake >= 3.18 (the login node's
# /usr/bin/cmake is 3.16).
module load cmake/3.26.3

export CC=$(which gcc)
export CXX=$(which g++)
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=$CXX
# The conda Open MPI wrappers (mpicc/mpicxx) are configured for the conda
# compiler x86_64-conda-linux-gnu-cc, which is not installed in the sph3 env.
# Open MPI honours OMPI_CC/OMPI_CXX, so route the wrappers to the module gcc.
export OMPI_CC=$CC
export OMPI_CXX=$CXX

# CUDA architectures to compile for: A100 = sm_80.
# Add ";90" for H100 or ";75" for Turing cards if the binary must run there too.
CUDA_ARCH_LIST=${CUDA_ARCH_LIST:-80}

NJOBS=${NJOBS:-16}

echo "nvcc:  $(nvcc --version | tail -1)"
echo "gcc:   $($CXX --version | head -1)"
echo "mpicc: $(which mpicc)  -> $(mpicc --version | head -1)"
mpicc --version >/dev/null 2>&1 || { echo "mpicc wrapper is not working (check OMPI_CC / conda env sph3 active)"; exit 1; }
echo "cmake: $(cmake --version | head -1)"
echo "CUDA_ARCH_LIST=$CUDA_ARCH_LIST"

# ── Dependency (gsd) ───────────────────────────────────────────────────────
# pgsd is NOT built here: hoomd-blue compiles pgsd.c directly via the symlinks
# in hoomd/extern (created by link_pgsd_module.sh), so the standalone pgsd
# build is unnecessary and its cmake step trips over the conda MPI wrappers.
if [[ "${1:-}" != "--hoomd" ]]; then
    # the symlink step reports 'File exists' on a re-run; that is harmless
    ./link_pgsd_module.sh || true

    cd "$GIT_SRC/dependencies/gsd-sph/gsd/"
    rm -rf build && mkdir build && cd build
    cmake ..
    make -j"$NJOBS"
fi

# ── hoomd-blue with the sph component, GPU enabled ────────────────────────
cd "$GIT_SRC/hoomd-blue/"
rm -rf build && mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_GPU=ON \
    -DENABLE_MPI=ON \
    -DCUDA_ARCH_LIST="$CUDA_ARCH_LIST" \
    -DMPI_C_COMPILER=$(which mpicc) \
    -DMPI_CXX_COMPILER=$(which mpicxx)
make -j"$NJOBS"

echo
echo "Build finished. Sanity check on a GPU node, e.g.:"
echo "  srun -p gpu --gres=gpu:A100:1 -c 8 --pty bash"
echo "  cd $GIT_SRC/hoomd-blue/build && python3 -c 'import hoomd,hoomd.sph; print(hoomd.version.gpu_enabled); hoomd.device.GPU()'"
