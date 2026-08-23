#!/bin/bash
#SBATCH --job-name=ddpm_A_perm            # Job name
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --nodes=1
#SBATCH --ntasks=4                        # MPI ranks
#SBATCH --cpus-per-task=8                 # OpenMP threads per rank
#SBATCH --hint=nomultithread              # physical cores only (nodes: 128c/256t)
#SBATCH --mem=32G
#SBATCH --time=04:00:00                   # Wall time limit (days-hrs:min:sec)
#SBATCH --output=ddpm_A_perm_%j.log
#SBATCH --error=ddpm_A_perm_%j.err
#SBATCH --partition=cpu

# Stage A permeability pre-runs (cylinder-array unit cells, n_d = 10 and 20).
# Writes porous_perm_*_permeability_summary.dat, which the fingering jobs
# (02/03) read to obtain the measured K -- submit this job FIRST.

# Modules stay purged: the sph3 conda env is self-contained (OpenMPI,
# compilers, python); inherited module LD_LIBRARY_PATHs could shadow it.
module purge 2>/dev/null || true

# wolfgang layout (adjust if the clone moves)
REPO=/data/work/ac126015/dd_pm/hoomd-sph3
CASE="$REPO/sph-simulations/02_densitygradient_driven_flow/04_porous_fingering_single_mode"
export PYTHONPATH="$REPO/hoomd-blue/build:$REPO/dependencies/gsd-sph/gsd/build:$REPO/helper_modules"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

MPIRUN=/home/ac126015/software/miniconda3/envs/sph3/bin/mpirun
PY=/home/ac126015/software/miniconda3/envs/sph3/bin/python3

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo "JobID = $SLURM_JOB_ID"

cd "$CASE"

# Launch with the SAME env's mpirun (PRRTE brings its own PMIx; no SLURM
# plugin or system module needed). -x forwards PYTHONPATH/OMP_NUM_THREADS
# to every rank. A mismatched launcher fractures the MPI world into
# singletons (the "GSD: Not a GSD file" race).

for ND in 10 20; do
    INIT=$(ls porous_perm_${ND}_dx_*_init.gsd 2>/dev/null | head -1)
    if [ -z "$INIT" ]; then
        "$PY" create_input_geometry.py "$ND" --perm
        INIT=$(ls porous_perm_${ND}_dx_*_init.gsd | head -1)
    fi
    echo "── Stage A permeability, n_d = $ND, init = $INIT ──"
    # fx = 1e-4 m/s^2 (Darcy regime), 20001 steps, 2000-step ramp
    "$MPIRUN" -np $SLURM_NTASKS --map-by slot:PE=$SLURM_CPUS_PER_TASK --bind-to core \
        -x PYTHONPATH -x OMP_NUM_THREADS \
        "$PY" "$CASE/run_permeability.py" "$INIT" 1e-4 20001 2000
done
