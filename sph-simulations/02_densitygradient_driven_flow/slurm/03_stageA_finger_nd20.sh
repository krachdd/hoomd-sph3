#!/bin/bash
#SBATCH --job-name=ddpm_A_f20             # Job name
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --nodes=2
#SBATCH --ntasks=32                       # MPI ranks (16/node x 8 = 2 full 128c nodes)
#SBATCH --cpus-per-task=8                 # OpenMP threads per rank
#SBATCH --hint=nomultithread              # physical cores only (nodes: 128c/256t)
#SBATCH --mem=128G                        # per node
#SBATCH --time=2-00:00:00                 # Wall time limit (cpu partition cap)
#SBATCH --output=ddpm_A_f20_%j.log
#SBATCH --error=ddpm_A_f20_%j.err
#SBATCH --partition=cpu

# Stage A convergence member, n_d = 20 (~2.83 M particles, dt halves vs n_d=10):
#   unstable single mode, 330 000 steps (same physical horizon sigma*t ~ 3).
# This is the resolution-convergence run for the paper (Table: sigma vs n_d).
# SIZING: cpu-long is drained, so this runs on cpu (2-day cap) with TWO
# nodes (32 x 8 = 256 cores, ~88 k particles/rank) -- estimated ~1.5 d,
# margin within the cap. Emergency fallback if even that is too slow:
# STEPS=220000 (sigma*t ~ 2, still enough for the growth-rate fit).
# When cpu-long returns: --nodes=1 --ntasks=16 --time=3-00:00:00
# --partition=cpu-long works too.
# Requires the n_d=20 permeability summary from job 01.
# Submit with:  sbatch --dependency=afterok:<jobid_01> 03_stageA_finger_nd20.sh

module purge 2>/dev/null || true

REPO=/data/work/ac126015/imb_test_conda/hoomd-sph3
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

ND=20
INIT=$(ls porous_finger_${ND}_dx_*_init.gsd 2>/dev/null | head -1)
if [ -z "$INIT" ]; then
    "$PY" create_input_geometry.py "$ND"
    INIT=$(ls porous_finger_${ND}_dx_*_init.gsd | head -1)
fi

SUMMARY=$(ls porous_perm_${ND}_dx_*_permeability_summary.dat 2>/dev/null | head -1)
if [ -z "$SUMMARY" ]; then
    echo "ERROR: no permeability summary for n_d=$ND -- run 01_stageA_perm.sh first." >&2
    exit 1
fi
K_M2=$(awk 'END{print $2}' "$SUMMARY")
echo "Using measured K = $K_M2 m^2  (from $SUMMARY)"

STEPS=330000   # sigma*t ~ 3;  cpu-partition fallback: 220000 (sigma*t ~ 2)

"$MPIRUN" -np $SLURM_NTASKS --map-by slot:PE=$SLURM_CPUS_PER_TASK --bind-to core \
    -x PYTHONPATH -x OMP_NUM_THREADS \
    "$PY" "$CASE/run_fingering_single_mode.py" "$INIT" "$K_M2" "$STEPS"
