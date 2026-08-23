#!/bin/bash
#SBATCH --job-name=ddpm_A_f10             # Job name
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --nodes=1
#SBATCH --ntasks=8                        # MPI ranks
#SBATCH --cpus-per-task=8                 # OpenMP threads per rank
#SBATCH --hint=nomultithread              # physical cores only (nodes: 128c/256t)
#SBATCH --mem=64G
#SBATCH --time=24:00:00                   # Wall time limit (days-hrs:min:sec)
#SBATCH --output=ddpm_A_f10_%j.log
#SBATCH --error=ddpm_A_f10_%j.err
#SBATCH --partition=cpu

# Stage A production, n_d = 10 (353 600 particles, dt ~ 3.6e-6 s):
#   run 1: unstable single mode, 165 000 steps (sigma*t ~ 3)
#   run 2: stable control (--stable),  50 000 steps
# Requires the measured K from job 01 (permeability summary file).
# Submit with:  sbatch --dependency=afterok:<jobid_01> 02_stageA_finger_nd10.sh

module purge 2>/dev/null || true

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

ND=10
INIT=$(ls porous_finger_${ND}_dx_*_init.gsd 2>/dev/null | head -1)
if [ -z "$INIT" ]; then
    "$PY" create_input_geometry.py "$ND"
    INIT=$(ls porous_finger_${ND}_dx_*_init.gsd | head -1)
fi

# Measured permeability from job 01 (last line, column 2 of the summary).
SUMMARY=$(ls porous_perm_${ND}_dx_*_permeability_summary.dat 2>/dev/null | head -1)
if [ -z "$SUMMARY" ]; then
    echo "ERROR: no permeability summary for n_d=$ND -- run 01_stageA_perm.sh first." >&2
    exit 1
fi
K_M2=$(awk 'END{print $2}' "$SUMMARY")
echo "Using measured K = $K_M2 m^2  (from $SUMMARY)"

echo "── unstable single mode, 165000 steps ──"
"$MPIRUN" -np $SLURM_NTASKS --map-by slot:PE=$SLURM_CPUS_PER_TASK --bind-to core \
    -x PYTHONPATH -x OMP_NUM_THREADS \
    "$PY" "$CASE/run_fingering_single_mode.py" "$INIT" "$K_M2" 165000

echo "── stable control, 50000 steps ──"
"$MPIRUN" -np $SLURM_NTASKS --map-by slot:PE=$SLURM_CPUS_PER_TASK --bind-to core \
    -x PYTHONPATH -x OMP_NUM_THREADS \
    "$PY" "$CASE/run_fingering_single_mode.py" "$INIT" "$K_M2" 50000 --stable
