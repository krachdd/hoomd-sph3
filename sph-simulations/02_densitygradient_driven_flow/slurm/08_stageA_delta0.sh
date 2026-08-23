#!/bin/bash
#SBATCH --job-name=ddpm_A_d0              # Job name
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --nodes=1
#SBATCH --ntasks=8                        # MPI ranks
#SBATCH --cpus-per-task=8                 # OpenMP threads per rank
#SBATCH --hint=nomultithread              # physical cores only (nodes: 128c/256t)
#SBATCH --mem=64G
#SBATCH --time=24:00:00                   # Wall time limit (days-hrs:min:sec)
#SBATCH --output=ddpm_A_d0_%j.log
#SBATCH --error=ddpm_A_d0_%j.err
#SBATCH --partition=cpu

# Stage A linearity check, n_d = 10: rerun the unstable single mode with
# HALF the perturbation amplitude (delta0 = 0.25 R = 50 um instead of the
# default 0.5 R). In the linear regime sigma must be independent of delta0;
# compare the fitted sigma against job 02's default-amplitude run.
# Output files carry a _d00.25 tag, so nothing collides with job 02.
# Submit with:  sbatch --dependency=afterok:<jobid_01> 08_stageA_delta0.sh

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
INIT=$(ls porous_finger_${ND}_dx_*_init.gsd 2>/dev/null | grep -v '_w' | head -1)
if [ -z "$INIT" ]; then
    "$PY" create_input_geometry.py "$ND"
    INIT=$(ls porous_finger_${ND}_dx_*_init.gsd | grep -v '_w' | head -1)
fi

SUMMARY=$(ls porous_perm_${ND}_dx_*_permeability_summary.dat 2>/dev/null | head -1)
if [ -z "$SUMMARY" ]; then
    echo "ERROR: no permeability summary for n_d=$ND -- run 01_stageA_perm.sh first." >&2
    exit 1
fi
K_M2=$(awk 'END{print $2}' "$SUMMARY")
echo "Using measured K = $K_M2 m^2, delta0 = 0.25 R"

"$MPIRUN" -np $SLURM_NTASKS --map-by slot:PE=$SLURM_CPUS_PER_TASK --bind-to core \
    -x PYTHONPATH -x OMP_NUM_THREADS \
    "$PY" "$CASE/run_fingering_single_mode.py" "$INIT" "$K_M2" 165000 --d0frac 0.25
