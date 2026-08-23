#!/bin/bash
#SBATCH --job-name=ddpm_A_sk              # Job name
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --nodes=1
#SBATCH --ntasks=16                       # MPI ranks (16 x 8 = full 128c node)
#SBATCH --cpus-per-task=8                 # OpenMP threads per rank
#SBATCH --hint=nomultithread              # physical cores only (nodes: 128c/256t)
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00                 # Wall time limit (cpu partition cap)
#SBATCH --array=0-1                       # 0: Lx=3a (w3), 1: Lx=7a (w7)
#SBATCH --output=ddpm_A_sk_%A_%a.log
#SBATCH --error=ddpm_A_sk_%A_%a.err
#SBATCH --partition=cpu

# Stage A sigma(k) dispersion study, n_d = 10: two additional domain widths
# around the default Lx = 5a (job 02), giving three wavenumbers total:
#   w3: Lx = 2.4 mm (~212 k particles, sigma ~ 7.4 1/s, 120 000 steps)
#   w7: Lx = 5.6 mm (~495 k particles, sigma ~ 4.0 1/s, 210 000 steps)
# Together with job 02's Lx = 4.0 mm these trace the dispersion relation
#   sigma(k) = U_b*k/(2*phi) - D*k^2
# and probe where Darcy-scale theory breaks as d/lambda grows.
# K is the same lattice as job 02 -> reuses the n_d=10 permeability summary.
# Submit with:  sbatch --dependency=afterok:<jobid_01> 07_stageA_sigma_k.sh

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
echo "JobID = $SLURM_JOB_ID  ArrayTask = $SLURM_ARRAY_TASK_ID"

cd "$CASE"

ND=10
case $SLURM_ARRAY_TASK_ID in
    0) W=3; STEPS=120000 ;;
    1) W=7; STEPS=210000 ;;
esac

INIT=$(ls porous_finger_w${W}_${ND}_dx_*_init.gsd 2>/dev/null | head -1)
if [ -z "$INIT" ]; then
    "$PY" create_input_geometry.py "$ND" --cells-x "$W"
    INIT=$(ls porous_finger_w${W}_${ND}_dx_*_init.gsd | head -1)
fi

SUMMARY=$(ls porous_perm_${ND}_dx_*_permeability_summary.dat 2>/dev/null | head -1)
if [ -z "$SUMMARY" ]; then
    echo "ERROR: no permeability summary for n_d=$ND -- run 01_stageA_perm.sh first." >&2
    exit 1
fi
K_M2=$(awk 'END{print $2}' "$SUMMARY")
echo "sigma(k) member w$W: init = $INIT, K = $K_M2 m^2, steps = $STEPS"

"$MPIRUN" -np $SLURM_NTASKS --map-by slot:PE=$SLURM_CPUS_PER_TASK --bind-to core \
    -x PYTHONPATH -x OMP_NUM_THREADS \
    "$PY" "$CASE/run_fingering_single_mode.py" "$INIT" "$K_M2" "$STEPS"
