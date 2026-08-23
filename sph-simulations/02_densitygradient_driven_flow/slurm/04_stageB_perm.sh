#!/bin/bash
#SBATCH --job-name=ddpm_B_perm            # Job name
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --nodes=1
#SBATCH --ntasks=4                        # MPI ranks
#SBATCH --cpus-per-task=8                 # OpenMP threads per rank
#SBATCH --hint=nomultithread              # physical cores only (nodes: 128c/256t)
#SBATCH --mem=32G
#SBATCH --time=06:00:00                   # Wall time limit (days-hrs:min:sec)
#SBATCH --output=ddpm_B_perm_%j.log
#SBATCH --error=ddpm_B_perm_%j.err
#SBATCH --partition=cpu

# Stage B permeability pre-runs (disordered-pack periodic cells, seed 42,
# n_d = 10 and 20). Writes porous_ra_perm_*_permeability_summary.dat, which
# the Ra-sweep jobs (05/06) read -- submit this job FIRST.

module purge 2>/dev/null || true

REPO=/data/work/ac126015/dd_pm/hoomd-sph3
CASE="$REPO/sph-simulations/02_densitygradient_driven_flow/05_porous_fingering_rayleigh"
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

SEED=42
for ND in 10 20; do
    INIT=porous_ra_perm_${ND}_seed${SEED}_init.gsd
    if [ ! -f "$INIT" ]; then
        "$PY" create_input_geometry.py "$ND" "$SEED" --perm
    fi
    echo "── Stage B permeability, n_d = $ND, init = $INIT ──"
    "$MPIRUN" -np $SLURM_NTASKS --map-by slot:PE=$SLURM_CPUS_PER_TASK --bind-to core \
        -x PYTHONPATH -x OMP_NUM_THREADS \
        "$PY" "$CASE/run_permeability.py" "$INIT" 1e-4 20001 2000
done
