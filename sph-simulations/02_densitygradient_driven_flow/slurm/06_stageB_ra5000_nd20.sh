#!/bin/bash
#SBATCH --job-name=ddpm_B_ra5k            # Job name
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --nodes=3
#SBATCH --ntasks=48                       # MPI ranks (16/node x 8 = 3 full 128c nodes)
#SBATCH --cpus-per-task=8                 # OpenMP threads per rank
#SBATCH --hint=nomultithread              # physical cores only (nodes: 128c/256t)
#SBATCH --mem=192G                        # per node
#SBATCH --time=2-00:00:00                 # Wall time limit (cpu partition cap)
#SBATCH --output=ddpm_B_ra5k_%j.log
#SBATCH --error=ddpm_B_ra5k_%j.err
#SBATCH --partition=cpu

# Stage B high-Ra member: Ra = 5000 at n_d = 20 (~8.4 M particles,
# dt ~ 1.4e-6 s, 200 000 steps = onset horizon). n_d = 20 is REQUIRED here:
# Pe_grid = Ra*phi*dx/H would be ~12 at n_d = 10 (under-resolved fronts).
#
# SIZING: cpu caps at 2 days and cpu-long is drained, so this heaviest run
# uses THREE nodes (48 ranks x 8 threads = 384 cores, ~175 k particles/rank)
# -- estimated ~1.3 d, margin within the cap. Drop to --nodes=2/--ntasks=32
# if queueing three nodes takes too long (~2 d, tight).
# Cross-node MPI is the standard HOOMD path but untested for this case at
# scale: shake it down first with a short run (STEPS=2000) before
# committing the full horizon. There is no restart mechanism yet -- if the
# job cannot fit the limit, split into legs re-initialised from the last
# dump frame.
# Requires the n_d=20 permeability summary from job 04.
# Submit with:  sbatch --dependency=afterok:<jobid_04> 06_stageB_ra5000_nd20.sh

module purge 2>/dev/null || true

REPO=/data/work/ac126015/imb_test_conda/hoomd-sph3
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
ND=20
RA=5000
INIT=porous_ra_pack_${ND}_seed${SEED}_init.gsd
if [ ! -f "$INIT" ]; then
    "$PY" create_input_geometry.py "$ND" "$SEED"
fi

SUMMARY=porous_ra_perm_${ND}_seed${SEED}_permeability_summary.dat
if [ ! -f "$SUMMARY" ]; then
    echo "ERROR: $SUMMARY missing -- run 04_stageB_perm.sh first." >&2
    exit 1
fi
K_M2=$(awk 'END{print $2}' "$SUMMARY")
PHI=$(awk  'END{print $4}' "$SUMMARY")
D=$("$PY" -c "print(f'{50*9.81*${K_M2}*0.012/(${PHI}*1e-3*${RA}):.4e}')")

echo "Using measured K = $K_M2 m^2, phi = $PHI  ->  Ra = $RA, D = $D m^2/s"

STEPS=200000   # onset horizon; set 2000 for a multi-node shake-down first

"$MPIRUN" -np $SLURM_NTASKS --map-by slot:PE=$SLURM_CPUS_PER_TASK --bind-to core \
    -x PYTHONPATH -x OMP_NUM_THREADS \
    "$PY" "$CASE/run_fingering_ra.py" "$INIT" "$K_M2" "$D" "$STEPS"
