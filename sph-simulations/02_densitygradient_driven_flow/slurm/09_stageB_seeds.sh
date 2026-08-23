#!/bin/bash
#SBATCH --job-name=ddpm_B_seed            # Job name
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --nodes=1
#SBATCH --ntasks=16                       # MPI ranks (16 x 8 = full 128c node)
#SBATCH --cpus-per-task=8                 # OpenMP threads per rank
#SBATCH --hint=nomultithread              # physical cores only (nodes: 128c/256t)
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00                 # Wall time limit (cpu partition cap)
#SBATCH --array=0-1                       # 0: seed 43, 1: seed 44
#SBATCH --output=ddpm_B_seed_%A_%a.log
#SBATCH --error=ddpm_B_seed_%A_%a.err
#SBATCH --partition=cpu

# Stage B disorder-statistics replicas, n_d = 10: repeat the Ra = 1000
# member of job 05 with two additional pack realisations (seeds 43, 44).
# Onset time and dominant wavelength must be statistically consistent
# across seeds -- the robustness check referees expect for a disordered
# medium. Each member is SELF-CONTAINED for its seed: it builds the
# periodic permeability cell, measures K (each realisation has its own K!),
# builds the pack, and runs the fingering case. No dependency on other jobs.
# Budget: perm ~1-2 h + fingering 100 k steps ~1.5 d -> fits the 2-day cap;
# fallback if tight: STEPS=80000.
# Submit with:  sbatch 09_stageB_seeds.sh

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
echo "JobID = $SLURM_JOB_ID  ArrayTask = $SLURM_ARRAY_TASK_ID"

cd "$CASE"

ND=10
RA=1000
STEPS=100000
SEED=$((43 + SLURM_ARRAY_TASK_ID))

# 1) permeability of THIS realisation (each seed has its own K)
PERM_INIT=porous_ra_perm_${ND}_seed${SEED}_init.gsd
SUMMARY=porous_ra_perm_${ND}_seed${SEED}_permeability_summary.dat
if [ ! -f "$SUMMARY" ]; then
    [ -f "$PERM_INIT" ] || "$PY" create_input_geometry.py "$ND" "$SEED" --perm
    "$MPIRUN" -np $SLURM_NTASKS --map-by slot:PE=$SLURM_CPUS_PER_TASK --bind-to core \
        -x PYTHONPATH -x OMP_NUM_THREADS \
        "$PY" "$CASE/run_permeability.py" "$PERM_INIT" 1e-4 20001 2000
fi
K_M2=$(awk 'END{print $2}' "$SUMMARY")
PHI=$(awk  'END{print $4}' "$SUMMARY")

# 2) fingering at Ra = 1000 on this seed's pack
INIT=porous_ra_pack_${ND}_seed${SEED}_init.gsd
[ -f "$INIT" ] || "$PY" create_input_geometry.py "$ND" "$SEED"

D=$("$PY" -c "print(f'{50*9.81*${K_M2}*0.012/(${PHI}*1e-3*${RA}):.4e}')")
echo "seed $SEED: K = $K_M2 m^2, phi = $PHI  ->  Ra = $RA, D = $D m^2/s, steps = $STEPS"

"$MPIRUN" -np $SLURM_NTASKS --map-by slot:PE=$SLURM_CPUS_PER_TASK --bind-to core \
    -x PYTHONPATH -x OMP_NUM_THREADS \
    "$PY" "$CASE/run_fingering_ra.py" "$INIT" "$K_M2" "$D" "$STEPS"
