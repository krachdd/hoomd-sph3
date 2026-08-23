#!/bin/bash
#SBATCH --job-name=ddpm_B_ra              # Job name
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --nodes=1
#SBATCH --ntasks=16                       # MPI ranks (16 x 8 = full 128c node)
#SBATCH --cpus-per-task=8                 # OpenMP threads per rank
#SBATCH --hint=nomultithread              # physical cores only (nodes: 128c/256t)
#SBATCH --mem=64G
#SBATCH --time=48:00:00                   # Wall time limit (days-hrs:min:sec)
#SBATCH --array=0-2                       # 0: Ra=200, 1: Ra=1000, 2: stable ctrl
#SBATCH --output=ddpm_B_ra_%A_%a.log
#SBATCH --error=ddpm_B_ra_%A_%a.err
#SBATCH --partition=cpu

# Stage B production Ra sweep, n_d = 10, full pack (~1.05 M particles,
# dt ~ 2.9e-6 s). Array members:
#   0: Ra = 200   unstable      100 000 steps (onset horizon, ~0.13 t_c)
#   1: Ra = 1000  unstable      100 000 steps
#   2: Ra = 1000  stable ctrl    60 000 steps
# For the mixing-width scaling W(t), extend to 237 000 steps (~0.3 t_c)
# and raise --time accordingly. Ra = 5000 needs n_d = 20 (job 06):
# Pe_grid = Ra*phi*dx/H caps n_d=10 at Ra ~ 1600.
# Requires the n_d=10 permeability summary from job 04.
# Submit with:  sbatch --dependency=afterok:<jobid_04> 05_stageB_ra_sweep.sh

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

SEED=42
ND=10
INIT=porous_ra_pack_${ND}_seed${SEED}_init.gsd
if [ ! -f "$INIT" ]; then
    # geometry creation is cheap but must not race between array members:
    # member 0 creates, the others wait for the file
    if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
        "$PY" create_input_geometry.py "$ND" "$SEED"
    else
        until [ -f "$INIT" ]; do sleep 30; done
        sleep 60   # let the writer finish
    fi
fi

SUMMARY=porous_ra_perm_${ND}_seed${SEED}_permeability_summary.dat
if [ ! -f "$SUMMARY" ]; then
    echo "ERROR: $SUMMARY missing -- run 04_stageB_perm.sh first." >&2
    exit 1
fi
K_M2=$(awk 'END{print $2}' "$SUMMARY")
PHI=$(awk  'END{print $4}' "$SUMMARY")

# D from the Rayleigh-Darcy number: D = drho*g*K*H/(phi*mu*Ra)
# drho = 50 kg/m^3, g = 9.81, H = 0.012 m, mu = 1e-3 Pa s (run-script defaults)
ra_to_D () { "$PY" -c "print(f'{50*9.81*${K_M2}*0.012/(${PHI}*1e-3*$1):.4e}')"; }

case $SLURM_ARRAY_TASK_ID in
    0) RA=200;  D=$(ra_to_D $RA); STEPS=100000; EXTRA="" ;;
    1) RA=1000; D=$(ra_to_D $RA); STEPS=100000; EXTRA="" ;;
    2) RA=1000; D=$(ra_to_D $RA); STEPS=60000;  EXTRA="--stable" ;;
esac

echo "Using measured K = $K_M2 m^2, phi = $PHI  ->  Ra = $RA, D = $D m^2/s, steps = $STEPS $EXTRA"

"$MPIRUN" -np $SLURM_NTASKS --map-by slot:PE=$SLURM_CPUS_PER_TASK --bind-to core \
    -x PYTHONPATH -x OMP_NUM_THREADS \
    "$PY" "$CASE/run_fingering_ra.py" "$INIT" "$K_M2" "$D" "$STEPS" $EXTRA
