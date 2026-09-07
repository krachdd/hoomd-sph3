#!/bin/bash
#SBATCH --job-name=caprise_gpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --output=caprise_gpu_%j.log
#SBATCH --error=caprise_gpu_%j.err
#
# Capillary rise on one A100 (GPU build of hoomd-sph3, feature/gpu).
#   sbatch caprise_gpu.sh [NL] [CASE] [STEPS]
#     NL    particles across R_cap (default 10 -> caprise_40_124_40 geometry)
#     CASE  row of capillary_rise_params.txt (default 0, theta=30 deg)
#     STEPS simulation steps (default 100001, as the CPU drivers)
# The environment follows buildall_wolfgang.sh: conda env sph3 (own Open MPI,
# launched with mpirun -np 1 -- direct srun aborts in MPI_Init), plus the
# gcc/12.2.0 and cuda/12.3 modules for the runtime libraries.

NL=${1:-10}
CASE=${2:-0}
STEPS=${3:-100001}

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "JobID             = $SLURM_JOB_ID"
echo "GPU               = $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
echo "NL=${NL}  CASE=${CASE}  STEPS=${STEPS}"

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate sph3
module purge
module load gcc/12.2.0 cuda/12.3
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

# SLURM executes a copy of this script from /var/spool/slurmd, so locate the
# case directory via the submission directory (submit from this directory).
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
pushd "${REPO_ROOT}" > /dev/null
source addsphpath.sh
popd > /dev/null

pushd "${SCRIPT_DIR}" > /dev/null
INIT=$(ls -t caprise_*_init.gsd 2>/dev/null | head -1)
if [[ -z "${INIT}" ]]; then
    echo "── Creating capillary geometry (NL=${NL}) ──"
    mpirun -np 1 python3 create_capillary_geometry.py ${NL}
    INIT=$(ls -t caprise_*_init.gsd | head -1)
fi
echo "── Case ${CASE}  NL=${NL}  steps=${STEPS}  init=${INIT} ──"
mpirun -np 1 --bind-to none python3 run_capillary_rise_gpu.py ${NL} "${INIT}" ${CASE} ${STEPS}
popd > /dev/null
