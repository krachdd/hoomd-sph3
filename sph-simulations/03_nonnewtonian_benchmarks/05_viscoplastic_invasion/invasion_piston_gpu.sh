#!/bin/bash
#SBATCH --job-name=inv_piston_gpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=inv_piston_gpu_%j.log
#SBATCH --error=inv_piston_gpu_%j.err
#
# WP2 piston-driven invasion on one A100 (single rank, SPH_DEVICE=gpu).
#   sbatch invasion_piston_gpu.sh [TAU_Y] [U_P] [SIGMA] [RES] [STEPS] [EXTRA...]
#     TAU_Y  defender yield stress [Pa]   (default 0 = Newtonian)
#     U_P    piston velocity [m/s]         (default 0.04)
#     SIGMA  surface tension [N/m]         (default 0.004)
#     RES    throat resolution [particles] (default 12, WP2 sweep)
#     STEPS  hard step budget              (default 92000)
#     EXTRA  passed through to run_piston_invasion.py (e.g. --ramp 2000 --resA_mm 12)
# The domain piston_domain_res${RES}_init.gsd is created if missing.

TAU_Y=${1:-0}; U_P=${2:-0.04}; SIGMA=${3:-0.004}; RES=${4:-12}; STEPS=${5:-92000}
shift 5 2>/dev/null; EXTRA="$@"

echo "Date=$(date)  Host=$(hostname -s)  JobID=$SLURM_JOB_ID"
echo "GPU=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
echo "tau_y=${TAU_Y} U_p=${U_P} sigma=${SIGMA} res=${RES} steps=${STEPS} extra='${EXTRA}'"

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate sph3
module purge
module load gcc/12.2.0 cuda/12.3
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export SPH_DEVICE=gpu

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
pushd "${REPO_ROOT}" > /dev/null; source addsphpath.sh; popd > /dev/null

pushd "${SCRIPT_DIR}" > /dev/null
INIT=piston_domain_res${RES}_init.gsd
if [[ ! -f "${INIT}" ]]; then
    echo "── Creating piston domain (res=${RES}, U_p=${U_P}) ──"
    mpirun -np 1 python3 create_piston_domain.py --res ${RES} --U_p ${U_P} --out "${INIT}"
fi
mpirun -np 1 --bind-to none python3 run_piston_invasion.py "${INIT}" ${TAU_Y} \
    --U_p ${U_P} --sigma ${SIGMA} --res ${RES} --steps ${STEPS} ${EXTRA}
popd > /dev/null
