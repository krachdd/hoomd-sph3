#!/bin/bash
#SBATCH --job-name=inv_bodyforce_gpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=inv_bodyforce_gpu_%j.log
#SBATCH --error=inv_bodyforce_gpu_%j.err
#
# WP2 body-force-driven invasion on one A100 (single rank, SPH_DEVICE=gpu).
#   sbatch invasion_bodyforce_gpu.sh [TAU_Y] [G] [SIGMA] [RES] [STEPS] [EXTRA...]
#     TAU_Y  defender yield stress [Pa]       (default 150 -> Y=12)
#     G      body force along +x [m/s^2]      (default 168.75 = 0.9 g_crit at Y=12)
#     SIGMA  surface tension [N/m]             (default 0.01)
#     RES    throat resolution [particles]     (default 12, WP2 sweep)
#     STEPS  step budget                       (default 50000)
#     EXTRA  passed through to run_bodyforce_invasion.py (e.g. --vcap 1.0 --damp 1000)
# The piston-free domain bodyforce_domain_res${RES}_init.gsd is created if missing.

TAU_Y=${1:-150}; G=${2:-168.75}; SIGMA=${3:-0.01}; RES=${4:-12}; STEPS=${5:-50000}
shift 5 2>/dev/null; EXTRA="$@"

echo "Date=$(date)  Host=$(hostname -s)  JobID=$SLURM_JOB_ID"
echo "GPU=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
echo "tau_y=${TAU_Y} g=${G} sigma=${SIGMA} res=${RES} steps=${STEPS} extra='${EXTRA}'"

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
INIT=bodyforce_domain_res${RES}_init.gsd
if [[ ! -f "${INIT}" ]]; then
    echo "── Creating piston-free domain (res=${RES}) ──"
    mpirun -np 1 python3 create_piston_domain.py --res ${RES} --no-piston --out "${INIT}"
fi
mpirun -np 1 --bind-to none python3 run_bodyforce_invasion.py "${INIT}" ${TAU_Y} \
    --g ${G} --sigma ${SIGMA} --res ${RES} --steps ${STEPS} ${EXTRA}
popd > /dev/null
