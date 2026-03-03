#!/bin/bash
#SBATCH --job-name=bcc_1em06             # Job name
#SBATCH --mail-type=BEGIN,END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de    # Where to send mail.  Set this to your email address
#SBATCH --ntasks=27                            # Number of MPI tasks (i.e. processes)
#SBATCH --nodes=1                             # Maximum number of nodes to be allocated
#SBATCH --distribution=cyclic:cyclic          # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --mem=50G
#SBATCH --time=1-24:00:00                       # Wall time limit (days-hrs:min:sec)
#SBATCH --output=bcc_1em06_%j.log         # Path to the standard output and error files relative to the working directory
#SBATCH --error=bcc_1em06_%j.err          # Path to the standard error file
#SBATCH --partition=cpu                       # put the job into the cpu partition

# Ensure that all of the cores are on the same Inifniband network
#SBATCH --contiguous

module load openmpi/5.0.1_gcc-12.2_cuda-12.3
module load gcc/12.2.0

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo "JobID = $SLURM_JOB_ID"

STEPS=${1:-100001}
DAMP=${2:-5000}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INIT="bcc100_init.gsd"
FX="1.0e-6"

# Create geometry if not yet present
if [[ ! -f "${INIT}" ]]; then
    pushd "${SCRIPT_DIR}" > /dev/null
    python3 create_bcc_geometry.py
    popd > /dev/null
else
    echo "  bcc100_init.gsd already exists – skipping geometry creation."
fi

echo "── fx = ${FX} m/s²   (Re_Darcy ≈ 0.015) ──"
/usr/local.nfs/software/openmpi/5.0.1_gcc-12.2_cuda-12.3/bin/mpirun -np 8 ./run_bcc_permeability.py \
    "${INIT}" "${FX}" "${STEPS}" "${DAMP}"
