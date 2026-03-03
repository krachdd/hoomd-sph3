#!/bin/bash
#SBATCH --job-name=caprise_c1_t045              # Job name
#SBATCH --mail-type=BEGIN,END,FAIL        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de
#SBATCH --ntasks=20                     # Number of MPI tasks (i.e. processes)
#SBATCH --nodes=1                         # Maximum number of nodes to be allocated
#SBATCH --distribution=cyclic:cyclic      # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --mem=60G
#SBATCH --time=1-24:00:00                   # Wall time limit (days-hrs:min:sec)
#SBATCH --output=caprise_c1_t045_%j.log         # Path to the standard output and error files relative to the working directory
#SBATCH --error=caprise_c1_t045_%j.err          # Path to the standard error file
#SBATCH --partition=cpu                   # put the job into the cpu partition

# Ensure that all of the cores are on the same Inifniband network
#SBATCH --contiguous

module load openmpi/4.1.4_gcc-11.3_cuda-11.7
module load gcc/11.3.0

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo "JobID = $SLURM_JOB_ID"

NL=${1:-10}
STEPS=${2:-100001}
echo "── Case 2  NL=${NL}  steps=${STEPS} ──"
/usr/local.nfs/software/openmpi/4.1.4_gcc-11.3_cuda-11.7/bin/mpirun -np 4 python3 run_capillary_rise.py ${NL} caprise_40_124_40_vs_0.000100_init.gsd 1 ${STEPS}
