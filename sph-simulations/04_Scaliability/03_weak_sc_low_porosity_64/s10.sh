#!/bin/bash
#SBATCH --job-name=lp_s10      # Job name
#SBATCH --mail-type=BEGIN,END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de    # Where to send mail.  Set this to your email address
#SBATCH --ntasks=640                  # Number of MPI tasks (i.e. processes)
#SBATCH --nodes=10                    # Maximum number of nodes to be allocated
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --mem=MaxMemPerNode          
#SBATCH --time=2:00:00              # Wall time limit (days-hrs:min:sec)
#SBATCH --output=lp_s10%j.log     # Path to the standard output and error files relative to the working directory
#SBATCH --partition=cpu              # put the job into the cpu partition



# OUTPUT und FEHLER Dateien. %j wird durch job id ersetzt.
# SBATCH -o lp_run_1.out # File to which STDOUT will be written
#SBATCH -e lp_s10%j.err # File to which STDERR will be written

module load openmpi/4.1.4_gcc-11.3_cuda-11.7
module load gcc/11.3.0 

# Exports
PYTHON=/home/ac126015/software/miniconda3/envs/sph3/bin/python3.11
PY_LOCAL_LIB=/home/ac126015/software/miniconda3/envs/sph3/lib
PY_LOCAL_BIN=/home/ac126015/software/miniconda3/envs/sph3/bin
PY_LOCAL_INC=/home/ac126015/software/miniconda3/envs/sph3/include
PY_LOCAL_SIT=/home/ac126015/software/miniconda3/envs/sph3/lib/python3.11/site-packages

export PYTHONPATH=$PY_LOCAL_SIT:$PYTHONPATH
export PATH=$PY_LOCAL_BIN:$PATH

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo "JobID = $SLURM_JOB_ID"

/usr/local.nfs/software/openmpi/4.1.4_gcc-11.3_cuda-11.7/bin/mpirun -np $SLURM_NTASKS ./lp_sc.py input_lp10.txt 10