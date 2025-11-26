#!/bin/bash
#SBATCH --job-name=ddlomg_      # Job name
#SBATCH --mail-type=BEGIN,END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=david.krach@mib.uni-stuttgart.de    # Where to send mail.  Set this to your email address
#SBATCH --ntasks=80                 # Number of MPI tasks (i.e. processes)
#SBATCH --nodes=1                    # Maximum number of nodes to be allocated
#SBATCH --mem=100G          
#SBATCH --time=01-00:00:00              # Wall time limit (days-hrs:min:sec)
#SBATCH --output=ddlomg_%j.log     # Path to the standard output and error files relative to the working directory
#SBATCH --partition=cpu             # put the job into the cpu partition

# OUTPUT und FEHLER Dateien. %j wird durch job id ersetzt.
# SBATCH -o mpi_test_job.out # File to which STDOUT will be written
#SBATCH -e ddlomg_%j.err # File to which STDERR will be written

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

/usr/local.nfs/software/openmpi/4.1.4_gcc-11.3_cuda-11.7/bin/mpirun -np $SLURM_NTASKS  ./run_dd_tv.py --resolution=100 --initgsd="dd_208_50_11_vs_5e-05_re_1_init.gsd"