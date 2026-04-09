#!/bin/bash
#SBATCH --job-name=MT_POD_DD
#SBATCH --partition=mpi64c
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=4
#SBATCH --nodes=8
#SBATCH --exclusive
#SBATCH --mem-per-cpu=3G
#SBATCH --time=1-00:00:00
#SBATCH --output=jobs/mt_test_%j.out
#SBATCH --error=jobs/mt_test_%j.err

echo "Job started at: $(date)"
echo "Running on nodes: $SLURM_NODELIST"
cd $SLURM_SUBMIT_DIR

########## ENVIRONMENT ##########
export PATH=/clonetroop/dasilva/libs/bin:$PATH
export LD_LIBRARY_PATH=/clonetroop/dasilva/libs/lib:$LD_LIBRARY_PATH

########## THREADING CONFIGURATION ##########
# Disable OpenMP in your code
export OMP_NUM_THREADS=1

# Use 4 threads for MKL per MPI rank
export MKL_NUM_THREADS=4
export MKL_DYNAMIC=false
export KMP_AFFINITY=compact

########## RUN ##########
mpirun -np 40 ./max3d-par

echo "Job finished at: $(date)"








