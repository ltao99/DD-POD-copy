#!/bin/bash
#PBS -N max3d-par
#PBS -l walltime=720:00:00 
#PBS -l nodes=n0017:ppn=50
#PBS -o stdout_serial
#PBS -e stderr_serial
#PBS -V
#PBS -m abe

cd $PBS_O_WORKDIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/intel/mkl/lib/intel64/

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/compilers_and_libraries/linux/lib/intel64

cd $PBS_O_WORKDIR
export NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`

export OMP_NUM_THREADS=1
export  MKL_SERIAL=YES


mpirun  -np  $NPROCS  ./a.out

#mpirun  -n 1 ./a.out