#!/bin/bash
#PBS -l nodes=1:ppn=16,walltime=10:00:00
#PBS -N CTQMC_SL
#PBS -q edu_shared
#PBS -m abe
#PBS -M emasco2@uic.edu
#PBS -e CTQMC_SL.err
#PBS -o CTQMC_SL.out

module load compilers/intel-2015
module load tools/gsl-1.16-intel

cd $PBS_O_WORKDIR
mkdir $PBS_JOBID
cd $PBS_JOBID
cp ../* .

mpirun -machinefile $PBS_NODEFILE -np $PBS_NP ./ctqmc > out 2> err
