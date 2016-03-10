#!/bin/bash
#PBS -l nodes=1:ppn=20,walltime=10:00:00
#PBS -N CTQMC
#PBS -q edu_shared
#PBS -m abe
#PBS -M netid@uic.edu
#PBS -e CTQMC.err
#PBS -o CTQMC.out

module load compilers/intel-2015
module load tools/gsl-1.16-intel

cd $PBS_O_WORKDIR
mkdir $PBS_JOBID
cd $PBS_JOBID
cp ../* .

mpirun -machinefile $PBS_NODEFILE -np $PBS_NP ./ctqmc > out 2> err
