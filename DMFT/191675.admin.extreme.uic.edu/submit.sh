#!/bin/bash
#PBS -l nodes=1:ppn=20,walltime=10:00:00
#PBS -N DMFT 
#PBS -q edu_shared
#PBS -m abe
#PBS -M emasco2@uic.edu

module load compilers/intel-2015
module load tools/gsl-1.16-intel

cd $PBS_O_WORKDIR
mkdir $PBS_JOBID
cd $PBS_JOBID
cp ../* .

echo mpirun -machinefile $PBS_NODEFILE -np $PBS_NP > para_com.dat
python dmft.py > out 2> err

