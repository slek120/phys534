#!/bin/bash
#PBS -N espressojob
#PBS -l nodes=1:ppn=20
#PBS -l walltime=3:00:00
#PBS -j oe
#PBS -q edu_shared 


module load apps/espresso-5.3.0-intel
NN=$(cat ${PBS_NODEFILE} | wc -l)
cd $PBS_O_WORKDIR
mkdir $PBS_JOBID
cd $PBS_JOBID
cp ../* .
mpirun -machinefile $PBS_NODEFILE -np $NN pw.x -i graphene.scf.in > out 2> err
