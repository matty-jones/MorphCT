#!/bin/bash -l
#SBATCH -p batch 
#SBATCH -J JOBNAME
#SBATCH -o OUTFILE
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=All
#SBATCH --mail-user=CHANGEME@boisestate.edu
#SBATCH -t 12:00:00
#SBATCH --exclusive
#SBATCH --gres=gpu:1

cd /scratch/${USER}
#srun --kill-on-bad-exit --mpi=pmi2 --cpu_bind=map_cpu:8 -n1 hoomd 2>&1 b.py --mode=gpu --gpu=0 | tee b.o &
srun --kill-on-bad-exit --mpi=pmi2 --cpu_bind=map_cpu:8 -n1 hoomd myfile.py --mode=gpu --gpu=0
#cp files you'd like to move off of scratch
#mv files that you'd like moved off of scratch
