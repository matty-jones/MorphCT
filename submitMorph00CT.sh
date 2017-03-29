#!/bin/bash -l
#SBATCH -p batch 
#SBATCH -J P3HTC60
#SBATCH -o log_Morph00CT.log
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-type=All
#SBATCH --mail-user=mattyjones@boisestate.edu
#SBATCH -t 48:00:00
#SBATCH --exclusive
#SBATCH --gres=gpu:2

cd /scratch/erjank_project/mattyMorphCT
on-conda
removePathDuplicates
source activate hoomd1.3
hoomd -u par00.py
