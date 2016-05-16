#!/bin/bash -l
#SBATCH -p batch 
#SBATCH -J morph0CT
#SBATCH -o log_Morph0CT.log
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-type=All
#SBATCH --mail-user=mattyjones@boisestate.edu
#SBATCH -t 48:00:00
#SBATCH --exclusive
#SBATCH --gres=gpu:2

cd /scratch/erjank_project/mattyMorphCT
hoomd -u runMorphCT.py 0
#cp files you'd like to move off of scratch
#mv files that you'd like moved off of scratch
