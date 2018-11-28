#!/bin/bash -l
#SBATCH -p batch 
#SBATCH -J <INPUTXML>
#SBATCH -o log_<INPUTXML>.log
#SBATCH -N 1
#SBATCH -n 7
#SBATCH -t 48:00:00

source activate morphct

hoomd -u <INPUTPAR>
