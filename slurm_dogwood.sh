#!/bin/bash
#SBATCH -p debug_queue
#SBATCH -J 152Sm0.5
#SBATCH -o out_%A.out
#SBATCH -N 2
#SBATCH --ntasks-per-node=44 
#SBATCH --time=00-04:00:00   #format days-hh:mm:ss

ulimit -s unlimited
mpirun python ./run_lpfam.py


