#!/bin/bash

#SBATCH --job-name int
#SBATCH --nodes 1
#SBATCH --tasks-per-node 12
#SBATCH --cpus-per-task 1
#SBATCH --mem 32gb
#SBATCH --time 72:00:00
#SBATCH --gpus-per-node v100:1

module purge
module add anaconda3
module add cuda/11.8.0
module add openmpi/5.0.1

#mpirun -np 12 ./lmp_gpu -sf gpu -in interstitial.in -log interstitial.log
mpirun -np 12 ./lmp_gpu -in interstitial.in -log interstitial.log