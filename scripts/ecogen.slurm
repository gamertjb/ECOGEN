#!/bin/bash -xv
#SBATCH --job-name=ECOGEN_Test
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=02:00:00
#SBATCH --output=ecogen_output.log
#SBATCH --partition=short

module swap gnu9/9.4.0 gcc/12.2.0
module load openmpi/4.1.5
module load pmix/4.2.9

# Recompile the program if necessary
#make clean
#make -j 64

mpirun -np 64 ./ECOGEN
