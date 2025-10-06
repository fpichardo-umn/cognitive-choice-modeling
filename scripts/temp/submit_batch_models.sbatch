#!/bin/bash -l
#SBATCH -J batch_fit
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --tmp=32gb
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL
#SBATCH -o log_files/fit/%x_%j.out
#SBATCH -e log_files/fit/%x_%j.err

# Load necessary modules
module load R/4.2.0-rocky8
module load gcc/13.1.0-mptekim
module load cmake/3.29.2-rocky8 
module load boost/1.82.0-gcc-13.1.0-sox6zlc


./scripts/parameter_recovery/run_full_PR.sh -k igt -m dual_orl -s luc --ses 00 -g sing --subjects 1-182 --n-subjects 182 --check_iter 15000 --n-iter 10000 --n-warmup 3000 all