#!/bin/bash


#SBATCH --job-name=NZ_b20_tanh_res #job name

# 80 MPI tasks in total
# Stallo has 16 or 20 cores/node and therefore we take
# a number that is divisible by both
#SBATCH --ntasks=1
##SBATCH --ntasks-per-node=40
#SBATCH --partition quiet #queue names: quiet, bigmem, debug

# run for five minutes format: d-hh:mm:ss #day, hour, minutes, second
#SBATCH --time=30-00:00:00

# 500MB memory per core : this is a hard limit
#SBATCH --mem-per-cpu=500MB

# turn on all mail notification
#SBATCH --mail-user=duverne@apc.in2p3.fr
#SBATCH --mail-type=FAIL

python3 SMBBH_inference.py \
    --path ${fullpath} \
    --skip-inference \
    --cut ${cut}