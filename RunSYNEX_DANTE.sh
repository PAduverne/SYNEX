#!/bin/bash -l

#SBATCH --job-name=SYNEX_TestRun

# 80 MPI tasks in total
# Stallo has 16 or 20 cores/node and therefore we take
# a number that is divisible by both
#SBATCH --ntasks=16

# run for a short time in format: d-hh:mm:ss
# #SBATCH --time=0-00:30:00

# 500MB memory per core - but why do we specify this here?
# this is a hard limit. This will give you --nodes x --ntasks-per-node x --cpus-per-task x 500MB total memory...
# Once your test job is done, check the slurm-xxx.out file for information of memory, cores used etc so you can tune
# your requests in future runs and not clog up the queues.
#SBATCH --mem-per-cpu=500MB

# turn on all mail notification
#SBATCH --mail-type=ALL

# you may not place bash commands before the last SBATCH directive

# Exort the python path etc
export PATH=/soft/openmpi/openmpi-1/bin/:/usr/lib/python3/dist-packages:/home/baird/.local/lib/python3.7/site-packages:/home/baird/.local/bin:$PATH
export PYTHONPATH=/home/baird/:/home/baird/SYNEX:/home/baird/SYNEX/lisabeta:/home/baird/SYNEX/lisabeta/lisabeta:/home/baird/lisabeta/lisabeta/struct:/home/baird/SYNEX/lisabeta/lisabeta/lisa:/home/baird/SYNEX/acor:/home/baird/SYNEX/ptemcee:/home/baird/.local/bin:/home/baird/.local/lib/python3.7/site-packages:/usr/lib/python3/dist-packages:/home/baird/SYNEX/SYNEX:$PYTHONPATH
export LD_LIBRARY_PATH=/home/baird/:/usr/local/lib:/soft/openmpi/openmpi-1/lib/:/soft/openmpi/openmpi-1/lib/openmpi/:/usr/local/openmpi/lib:/home/baird/.local/bin:/usr/local/openmpi/lib/openmpi/:/usr/lib/python3/dist-packages/:/home/baird/.local/lib/python3.7/site-packages:$LD_LIBRARY_PATH

# define and create a unique scratch directory
SCRATCH_DIRECTORY=/scratch/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

# we copy everything we need to the scratch directory
# ${SLURM_SUBMIT_DIR} points to the path where this script was submitted from
# do we need to copy just the SYNEX2.py script or the whole SYNEX folder accross?
cp -R ${SLURM_SUBMIT_DIR}/SYNEX ${SCRATCH_DIRECTORY}

# we execute the job and time it
time mpirun -np $SLURM_NTASKS python3 ./SYNEX/SYNEX2.py > my_output

# after the job is done we copy our output back to $SLURM_SUBMIT_DIR
cp -R ${SCRATCH_DIRECTORY}/my_output ${SLURM_SUBMIT_DIR}

# we step out of the scratch directory and remove it
cd ${SLURM_SUBMIT_DIR}
# rm -rf ${SCRATCH_DIRECTORY}

# happy end
exit 0
