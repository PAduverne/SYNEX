#!/bin/bash -l
#SBATCH --job-name=RandomizedSYNEX
#SBATCH --ntasks=32
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=500MB
#SBATCH --mail-type=ALL

# Activate conda env to use mpi4py etc
module purge
eval "$(/soft/anaconda3/bin/conda shell.bash hook)"
module load openmpi/3.1.6
conda activate fast-mpi4py

# Exort the python path etc
export PYTHONPATH=/usr/lib/python3/dist-packages/:$PYTHONPATH

# define and create a unique scratch directory
SCRATCH_DIRECTORY=/scratch/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

# we execute the job and time it
time mpirun -np $SLURM_NTASKS python3 ${SLURM_SUBMIT_DIR}/SYNEX/SYNEX_randomized.py > ${SLURM_SUBMIT_DIR}/RandomizedSYNEX

# we step out of the scratch directory and remove it
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}

# happy end
exit 0 
