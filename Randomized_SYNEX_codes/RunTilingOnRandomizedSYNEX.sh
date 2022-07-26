#!/bin/bash -l
#SBATCH --job-name=TilingRandomizedSYNEX
#SBATCH --partition=bigmem # quiet
#SBATCH --ntasks=32
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=3500MB # 1500MB
#SBATCH --mail-type=ALL

# Activate conda env to use mpi4py etc
eval "$(/soft/anaconda3/bin/conda shell.bash hook)"
module load openmpi/3.1.6
conda activate fast-mpi4py

# Exort the python path etc
export PATH=/soft/anaconda3/bin:/soft/anaconda3/condabin:/home/baird/.local/bin:$PATH
export PYTHONPATH=/soft/anaconda3/lib:/home/baird/.local/lib:/usr/lib:$PYTHONPATH

# Run the job
time mpirun -np $SLURM_NTASKS python3 ${SLURM_SUBMIT_DIR}/SYNEX/SYNEX_TileRandomized.py > ${SLURM_SUBMIT_DIR}/TiledRandomizedSYNEX.txt

# happy end
exit 0
