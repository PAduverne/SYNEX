#!/bin/bash -l
#SBATCH --job-name=TilingRandomizedSYNEX
#SBATCH --partition=bigmem
#SBATCH --ntasks=32
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=3500MB # 10 systems all Tcut, otherwise errors...
#SBATCH --mail-type=ALL

# Activate conda env to use mpi4py etc
eval "$(/soft/anaconda3/bin/conda shell.bash hook)"
module load openmpi/3.1.6
conda activate fast-mpi4py

# Export PATH and PYTHONPATH
export PATH=/soft/anaconda3/condabin:/home/baird/.local/bin:$PATH
export PATH=/soft/anaconda3/bin:$PATH
export PYTHONPATH=/home/baird/.local/lib:/usr/lib:$PYTHONPATH
export PYTHONPATH=/soft/anaconda3/lib:$PYTHONPATH

# gwemopt command to run tiling -- change this if submitting elsewhere
LAUNCH_INFER=${SLURM_SUBMIT_DIR}/SYNEX/SYNEX_TileRandomized.py

# Output file
OUT_FILE=TiledRandomizedSYNEX.txt

# Run the job
time mpirun -np $SLURM_NTASKS python3 ${LAUNCH_INFER} > $OUT_FILE

# happy end
exit 0
