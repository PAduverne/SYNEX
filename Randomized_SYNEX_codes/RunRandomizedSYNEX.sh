#!/bin/bash -l
#SBATCH --job-name=RandomizedSYNEX
#SBATCH --partition=quiet
#SBATCH --ntasks=32
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=500MB
#SBATCH --mail-type=ALL

# Activate conda env to use mpi4py etc
eval "$(/soft/anaconda3/bin/conda shell.bash hook)"
module load openmpi/3.1.6
conda activate fast-mpi4py

# Export PATH and PYTHON_PATH
export PATH=/soft/anaconda3/condabin:/home/baird/.local/bin:$PATH
export PATH=/soft/anaconda3/bin:$PATH
export PYTHONPATH=/home/baird/.local/lib:/usr/lib:$PYTHONPATH
export PYTHONPATH=/soft/anaconda3/lib:$PYTHONPATH

# Get list of json files to run
JSONFILE_LIST=(`python3 /SYNEX_randomized.py`)
len=${#JSONFILE_LIST[@]}

# Lisabeta command to run inference -- change this if submitting elsewhere
LAUNCH_INFER=${SLURM_SUBMIT_DIR}/../lisabeta/lisabeta/inference/ptemcee_smbh.py

# Output file
OUT_FILE = RandomizedSYNEX.txt

# Execute each json one-by-one
for (( i=0; i<$len; i++ ))
do
  JSONFILE=${JSONFILE_LIST[$i]}
  echo $JSONFILE

  time mpirun -np $SLURM_NTASKS python3 ${LAUNCH_INFER} $JSONFILE > ${OUT_FILE}
  echo "Copying data to apcssh"
  echo "..."
  HFILETMP=${JSONFILE/%.json/.h5}
  HFILE=${HFILETMP/inference_param_files/inference_data}
  echo $i, $len, $JSONFILE, $HFILE
  scp -i ~/.ssh/id_rsa $JSONFILE baird@apcssh.in2p3.fr:/home/baird/
  scp -i ~/.ssh/id_rsa $HFILE baird@apcssh.in2p3.fr:/home/baird/
  echo "DONE"
  echo " "
  echo "Deleting json and h5 files on cluster"
  echo "..."
  rm {$JSONFILE,$HFILE}
  echo " "
  echo " ......... DONE ......... "
  echo " "
  echo " "
done

echo " -- FINISHED -- "

# happy end
exit 0
