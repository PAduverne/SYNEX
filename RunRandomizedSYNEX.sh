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

# Exort the python path etc
export PATH=/soft/anaconda3/bin:/soft/anaconda3/condabin:/home/baird/.local/bin:$PATH
export PYTHONPATH=/soft/anaconda3/lib:/home/baird/.local/lib:/usr/lib:$PYTHONPATH

# we execute each job and time one-by-one
JSONFILE_LIST=(`python3 ${SLURM_SUBMIT_DIR}/SYNEX/SYNEX_randomized.py`)
len=${#JSONFILE_LIST[@]}
for (( i=0; i<$len; i++ ))
do
  JSONFILE=${JSONFILE_LIST[$i]}
  echo $JSONFILE
  time mpirun -np $SLURM_NTASKS python3 ${SLURM_SUBMIT_DIR}/SYNEX/lisabeta/lisabeta/inference/ptemcee_smbh.py $JSONFILE > ${SLURM_SUBMIT_DIR}/RandomizedSYNEX.txt
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
