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

# Export SYNEX directory as seperate variable - may need to be changed
# depending on your installation
export SYNEX_DIR=~/SYNEX

##########
###
### Enough memory for 90 systems at a time - see what needs doing ###
###
##########

# Is inference_param_files directory empty?
JSONFILE_LIST_CLUST=(ls ${SYNEX_DIR}/inference_param_files/Randomized_*.json)
len_json_list=${#JSONFILE_LIST_CLUST[@]}

# If so, then grab up to 10 more sources
if [len_json_list -eq 0]
then
  # Some useful commands
  SSH_COMM = "scp -i ~/.ssh/id_rsa baird@apcssh.in2p3.fr:/home/baird"
  SCP_COMM = "scp -i ~/.ssh/id_rsa baird@apcssh.in2p3.fr:/home/baird"

  # Get list of all systems
  H5FILE_ssh_LIST=($SSH_COMM/ 'ls Randomized_*.h5')
  len_H5FILE_ssh_LIST=${#H5FILE_ssh_LIST[@]}

  # Get list of all source save files
  SAVEFILE_LIST=($SSH_COMM/ 'ls sources/*.dat')

  # Initiate container for files to transfer across
  FILES_TO_TRANSFER=()

  # Remove from lisabeta data list all completed systems
  for H5FILE in ${H5FILE_ssh_LIST} ; do
    SaveFileTMP = "$H5FILE" | sed 's/.h5//'
    # SaveFileTMP=${H5FILE/%.h5/.dat}
    for file in ${SAVEFILE_LIST[@]} ; do
      file2 = "$file" | sed 's/.dat//'
      if [[ $file2 == $SaveFileTMP ]] ; then
        break
      fi
    done
    FILES_TO_TRANSFER+=("$SaveFileTMP")
    # Stop when we have 10 sources to transfer
    if [[ ${#FILES_TO_TRANSFER[@]} -eq 10 ]] ; then
      break
    fi
  done

  # Transfer all files
  SCP_COMM/${FILES_TO_TRANSFER[@]}.h5 ${SYNEX_DIR}/inference_data/
  SCP_COMM/${FILES_TO_TRANSFER[@]}.json ${SYNEX_DIR}/inference_param_files/

  # # Now set flags for tiling command options
  # USETRACK=False # This will be created using new sources
  # MAKE_SOURCES=True # Make new sources for h5 data
fi

# gwemopt command to run tiling -- change this if submitting elsewhere
LAUNCH_INFER=${SLURM_SUBMIT_DIR}/SYNEX/SYNEX_TileRandomized.py

# Output file
OUT_FILE=TiledRandomizedSYNEX.txt

# Run the job
time mpirun -np $SLURM_NTASKS python3 ${LAUNCH_INFER} > $OUT_FILE

# happy end
exit 0
