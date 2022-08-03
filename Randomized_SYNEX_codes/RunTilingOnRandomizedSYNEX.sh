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

# Directory to check for existing data on cluster
CLUST_JSON_DIR=${SYNEX_DIR}/inference_param_files/

# If Is inference_param_files directory empty, grab 10 more sources to transfer
if [ ! "$(ls -A $CLUST_JSON_DIR)" ] ### JSON DIR lust be empty for this to work.
then
  # Some useful commands
  ssh_home=baird@apcssh.in2p3.fr
  ssh_home2=baird@apcssh.in2p3.fr:/home/baird/

  # Get list of all systems
  H5FILE_ssh_LIST=($(ssh -i ~/.ssh/id_rsa ${ssh_home} 'ls Randomized_*.h5'))
  len_H5FILE_ssh_LIST=${#H5FILE_ssh_LIST[@]}

  # check things are ok
  echo "Test test test:" $len_H5FILE_ssh_LIST ${H5FILE_ssh_LIST[1]}

  # Get list of all source save files
  SAVEFILE_LIST=($(ssh -i ~/.ssh/id_rsa ${ssh_home} 'ls sources/'))

  # Initiate container for files to transfer across
  TRANSFERFILES=()

  # Remove from lisabeta data list all completed systems
  for H5FILE in ${H5FILE_ssh_LIST[@]} ; do
    SaveFileTMP="$H5FILE" | sed 's/.h5//'
    # SaveFileTMP=${H5FILE/%.h5/.dat}
    for file in ${SAVEFILE_LIST[@]} ; do
      file2="$file" | sed 's/.dat//'
      if [[ $file2 == $SaveFileTMP ]] ; then
        break
      fi
    done
    TRANSFERFILES+=("$SaveFileTMP")
    # Stop when we have 10 sources to transfer
    if [[ ${#TRANSFERFILES[@]} -eq 10 ]] ; then
      break
    fi
  done

  # Transfer all files
  scp -i ~/.ssh/id_rsa ${ssh_home2}/${TRANSFERFILES[@]}.h5 ${SYNEX_DIR}/inference_data/
  scp -i ~/.ssh/id_rsa ${ssh_home2}/${TRANSFERFILES[@]}.json ${SYNEX_DIR}/inference_param_files/

  # # Now set flags for tiling command options
  # USETRACK=False # This will be created using new sources
  # MAKE_SOURCES=True # Make new sources for h5 data
fi

# gwemopt command to run tiling -- change this if submitting elsewhere
LAUNCH_TILING=${SYNEX_DIR}/Randomized_SYNEX_codes/SYNEX_TileRandomized.py

# Output file
OUT_FILE=TiledRandomizedSYNEX.txt

# Run the job
time mpirun -np $SLURM_NTASKS python3 ${LAUNCH_TILING} > $OUT_FILE

# happy end
exit 0
