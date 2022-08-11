#!/bin/bash -l
#SBATCH --job-name=TilingRandomizedSYNEX
#SBATCH --partition=bigmem
#SBATCH --ntasks=32
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=3500MB # 40 source H5 data files, otherwise out of memory
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
### Enough memory for 40 systems at a time - see what needs doing ###
###
##########

# Directory to check for existing data on cluster
CLUST_JSON_DIR=${SYNEX_DIR}/inference_param_files/

### Include an else if the JSON dir isn't empty in which caase we send the options
### to the tiling function to use the track file and complete the stuff still
### there. Can we maybe also include a check if the fil is empty? In which case
### we we proceed to Tiling without source creation etc?

# If Is inference_param_files directory empty, grab 10 more sources to transfer
if [ ! "$(ls -A $CLUST_JSON_DIR)" ] ### JSON DIR lust be empty for this to work.
then
  # Some useful directories and files
  ssh_home=baird@apcssh.in2p3.fr
  ssh_home2=baird@apcssh.in2p3.fr:/home/baird
  rsa=~/.ssh/id_rsa

  # Get list of all systems
  H5FILE_ssh_LIST=($(ssh -i $rsa ${ssh_home} 'ls Randomized_*.h5'))
  len_H5FILE_ssh_LIST=${#H5FILE_ssh_LIST[@]}

  # Get list of all source save files
  SAVEFILE_LIST=($(ssh -i $rsa ${ssh_home} 'ls sources/'))

  # Initiate container for files to transfer across
  TRANSFERFILES=()

  # Remove from lisabeta data list all completed systems
  for H5FILE in ${H5FILE_ssh_LIST[@]} ; do
    ADDH5=1
    for file in ${SAVEFILE_LIST[@]} ; do
      if [[ $file == ${H5FILE/%.h5/.dat} ]] ; then
        ADDH5=0
        break
      fi
    done

    # Add H5 file if we are allowed
    if [[ ADDH5 -eq 1 ]] ; then
      TRANSFERFILES+=("${H5FILE/%.h5}")
    fi

    # Stop when we have 40 source data files to transfer
    if [[ ${#TRANSFERFILES[@]} -eq 40 ]] ; then
      break
    fi
  done

  # Transfer all files
  FILES=($(echo "${TRANSFERFILES[@]}" | tr ' ' ','))
  scp -i $rsa ${ssh_home2}/{${FILES[@]}}.h5 ${SYNEX_DIR}/inference_data/
  scp -i $rsa ${ssh_home2}/{${FILES[@]}}.json ${SYNEX_DIR}/inference_param_files/

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

# Copy everything back to APCSSH
SourceSaves=($(ls ${SYNEX_DIR}/Saved_Source_Dicts/Randomized_*.dat))
for SourceSave in ${SourceSaves[@]} ; do
  scp -i $rsa $SourceSave baird@apcssh.in2p3.fr:/home/baird/sources/
  # scp -i $rsa $SourceSave ${ssh_home2}/sources/ # ???
done

TelesSaves=($(ls ${SYNEX_DIR}/Saved_Telescope_Dicts/Randomized_*.dat))
# len=${#TelesSaves[@]}
for TelesSave in ${TelesSaves[@]} ; do
  scp -i $rsa $TelesSave baird@apcssh.in2p3.fr:/home/baird/telescopes/
  # scp -i $rsa $TelesSave ${ssh_home2}/telescopes/ # ???
done

# Clear folders
rm ${SYNEX_DIR}/Saved_Telescope_Dicts/Randomized_*
rm ${SYNEX_DIR}/Saved_Source_Dicts/Randomized_*
rm ${SYNEX_DIR}/inference_param_files/Randomized_*
rm ${SYNEX_DIR}/inference_data/Randomized_*

# happy end
exit 0
