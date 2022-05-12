import numpy as np
import os,sys

# Set up iter params for filenames
FileBeginning = "Randomized_angles_spins_MRat_" # "Randomized_SYNEX_" # "RedGrid_Mpi_" #
Cuts = ["_0cut","_5hr","_10hr","_1d","_3d","_1wk", "_2wk", "_3wk", "_1mon"] # ["0cut_","4hr_","1d_","1wk_"] # ["_1mon"] #
FileEnds = [str(ii) for ii in range(1+10,20+1)] # ["5","6","7","8","9"] # ["1"] #


# create home dir name for local host
# HomeDirLocal = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/Randomized_SYNEX2/"

# JsonFiles = ""
# DataFiles = ""
# DataRawFiles = ""
# for Cut in Cuts:
#     Filetmp = FileBeginning + Cut
#     for FileEnd in FileEnds:
#         Filetmp2 = Filetmp + FileEnd
#         JsonFiles = JsonFiles + Filetmp2 + ","
#         DataFiles = DataFiles + Filetmp2 + ","
#         DataRawFiles = DataRawFiles + Filetmp2 + "_raw,"

JsonFiles = ""
DataFiles = ""
DataRawFiles = ""
for FileEnd in FileEnds:
    Filetmp = FileBeginning + FileEnd
    for Cut in Cuts:
        Filetmp2 = Filetmp + Cut
        JsonFiles = JsonFiles + Filetmp2 + ","
        DataFiles = DataFiles + Filetmp2 + ","
        DataRawFiles = DataRawFiles + Filetmp2 + "_raw,"

DO_DANTE = False

if DO_DANTE:
    # transfer jsons
    os.system("scp baird@apcdante:/home/baird/SYNEX/inference_param_files/\{" + JsonFiles[:-1] + "\}.json inference_param_files/Randomized_SYNEX2/")

    # transfer data
    os.system("scp baird@apcdante:/home/baird/SYNEX/inference_data/\{" + DataFiles[:-1] + "\}.h5 inference_data/Randomized_SYNEX2/")

    # transfer raw data
    os.system("scp baird@apcdante:/home/baird/SYNEX/inference_data/\{" + DataRawFiles[:-1] + "\}.h5 inference_data/Randomized_SYNEX2/")
else:
    # transfer jsons
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/\{" + JsonFiles[:-1] + "\}.json inference_param_files/Randomized_SYNEX2/")

    # transfer data
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/\{" + DataFiles[:-1] + "\}.h5 inference_data/Randomized_SYNEX2/")

    # transfer raw data
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/\{" + DataRawFiles[:-1] + "\}.h5 inference_data/Randomized_SYNEX2/")
