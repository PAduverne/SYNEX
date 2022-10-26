import numpy as np
import os,sys

# What should we do?
DO_INFERENCE_FROM_DANTE = False
DO_INFERENCE_FROM_SSH = True
DO_INFERENCE_TO_SSH = False
DO_TILED_FROM_SSH = False




####
#
# code removed around line 60 -- if new code breaks, reimplement the block below.
#
###

# JsonFiles = ""
# DataFiles = ""
# NFiles=0
# for FileEnd in FileEnds:
#     Filetmp = FileBeginning + FileEnd
#     for Cut in Cuts:
#         NFiles+=1
#         Filetmp2 = Filetmp + Cut
#         JsonFiles = JsonFiles + Filetmp2 + ","
#         DataFiles = DataFiles + Filetmp2 + ","
# JsonFiles=JsonFiles[:-1]
# DataFiles=DataFiles[:-1]

####
#
# code removed around line 60 -- if new code breaks, reimplement the block above.
#
###





# begin transfering
if DO_INFERENCE_FROM_DANTE:
    # Get all json files
    JsonFiles=list(os.popen('ssh baird@apcdante:/home/baird/SYNEX/ "ls inference_param_files"').read())
    JsonFiles = str("".join(JsonFiles)).replace('\n','')
    JsonFiles = ",".join(list(JsonFiles.split(".json")))

    # Get all h5 files
    DataFiles=list(os.popen('ssh baird@apcdante:/home/baird/SYNEX/ "ls inference_data"').read())
    DataFiles = str("".join(DataFiles)).replace('\n','')
    DataFiles = ",".join(list(DataFiles.split(".h5")))

    # transfer jsons
    os.system("scp baird@apcdante:/home/baird/SYNEX/inference_param_files/\{" + JsonFiles[:-1] + "\}.json inference_param_files/Randomized_SYNEX2/")

    # transfer data
    os.system("scp baird@apcdante:/home/baird/SYNEX/inference_data/\{" + DataFiles[:-1] + "\}.h5 inference_data/Randomized_SYNEX2/")
elif DO_INFERENCE_FROM_SSH:
    # Get all json files
    JsonFiles=list(os.popen('ssh baird@apcssh.in2p3.fr "ls *.json"').read())
    JsonFiles = str("".join(JsonFiles)).replace('\n','')
    JsonFiles = ",".join(list(JsonFiles.split(".json")))

    # Get all h5 files
    DataFiles=list(os.popen('ssh baird@apcssh.in2p3.fr "ls *.h5"').read())
    DataFiles = str("".join(DataFiles)).replace('\n','')
    DataFiles = ",".join(list(DataFiles.split(".h5")))

    # transfer jsons
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/\{" + JsonFiles[:-1] + "\}.json inference_param_files/Randomized_SYNEX2/")

    # transfer data
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/\{" + DataFiles[:-1] + "\}.h5 inference_data/Randomized_SYNEX2/")
elif DO_INFERENCE_TO_SSH:
    # Set up params for filenames we want to transfer over
    FileBeginning = "Randomized_angles_spins_MRat_" #
    Cuts = ["_0cut","_5hr","_10hr","_1d","_3d","_1wk", "_2wk", "_3wk", "_1mon"] # ["_5hr","_10hr","_1d","_3d","_1wk", "_2wk", "_3wk", "_1mon"] #
    FileEnds = [str(ii) for ii in range(1,130+1)]

    # Iterate through each param, adding to the string
    Files = [FileBeginning + FileEnd + Cut for FileEnd in FileEnds for Cut in Cuts]
    NFiles=len(Files)
    bunch_limit=300

    # Transfer in bunches otherwise we reach memory issues (idk why)
    if NFiles>bunch_limit:
        file_lump=int(NFiles//bunch_limit+1)
        for ii in range(file_lump):
            start=bunch_limit*ii
            end=bunch_limit*(ii+1) if bunch_limit*(ii+1)<NFiles else NFiles
            # transfer jsons
            os.system("scp inference_param_files/Randomized_SYNEX2/{" + ",".join(Files[start:end]) + "}.json baird@apcssh.in2p3.fr:/home/baird/")
            # transfer data
            os.system("scp inference_data/Randomized_SYNEX2/{" + ",".join(Files[start:end]) + "}.h5 baird@apcssh.in2p3.fr:/home/baird/")
    else:
        # join all together
        TransferFiles = ",".join(Files)
        # transfer jsons
        os.system("scp inference_param_files/Randomized_SYNEX2/{" + TransferFiles + "}.json baird@apcssh.in2p3.fr:/home/baird/")
        # transfer data
        os.system("scp inference_data/Randomized_SYNEX2/{" + TransferFiles + "}.h5 baird@apcssh.in2p3.fr:/home/baird/")
elif DO_TILED_FROM_SSH:
    # Get all source files
    source_saves=list(os.popen('ssh baird@apcssh.in2p3.fr "ls sources"').read())
    source_saves = str("".join(source_saves)).replace('\n','')
    print(len(list(source_saves.split(".dat"))),"source saves...")
    source_saves = ",".join(list(source_saves.split(".dat")))

    # Get all telescope files
    telescope_saves=list(os.popen('ssh baird@apcssh.in2p3.fr "ls telescopes"').read())
    telescope_saves = str("".join(telescope_saves)).replace('\n','')
    len_tels=len(list(telescope_saves.split(".dat")))
    print(len_tels, len(telescope_saves[0:int(len_tels//2)]), len(telescope_saves[int(len_tels//2):len_tels]),"telescope saves...")
    telescope_saves = list(telescope_saves.split(".dat"))
    telescope_saves_1 = ",".join(telescope_saves[0:int(len_tels//4)])
    telescope_saves_2 = ",".join(telescope_saves[int(len_tels//4):2*int(len_tels//4)])
    telescope_saves_3 = ",".join(telescope_saves[2*int(len_tels//4):3*int(len_tels//4)])
    telescope_saves_4 = ",".join(telescope_saves[3*int(len_tels//4):len_tels])
    print(len_tels, len(telescope_saves_1), len(telescope_saves_2), len(telescope_saves_3), len(telescope_saves_4))

    # transfer source savefiles
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/sources/\{" + source_saves[:-1] + "\}.dat Saved_Source_Dicts/Randomized_SYNEX2/")

    # transfer telescope savefiles
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/telescopes/\{" + telescope_saves_1 + "\}.dat Saved_Telescope_Dicts/Randomized_SYNEX2/")
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/telescopes/\{" + telescope_saves_2 + "\}.dat Saved_Telescope_Dicts/Randomized_SYNEX2/")
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/telescopes/\{" + telescope_saves_3 + "\}.dat Saved_Telescope_Dicts/Randomized_SYNEX2/")
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/telescopes/\{" + telescope_saves_4[:-1] + "\}.dat Saved_Telescope_Dicts/Randomized_SYNEX2/")
