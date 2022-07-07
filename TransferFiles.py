import numpy as np
import os,sys

# Set up iter params for filenames
FileBeginning = "Randomized_angles_spins_MRat_" #
Cuts = ["_5hr","_10hr","_1d","_3d","_1wk", "_2wk", "_3wk", "_1mon"] # ["_0cut","_5hr","_10hr","_1d","_3d","_1wk", "_2wk", "_3wk", "_1mon"] # ["0cut_","4hr_","1d_","1wk_"] #
FileEnds = [str(ii) for ii in range(1,45+1)]

# What should we do?
DO_DANTE = False
TO_SSH = False
RETURN_TILED_FROM_SSH = True

if RETURN_TILED_FROM_SSH:
    # from fabric import Connection
    # source_saves_list = Connection('baird@apcssh.in2p3.fr').run("glob.glob('/home/baird/sources/Randomized_*')") # , hide=True)
    # telescope_saves_list = Connection('baird@apcssh.in2p3.fr').run("glob.glob('/home/baird/detectors/Randomized_*')") # , hide=True)
    source_saves_list=list(os.popen('ssh baird@apcssh.in2p3.fr "ls sources"').read()) ### CAN WE REMOVE THE "" xtirely since we only have one command and no piping to a second command?
    telescope_saves_list=list(os.popen('ssh baird@apcssh.in2p3.fr "ls detectors"').read())
    print("popen check:",len(source_saves_list),np.shape(source_saves_list),source_saves_list[:20])
    source_saves_list = str("".join(source_saves_list)).replace('\n','')
    telescope_saves_list = str("".join(telescope_saves_list)).replace('\n','')
    print("join check:",len(source_saves_list),np.shape(source_saves_list),source_saves_list[:40])
    source_saves = ",".join(list(source_saves_list.split(".dat")))
    telescope_saves = ",".join(list(telescope_saves_list.split(".dat")))
    print("rejoin check:",len(source_saves),np.shape(source_saves),source_saves[:20])
    print(source_saves)
    # NFiles=len(source_saves_list) ## In case we need this later...
    # source_saves=""
    # telescope_saves = ""
    # for el in source_saves_list: source_saves += el+"," ## Not sure why but using ",".join(list) throws an error when doing the system execution of scp... So we gotta do things manually.
    # for el in telescope_saves_list: telescope_saves += el+"," ## Not sure why but using ",".join(list) throws an error when doing the system execution of scp... So we gotta do things manually.
else:
    JsonFiles = ""
    DataFiles = ""
    DataRawFiles = ""
    NFiles=0
    for FileEnd in FileEnds:
        Filetmp = FileBeginning + FileEnd
        for Cut in Cuts:
            NFiles+=1
            Filetmp2 = Filetmp + Cut
            JsonFiles = JsonFiles + Filetmp2 + ","
            DataFiles = DataFiles + Filetmp2 + ","
            DataRawFiles = DataRawFiles + Filetmp2 + "_raw,"


if DO_DANTE:
    # transfer jsons
    os.system("scp baird@apcdante:/home/baird/SYNEX/inference_param_files/\{" + JsonFiles[:-1] + "\}.json inference_param_files/Randomized_SYNEX2/")

    # transfer data
    os.system("scp baird@apcdante:/home/baird/SYNEX/inference_data/\{" + DataFiles[:-1] + "\}.h5 inference_data/Randomized_SYNEX2/")

    # transfer raw data
    os.system("scp baird@apcdante:/home/baird/SYNEX/inference_data/\{" + DataRawFiles[:-1] + "\}.h5 inference_data/Randomized_SYNEX2/")
elif TO_SSH:
    bunch_limit=100
    if NFiles>bunch_limit:
        JsonFiles = JsonFiles[:-1].split(",")
        DataFiles = DataFiles[:-1].split(",")
        file_lump=int(NFiles//bunch_limit+1)
        for ii in range(file_lump):
            start=bunch_limit*ii
            end=bunch_limit*(ii+1) if bunch_limit*(ii+1)<NFiles else NFiles
            # transfer jsons
            os.system("scp inference_param_files/Randomized_SYNEX2/{" + ",".join(JsonFiles[start:end]) + "}.json baird@apcssh.in2p3.fr:/home/baird/")
            # transfer data
            # os.system("scp inference_data/Randomized_SYNEX2/{" + ",".join(DataFiles[start:end]) + "}.h5 baird@apcssh.in2p3.fr:/home/baird/")
    else:
        # transfer jsons
        os.system("scp inference_param_files/Randomized_SYNEX2/{" + JsonFiles[:-1] + "}.json baird@apcssh.in2p3.fr:/home/baird/")
        # transfer data
        os.system("scp inference_data/Randomized_SYNEX2/{" + DataFiles[:-1] + "}.h5 baird@apcssh.in2p3.fr:/home/baird/")
elif RETURN_TILED_FROM_SSH:
    # transfer source savefiles
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/sources/\{" + source_saves[:-1] + "\}.dat Saved_Source_Dicts/Randomized_SYNEX2/")

    # transfer telescope savefiles
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/detectors/\{" + telescope_saves[:-1] + "\}.dat Saved_Telescope_Dicts/Randomized_SYNEX2/")
else:
    # transfer jsons
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/\{" + JsonFiles[:-1] + "\}.json inference_param_files/Randomized_SYNEX2/")

    # transfer data
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/\{" + DataFiles[:-1] + "\}.h5 inference_data/Randomized_SYNEX2/")

    # transfer raw data
    os.system("scp baird@apcssh.in2p3.fr:/home/baird/\{" + DataRawFiles[:-1] + "\}.h5 inference_data/Randomized_SYNEX2/")
