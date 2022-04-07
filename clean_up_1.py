import glob
from SYNEX.SYNEX_Utils import SYNEX_PATH
H5_Files = glob.glob(SYNEX_PATH+"/inference_data/Randomized_SYNEX2/Randomized_*_raw.h5")
Total_Json_Files = glob.glob(SYNEX_PATH+"/inference_param_files/Randomized_SYNEX2/Randomized_*.json")

Json_Files = []
for H5_File in H5_Files:
    Json_File = H5_File.split("_raw.h5")[0]+".json"
    Json_File = Json_File.split("inference_data")[0] + "inference_param_files" + Json_File.split("inference_data")[-1]
    Json_Files+=[Json_File]

with open('clean_up_list.txt', 'w') as f:
    for Json_File in Json_Files: f.write(Json_File.split("/")[-1]+"\n") ### No substructure in data savefile locations on cluster...
