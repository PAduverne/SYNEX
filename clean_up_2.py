from SYNEX.SYNEX_Utils import SYNEX_PATH
import os

with open('clean_up_list.txt') as f:
    JsonFilesToCleanUp = f.readlines()

for JsonFile in JsonFilesToCleanUp:
    filename = SYNEX_PATH+"/inference_param_files/Randomized_SYNEX2/"+JsonFile[:-1] ### Because there is "\n at the end of each line"
    if os.path.exists(filename):
        print(filename, "found!")
        # os.remove(filename)
    else:
        print(filename, "does not exist...")
