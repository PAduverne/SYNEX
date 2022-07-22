from SYNEX import SYNEX_Detectors as SYDs
from SYNEX import SYNEX_Sources as SYSs
from SYNEX import SYNEX_Utils as SYU
from SYNEX.SYNEX_Utils import SYNEX_PATH
import glob,os

detector_files=glob.glob(SYNEX_PATH+"/Saved_Telescope_Dicts/Randomized_angles_spins_MRat_*")
detectors=[]
sources=[]

# Check which files are completed
for f in detector_files:
    detector = SYDs.Athena(**{"ExistentialFileName":f,"verbose":False})
    source = detector.detector_source_coverage["source save file"] if detector.detector_source_coverage!=None else None
    if source!=None and os.path.isfile(f):
        detectors.append(f)
        sources.append(source)

# Read tracking file
CloningTrackFile=SYNEX_PATH+"/TileTrackFiles/ProgressTrackFile.txt"
with open(CloningTrackFile, 'r') as f: lines=f.readlines()

# Reduce tracking file contents
for d,s in zip(detectors,sources): lines = [line for line in lines if (not d in line) and (not s in line)]

# Save reduced contents back to file
with open(CloningTrackFile, 'w') as f: f.writelines(lines)

# Drop duplicates
sources = list( dict.fromkeys(sources) )
detectors = list( dict.fromkeys(detectors) )

# Output for transfer in bash function
for f in sources+detectors: print(f)
