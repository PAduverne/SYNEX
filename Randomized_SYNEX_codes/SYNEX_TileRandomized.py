import numpy as np
import os,sys,glob
from astropy.time import Time
import time

# Define SYNEX_PATH variable, add to path if not already there
#  -- move to shell script?
try:
    from SYNEX.SYNEX_Utils import SYNEX_PATH
except:
    SYNEX_PATH=os.popen('pwd').read().split("SYNEX")[0]+"SYNEX"
    sys.path.insert(1, SYNEX_PATH)

from SYNEX import SYNEX_Detectors as SYDs
from SYNEX import SYNEX_Sources as SYSs
from SYNEX import SYNEX_Utils as SYU

# mpi stuff
try:
    from mpi4py import MPI
except ModuleNotFoundError:
    MPI = None

##############################
#
# Flags to be moved to user input args later
#
# Set verbosity inside SYNEX (making objects etc)
verbose = False
# Load combinations from track file or clone from scratch?
USETRACK=False
# Make sources from lisabeta h5 files i.e. do savefiles already exist?
MAKE_SOURCES = True
#
#
##############################

# Using MPI or not
if MPI is not None:
    MPI_size = MPI.COMM_WORLD.Get_size()
    MPI_rank = MPI.COMM_WORLD.Get_rank()
    comm = MPI.COMM_WORLD
else:
    MPI_size=1
    MPI_rank=0

# Make source save files from inference data ?
if MAKE_SOURCES==True:
    # Force off search in track file for systems
    USETRACK=False
    # Wild card string to search for inference files (assume no raw files)
    SearchString=SYNEX_PATH+"/inference_data/Randomized_angles_spins_MRat_*.h5"
    # Split up sources across nodes
    if MPI_rank==0:
        # Master node only gets list of all filenames
        source_infer_datafiles=glob.glob(SearchString)
        Nvals=len(source_infer_datafiles)
        nPerCore=int(Nvals//MPI_size)
        # Start scattering -- variable Nvals so can't use mpi.scatter ...
        ii_start=0
        ii_end=nPerCore if nPerCore<Nvals else Nvals
        source_infer_data = source_infer_datafiles[ii_start:ii_end] # master node
        for core_ii in range(1,MPI_size):
            ii_start=ii_end
            ii_end=nPerCore*(core_ii+1) if nPerCore*(core_ii+1)<Nvals else Nvals
            if core_ii==MPI_size-1: ii_end=Nvals # Extras from rounding errors
            comm.send(source_infer_datafiles[ii_start:ii_end], dest=core_ii)
    else:
        source_infer_data = comm.recv(source=0)

    # Create sources and explicitly save since writting
    # is restricted to save disk space on cluster
    if MPI_rank==0: print("Making source savefiles from inference data...")
    SourceSaveFolder=SYNEX_PATH+"/Saved_Source_Dicts/"
    for FileName in source_infer_data:
        ExName=FileName.split("/inference_data/")[-1].replace(".h5",".dat")
        ExName=SourceSaveFolder+ExName
        Merger=SYU.GetSourceFromLisabetaData(FileName,
                            **{"ExistentialFileName":ExName,"verbose":verbose})
        # Explicit save
        Merger.ExistentialCrisis()

# Line up all nodes before continuing
comm.Barrier()
if MPI_rank==0: print("Done.")

# Telescope base args -- move this to a default file and change by user input?
t0 = '2034-01-01T00:00:00.00' # YYYY-MM-DDTHH:mm:SS.MS 01/01/2034
t = Time(t0, format='isot', scale='utc').gps
ExName=SYNEX_PATH+"/Saved_Telescope_Dicts/Athena_base.dat"
Athena_kwargs={"ExistentialFileName":ExName,
                "verbose":verbose,
                "telescope":"Athena",
                "Tobs":np.array([0.,9.]), # pairs of [Tstart,Tend] (in DAYS).
                "tilesType" : "moc",
                "timeallocationType" : "powerlaw",
                "scheduleType" : "greedy",
                "doCalcTiles" : False,
                "Ntiles" : None,
                "frozenAthena" : False, # False,
                "exposuretime" : 10000.,
                "min_observability_duration" : None,
                "inc" : 60.,
                "MeanRadius" : 750000000.,
                "semi_maj" : 750000000.,
                "eccentricity" : 0.4,
                "ArgPeriapsis" : 20.,
                "AscendingNode" : -10.,
                "phi_0" : 10.,
                "period" : 90.,
                "gps_science_start" : t,
                "mission_duration" : 2.,
                "filt_change_time" : 0.,
                "overhead_per_exposure" : 0.,
                "latitude" : 0.,  ### None if we want a telesscopic orbit?
                "longitude" : 0., ### None if we want a telesscopic orbit?
                "elevation" : 0., ### None if we want a telesscopic orbit?
                "slew_rate" : 1.,
                "horizon" : 0.,   ### None if we want a telesscopic orbit?
                "doMinimalTiling" : True,
                "readout" : 0.0001,
                "doSingleExposure" : True,
                "iterativeOverlap" : 0.,
                "maximumOverlap" : 1.0,
                "sat_sun_restriction" : 5.,
                "sat_earth_constraint" : 5.,
                "sat_moon_constraint" : 5.,
               } ### "doPerturbativeTiling" : True ?

# Get all source savefiles to test
CutsToTest = ["5hr","10hr","1d","3d","1wk","2wk","3wk","1mon"]
SearchString=SYNEX_PATH+"/Saved_Source_Dicts/Randomized_angles_spins_MRat_*_"
SourceExNames = sorted([File for c in CutsToTest for File in
    glob.glob(SearchString+c+".dat")]) ### Makes "10" go in front of "1"...

# What to vary when cloning base detector args ---- to be user input
T_obs_array = [np.array([0.,1.]),np.array([0.,2.]),
                   np.array([0.,3.]),np.array([0.,4.])]
cloning_params={"Tobs":T_obs_array}

# Use tracking file or clone new stuff?
if USETRACK:
    # Systems will be loaded from file and cloning params etc ignored
    CloningTrackFile=SYNEX_PATH+"/TileTrackFiles/ProgressTrackFile.txt"
else:
    CloningTrackFile=None

# Tile all combinations of cloning params and sources
SaveInSubFile=None
SaveFileCommonStart="Randomized_angles_spins_MRat"
SourceIndexingByString="MRat_"
T0 = time.time()
SYU.TileSkyArea(CloningTrackFile=CloningTrackFile,sources=SourceExNames,
                    detectors=None,base_telescope_params=Athena_kwargs,
                    cloning_params=cloning_params,SaveInSubFile=SaveInSubFile,
                    SaveFileCommonStart=SaveFileCommonStart,
                    SourceIndexingByString=SourceIndexingByString,verbose=False)
T1 = time.time()
print("Time to finish:",T1-T0,"s")
