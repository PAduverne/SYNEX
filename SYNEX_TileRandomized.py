import numpy as np
import astropy.constants as const
import astropy.units as u
from astropy.cosmology import WMAP9 as cosmo
from astropy.utils.data import download_file
from astropy.io import fits
from astropy.time import Time

import lisabeta.lisa.lisatools as lisatools
import lisabeta.utils.plotutils as plotutils

from SYNEX import SYNEX_Detectors as SYDs
from SYNEX import SYNEX_Sources as SYSs
from SYNEX import SYNEX_Utils as SYU
from SYNEX.SYNEX_Utils import SYNEX_PATH

from numpy.random import rand

import time
import os
import json
import glob
import copy
import healpy as hp
from gwemopt import utils as gou
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pylab as pylab
from SYNEX.SYNEX_Utils import pylab_params
pylab.rcParams.update(pylab_params)
try:
    mpl.use('MacOSX')
except:
    a=1
















########################### Example - Tile randomized sources with several Athena params **ON ClUSTER** ###########################

# Using MPI or not
if MPI is not None:
    MPI_size = MPI.COMM_WORLD.Get_size()
    MPI_rank = MPI.COMM_WORLD.Get_rank()
    comm = MPI.COMM_WORLD
else:
    MPI_size=1
    MPI_rank=0

# See if we need to make souce files
USETRACK=False
MAKE_SOURCES = True
if MAKE_SOURCES==True:
    # Split up sources across nodes
    if MPI_size>1:
        source_infer_datafiles=glob.glob(SYNEX_PATH + "/inference_data/Randomized_angles_spins_MRat_*.h5")
        Nvals=len(source_infer_datafiles)
        nPerCore=int(Nvals//MPI_size)
        if MPI_rank==0:
            data_ii_start=0
            data_ii_end=nPerCore if nPerCore<Nvals else Nvals
            source_infer_data = source_infer_datafiles[data_ii_start:data_ii_end]
            for core_ii in range(1,MPI_size):
                data_ii_start=data_ii_end
                data_ii_end=nPerCore*(core_ii+1) if nPerCore*(core_ii+1)<Nvals else Nvals
                if core_ii==MPI_size-1: data_ii_end=Nvals ## Catch any extras left from rounding errors
                comm.send(source_infer_datafiles[data_ii_start:data_ii_end], dest=core_ii)
        else:
            source_infer_data = comm.recv(source=0)
    else:
        source_infer_data = glob.glob(SYNEX_PATH + "/inference_data/Randomized_angles_spins_MRat_*.h5")
        Nvals=len(source_infer_datafiles)

    # Turn data into savefiles
    USETRACK=False
    if MPI_rank==0: print("Making source savefiles from inference data...")
    source_infer_data=glob.glob(SYNEX_PATH + "/inference_data/Randomized_angles_spins_MRat_*.h5") # We are no longer saving raw files so this is sufficient when searching for data files
    for FileName in source_infer_data:
        ExName=SYNEX_PATH + "/Saved_Source_Dicts/" + FileName.split("/inference_data/")[-1].split(".h5")[0]
        Merger = SYU.GetSourceFromLisabetaData(FileName,**{"ExistentialFileName":ExName})
    if MPI_rank==0: print("Done.")

# Catch everything up
comm.Barrier()

# Set verbosity
verbose = False # Verbosity inside SYNEX (making objects etc)
verbose2 = True # Verbosity in this script alone (troubleshooting)

# Telescope args
t0 = '2034-01-01T00:00:00.00' # YYYY-MM-DDTHH:mm:SS.MS 01/01/2034
t = Time(t0, format='isot', scale='utc').gps
Athena_kwargs={"ExistentialFileName":SYNEX_PATH+"/Saved_Telescope_Dicts/Athena_base.dat",
                "verbose":verbose,
                "telescope":"Athena",
                "Tobs":np.array([0.,9.]), # pairs of [Tstart,Tend], for times in DAYS.
                "tilesType" : "moc", # "greedy", # "hierarchical", # "ranked", # moc/greedy/hierarchical/ranked/galaxy.
                "timeallocationType" : "powerlaw", # "absmag" / "powerlaw" / "waw" / "manual" / "pem"
                "scheduleType" : "greedy",
                "doCalcTiles" : False, # Calculates the no. of "hierarchical" tiles based on sky area prob threashold and FoV
                "Ntiles" : None, # 50, # Speific to tilesType=hierarchical and greedy. Needs to be set if "doCalcTiles"=False
                "frozenAthena" : False, # False,
                "exposuretime" : 10000.,
                "min_observability_duration" : None, # 10./3600., # in HOURS
                "inc" : 60., # 60., # In DEGREES, incline of orbital plane normal to Sun-Earth axis.
                "MeanRadius" : 750000000., # 750000000., # meters (from earth-orbit normal axis)
                "semi_maj" : 750000000., # 750000000., # equivalent to MeanRadius axis ONLY IF we say we are really orbiting the centre and not the focal point
                "eccentricity" : 0.4, # 0.8
                "ArgPeriapsis" : 20., # 0., # In DEGREES, angle of point of closest approach to FOCAL POINT IN ORBIT PLANE
                "AscendingNode" : -10., # 0., # In DEGREES
                "phi_0" : 10., # 0., # in DEGREES, initial phase of Athena when measurments start
                "period" : 90., # 180., # In days, for one complete halo orbit about L2
                "gps_science_start" : t, # 1703721618.0, # 01/01/2034 00:00:00.000 UTC -- gps start time of science meaasurements
                "mission_duration" : 2., # In years
                "filt_change_time" : 0., # In seconds?
                "overhead_per_exposure" : 0., # 0., #  In seconds?
                "latitude" : 0., # None, # 20.7204,       ### None if we want a telesscopic orbit?
                "longitude" : 0., # None, # -156.1552,    ### None if we want a telesscopic orbit?
                "elevation" : 0., # None, # 3055.0,       ### None if we want a telesscopic orbit? GWEMOPT uses these for airmass calcs... Ask to raise flag for this?
                "slew_rate" : 1., # None, # 1., # in deg/s -- Idea paper has 1 deg/sec
                "horizon" : 0., # 30.,                    ### None if we want a telesscopic orbit?
                "doMinimalTiling" : True, #  True,
                "readout" : 0.0001, # 0.0001, #
                "doSingleExposure" : True, # False
                "iterativeOverlap" : 0., # 1.0, # Maybe set this to 0? I can't find where this is used...
                "maximumOverlap" : 1.0,
                "sat_sun_restriction" : 5., # 45.,
                "sat_earth_constraint" : 5., # 30.,
                "sat_moon_constraint" : 5., # 20.0,
               } ### What about doPerturbativeTiling? ### "doPerturbativeTiling" : True

# Get sources to test
CutsToTest = ["5hr","10hr","1d","3d","1wk","2wk","3wk","1mon"]
SourceExNames = sorted([File for c in CutsToTest for File in glob.glob(SYNEX_PATH + "/Saved_Source_Dicts/Randomized_angles_spins_MRat_*_"+c+".dat")]) ### Makes "10" go in front of "1"... Problematic af.

# Create Athena Detectors by cloning a base detector
T_obs_array = [np.array([0.,1.]),np.array([0.,2.]),np.array([0.,3.]),np.array([0.,4.])]
cloning_params={"Tobs":T_obs_array}

# Use tracking file or clone new stuff?
CloningTrackFile=SYNEX_PATH+"/TileTrackFiles/ProgressTrackFile.txt" if USETRACK else None

# Tile all combinations of cloning params and sources
SaveInSubFile=None
SaveFileCommonStart="Randomized_angles_spins_MRat"
SourceIndexingByString="MRat_"
T0 = time.time()
SYU.TileSkyArea(CloningTrackFile=CloningTrackFile,sources=SourceExNames,detectors=None,base_telescope_params=Athena_kwargs,cloning_params=cloning_params,SaveInSubFile=SaveInSubFile,SaveFileCommonStart=SaveFileCommonStart,SourceIndexingByString=SourceIndexingByString,verbose=False)
T1 = time.time()
print("Time to finish:",T1-T0,"s")

# Reorder things a bit and make labels
# detectors_out = [[detector_out for detector_out,SysCut in zip(detectors_out,Systems_ID_Cut) if SysCut[1]==Cut] for Cut in CutsToTest]

# Plot
# SYU.PlotPhotonAccumulation(detectors_out, SaveFig=False, SaveFileName=None)
# SYU.PlotSourcePhotons(detectors_out, labels=CutsToTest, BoxPlot=True, SaveFig=False, SaveFileName=None)
