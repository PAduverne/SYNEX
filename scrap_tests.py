import numpy as np
import os
import lisabeta.utils.plotutils as plotutils
import time
import glob
import copy
from astropy.cosmology import Planck13, z_at_value # needed only to convert a given distance to redshift at initialization
from astropy.cosmology import WMAP9 as cosmo
from astropy.time import Time
import json
import healpy as hp
import gwemopt
from gwemopt import utils as gou
import SYNEX.SYNEX_Utils as SYU
import SYNEX.SYNEX_Detectors as SYDs
import SYNEX.SYNEX_Sources as SYSs
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pylab as pylab
from SYNEX.SYNEX_Utils import SYNEX_PATH
from SYNEX.SYNEX_Utils import pylab_params
pylab.rcParams.update(pylab_params)








########################### Example - Plots for Athena 2022 Conf ###########################

# Try all systems with 1d Tcut and any Tobs
# AthenaSaves=glob.glob(SYNEX_PATH+"/Saved_Telescope_Dicts/Tobs_Tests/Randomized_SYNEX2_System_*_1d_TileTime_*.h5")
# AthenaSaves=[SYNEX_PATH+"/Saved_Telescope_Dicts/Tobs_Tests/Randomized_SYNEX2_System_*_1mon_TileTime_4d.h5"]
AthenaSaves=glob.glob(SYNEX_PATH+"/Saved_Telescope_Dicts/Tobs_Tests/Randomized_SYNEX2_System_*_1mon_TileTime_4d.h5")
Athenas=[SYDs.Athena(**{"ExistentialFileName":f,"verbose":False}) for f in AthenaSaves]
Athenas=[Athena for Athena in Athenas if hasattr(Athena,"detector_source_coverage") and "source save file" in Athena.detector_source_coverage]
# Athenas=[Athena for Athena in Athenas if "cloning value" in Athena.detector_source_coverage]

# Get source savefiles
Mergers = [SYSs.SMBH_Merger(**{"ExistentialFileName":Athena.detector_source_coverage["source save file"],"verbose":False}) for Athena in Athenas]

# Retile one thing to check some new stuff
SaveInSubFile="Tobs_Tests"
SaveFileCommonStart="Randomized_angles_spins_MRat"
SourceIndexingByString="MRat_"
print("Pre-tiling check:",[Athenas[0].detector_source_coverage["Source tile timeranges (isot)"][:] if "Source tile timeranges (isot)" in Athenas[0].detector_source_coverage else None])
setattr(Athenas[0],"detector_source_coverage", None)
Athena_out=SYU.TileSkyArea(sources=Mergers[0],detectors=Athenas[0],base_telescope_params=None,cloning_params=None,SaveInSubFile=SaveInSubFile,SaveFileCommonStart=SaveFileCommonStart,SourceIndexingByString=SourceIndexingByString,verbose=False)
Athenas[0]=Athena_out[0]
print("Post-tiling check:",Athenas[0].detector_source_coverage["Source tile timeranges (isot)"])
print("Another check:",Athenas[0].detector_config_struct["cloning keys"], Athenas[0].detector_config_struct["cloning values"])
# Check things are as we expect...
# print("Load checks 1:",[-Merger.DeltatL_cut/(24*60*60) for Merger in Mergers], [Athena.detector_go_params["Tobs"] for Athena in Athenas])
# print("Load checks 2:",len([-Merger.DeltatL_cut/(24*60*60) for Merger in Mergers]), len([Athena.detector_go_params["Tobs"] for Athena in Athenas])) # , [Merger.ExistentialFileName for Merger in Mergers], [Athena.ExistentialFileName for Athena in Athenas])
# print("Load checks 3:",[Athena.detector_source_coverage.keys() for Athena in Athenas]) # , [Merger.ExistentialFileName for Merger in Mergers], [Athena.ExistentialFileName for Athena in Athenas])
print("Load checks 4:",Athenas[0].detector_source_coverage["Source tile timeranges (isot)"],Athenas[0].detector_source_coverage["Source tile timeranges (days)"])

# Plot an example Skymap
SYU.PlotSkyMapData(Mergers[0],SaveFig=False,plotName=None,DO_CONTOURS=False)

# Plot example orbit
# SYU.PlotOrbit(Athenas[0].detector_config_struct,MollProj=False,SaveFig=False)

# Plot source photons
SYU.PlotSourcePhotons(Athenas, labels=None, BoxPlot=True, SaveFig=False, SaveFileName=None)

# Run through some basic correlation analyses
Data = SYU.DecTreeReg(Athenas)


























# ########################### Example - Analyze some randomized results ###########################
#
# # Try all systems with 1d Tcut and any Tobs
# AthenaSaves=glob.glob(SYNEX_PATH+"/Saved_Telescope_Dicts/Tobs_Tests/Randomized_SYNEX2_System_*_1d_TileTime_*.h5")
# Athenas=[SYDs.Athena(**{"ExistentialFileName":f,"verbose":False}) for f in AthenaSaves]
#
# # Run through some basic correlation analyses
# Data = SYU.DecTreeReg(Athenas)










































# ########################### By-hand calculation of expected source photons ###########################
#
# # Source
# FileName = "IdeaPaperSystem_9d"
# Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/TestSystem_9d_base.dat"}
# Merger = SYU.GetSourceFromLisabetaData(FileName,**Merger_kwargs)
#
# # Check what's inside the ARF file
# from astropy.io import fits
# from SYNEX.SYNEX_Utils import SYNEX_PATH
# ARF_file=SYNEX_PATH+"/XIFU_CC_BASELINECONF_2018_10_10.arf"
# hdul = fits.open(ARF_file)
# # hdul.info()
# ARF1 = hdul[1].data[:]
# N = len(ARF1)
#
# # Calculate EM flux and CTR
# Merger.GenerateEMFlux()
# Merger.GenerateCTR(ARF_file)







































# ########################### Example - Load saved objects after tiling and plot useful results ###########################
#
# # Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_dev.dat"}
# # Merger=SYSs.SMBH_Merger(**Merger_kwargs)
#
# DetSaveFiles1=sorted(glob.glob("/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/ExpTime_doSingleExposure/Athena_base_exposuretime_*.dat"))
# DetSaveFiles2=sorted(glob.glob("/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/ExpTime_Not_doSingleExposure/Athena_base_exposuretime_*.dat"))
# DetSaveFiles3=sorted(glob.glob("/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/ExpTime_Overhead_doSingleExposure/Athena_base_exposuretime_*.dat"))
# DetSaveFiles4=sorted(glob.glob("/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/ExpTime_Overhead_Not_doSingleExposure/Athena_base_exposuretime_*.dat"))
# detectors1=[SYDs.Athena(**{"ExistentialFileName":DetSaveFile}) for DetSaveFile in DetSaveFiles1]
# detectors2=[SYDs.Athena(**{"ExistentialFileName":DetSaveFile}) for DetSaveFile in DetSaveFiles2]
# detectors3=[SYDs.Athena(**{"ExistentialFileName":DetSaveFile}) for DetSaveFile in DetSaveFiles3]
# detectors4=[SYDs.Athena(**{"ExistentialFileName":DetSaveFile}) for DetSaveFile in DetSaveFiles4]
# detectors=[detectors1,detectors2,detectors3,detectors4]
#
# labels=[r"No overhead, doSingleExposure=True",r"No overhead, doSingleExposure=False",
#         r"Overhead, doSingleExposure=True",r"Overhead, doSingleExposure=False"]
#
# # SYU.PlotPhotonAccumulation(detectors, SaveFig=False, SaveFileName=None)
# SYU.PlotSourcePhotons(detectors, labels=labels, SaveFig=False, SaveFileName=None)


























#########################################################################################################################
