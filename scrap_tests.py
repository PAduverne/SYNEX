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





########################### Example - Tile sources with several Tobs ###########################

# Set verbosity
verbose = False # Verbosity inside SYNEX (making objects etc)
verbose2 = True # Verbosity in this script alone

# Telescope args
t0 = '2034-01-01T00:00:00.00' # YYYY-MM-DDTHH:mm:SS.MS 01/01/2034
t = Time(t0, format='isot', scale='utc').gps
Athena_kwargs={"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_base.dat",
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
Athena_kwargs_WithNewEx = copy.deepcopy(Athena_kwargs)
Athena_kwargs_WithNewEx.update({"NewExistentialFileName":None})

# Get sources to test
CutsToTest = ["1d","1wk","3wk"]
FileNames = [File for c in CutsToTest for File in glob.glob(SYNEX_PATH + "/inference_data/Randomized_SYNEX2/Randomized_angles_spins_MRat_*_"+c+".h5")]
Systems_ID_Cut = [FileName[FileName.rfind('t')+2:FileName.rfind('.')] for FileName in FileNames] ### 't' is last letter before ID that does not occur after ID.
Systems_ID_Cut=[tuple(System_ID_Cut.split("_")) for System_ID_Cut in Systems_ID_Cut]
SystemIDs = [int(System_ID_Cut[0]) for System_ID_Cut in Systems_ID_Cut]
ExNames = [FileName.replace("/inference_data/", "/Saved_Source_Dicts/").replace(".h5", ".dat") for FileName in FileNames]
Mergers = [SYU.GetSourceFromLisabetaData(FileName,**{"ExistentialFileName":ExName,"verbose":verbose}) for FileName,ExName in zip(FileNames,ExNames)]

# Create Athena Detectors
T_obs_array = [np.array([0.,1.]),np.array([0.,2.]),np.array([0.,3.]),np.array([0.,4.])] ###
n_T_obs = len(T_obs_array)
Tobs_maxes = [T[1] for T in T_obs_array]
ExNames=[]
for System_ID_Cut in Systems_ID_Cut:
    ExNames.extend([SYNEX_PATH+"/Saved_Telescope_Dicts/Tobs_Tests/Randomized_SYNEX2_System_"+"_".join(System_ID_Cut)+"_TileTime_"+str(int(T_obs_array[i][1]))+"d.h5" for i in range(len(Tobs_maxes))])
Athena_kwargs_dicts=[dict(Athena_kwargs, **{"ExistentialFileName":ExName,"Tobs":np.array([0.,int(ExName.split("_TileTime_")[-1].split("d.")[0])])}) if os.path.isfile(ExName) else dict(Athena_kwargs_WithNewEx, **{"NewExistentialFileName":ExName, "Tobs":np.array([0.,int(ExName.split("_TileTime_")[-1].split("d.")[0])])}) for ExName in ExNames]
Athenas = [SYDs.Athena(**Athena_kwargs_dict) for Athena_kwargs_dict in Athena_kwargs_dicts]

print("-"*50)
for i in range(3):
    print(Mergers[i].ExistentialFileName)
    print(Athenas[n_T_obs*i].ExistentialFileName)
    print(Athenas[n_T_obs*i+1].ExistentialFileName)
    print(Athenas[n_T_obs*i+2].ExistentialFileName)
    print(Athenas[n_T_obs*i+3].ExistentialFileName)
    print("-"*50)

# Tile all telescope classes that haven't already been tiled -- skip those that have
tiling_t0=time.time()
detectors_out=[]
for iMerg,Merger in enumerate(Mergers):
    for Athena in Athenas[n_T_obs*iMerg:n_T_obs*(iMerg+1)]:
        if Athena.detector_source_coverage!=None:
            if verbose2: print("-"*20,"Detector already tiled","-"*20)
            print(Athena.detector_go_params["Tobs"], -Merger.DeltatL_cut/(24.*60.*60.))
            print(Athena.ExistentialFileName, Merger.ExistentialFileName, Merger.chi1)
            detectors_out += [Athena]
        else:
            if verbose2: print("-"*20,"Detector not yet tiled","-"*20)
            tiling_t0_a = time.time()
            print(Athena.detector_go_params["Tobs"], -Merger.DeltatL_cut/(24.*60.*60.))
            print(Athena.ExistentialFileName, Merger.ExistentialFileName, Merger.chi1)
            detectors_out += SYU.TileSkyArea(Merger,detectors=Athena,base_telescope_params=None,cloning_params=None,verbose=verbose)
            tiling_t0_b = time.time()
            if verbose2: print("Merger system:", iMerg+1, "/", len(Mergers), " in", tiling_t0_b-tiling_t0_a, "s")
        if verbose2: print("\n")
        if verbose2: print("-"*63)
tiling_t1=time.time()
if verbose2: print("Total time for",len(detectors_out),"sources:",tiling_t1-tiling_t0,"s")

# Reorder things a bit and make labels
detectors_out = [[detector_out for detector_out,SysCut in zip(detectors_out,Systems_ID_Cut) if SysCut[1]==Cut] for Cut in CutsToTest]

# Plot
# SYU.PlotPhotonAccumulation(detectors_out, SaveFig=False, SaveFileName=None)
SYU.PlotSourcePhotons(detectors_out, labels=CutsToTest, BoxPlot=True, SaveFig=False, SaveFileName=None)

































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
