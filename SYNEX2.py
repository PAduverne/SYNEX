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

from numpy.random import rand

import time
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






# ########################### Example - Create and save source object from lisabeta file ###########################
#
# # Merger args
# Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_base.dat"}
# FileName="IdeaPaperSystem_9d"
# Merger = SYU.GetSourceFromLisabetaData(FileName,**Merger_kwargs)
#
# MergerDict=dict(Merger.__dict__)
# for key,val in MergerDict.items(): print(key,val)










########################### Example - Tile using gwemopt on Cluster ###########################

# Merger args
Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_base.dat",
                 "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_dev.dat"}

# Base telescope args
t0 = '2034-01-01T00:00:00.00' # YYYY-MM-DDTHH:mm:SS.MS 01/01/2034
t = Time(t0, format='isot', scale='utc').gps
Athena_kwargs={"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_base.dat",
                "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_dev.dat",
                "orbitFile":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/orbit_files/Athena_20340601_728d_inc60_R750Mkm_ecc4_ArgPeri20_AscNode-10_phi020_P90_frozenFalse_base.dat",
                "NeworbitFile":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/orbit_files/Athena_20340601_728d_inc60_R750Mkm_ecc4_ArgPeri20_AscNode-10_phi020_P90_frozenFalse_dev.dat",
                "telescope":"Athena",
                "tilesType" : "moc", # "greedy", # "hierarchical", # "ranked", # moc/greedy/hierarchical/ranked/galaxy.
                "timeallocationType" : "powerlaw", # "absmag" / "powerlaw" / "waw" / "manual" / "pem"
                "scheduleType" : "greedy",
                "doCalcTiles" : False, # Calculates the no. of "hierarchical" tiles based on sky area prob threashold and FoV
                "Ntiles" : None, # 50, # Speific to tilesType=hierarchical and greedy. Needs to be set if "doCalcTiles"=False
                "frozenAthena" : False, # False,
                "exposuretime" : None, # 10000., # 6*60*60, # 60., #
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
                "overhead_per_exposure" : 0., # In seconds?
                "latitude" : 0., # None, # 20.7204,       ### None if we want a telesscopic orbit?
                "longitude" : 0., # None, # -156.1552,    ### None if we want a telesscopic orbit?
                "elevation" : 0., # None, # 3055.0,       ### None if we want a telesscopic orbit? GWEMOPT uses these for airmass calcs... Ask to raise flag for this?
                "slew_rate" : None, # 1., # in deg/s -- Idea paper has 1 deg/sec
                "horizon" : 0., # 30.,                    ### None if we want a telesscopic orbit?
                "doMinimalTiling" : True, #  True,
                "readout" : None, # 0.0001, # 6
                "doSingleExposure" : False, # False
                "iterativeOverlap" : 1.0,
                "maximumOverlap" : 1.0,
                "sat_sun_restriction" : 5., # 45.,
                "sat_earth_constraint" : 5., # 30.,
                "sat_moon_constraint" : 5., # 20.0,
               } ### What about doPerturbativeTiling? ### "doPerturbativeTiling" : True


# Check what's inside the ARF file
from astropy.io import fits
from SYNEX.SYNEX_Utils import SYNEX_PATH
ARF_file=SYNEX_PATH+"/XIFU_CC_BASELINECONF_2018_10_10.arf"
hdul = fits.open(ARF_file)
# hdul.info()
print("ARF file contents:",type(hdul))
ARF = hdul[1].data[:]
N = len(ARF)
print("ARF file contents:",type(ARF),np.shape(ARF),type(ARF[0]),np.shape(ARF[0]),ARF[0:3])

# Test tiling with detector cloning
ex_times=np.logspace(2,4,num=20,endpoint=True,base=10.)
cloning_params={"exposuretime":ex_times}
tiling_t0=time.time()
detectors = SYU.TileSkyArea(Merger_kwargs,detectors=None,base_telescope_params=Athena_kwargs,cloning_params=cloning_params)
tiling_t1=time.time()
print("Total time for",len(detectors),"detectors:",tiling_t1-tiling_t0, "s")
















# ########################### Example - Init source and telescope and tile using gwemopt ###########################
# # Source file we want to resurract - note we only specify everything afer .../inference_data or .../inference_param_files and don't need to give suffix either
# # FileName = "Randomized_SYNEX/RemakePaperPlots/Randomized_angles_spins_MRat_10_1wk" # "IdeaPaperSystem_9d"
#
# # Merger args
# Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_base.dat",
#                  "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_dev.dat"}
# # Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/TestSystem_9d_base.dat"}
# # Merger_kwargs = {"ExistentialFileName":"/home/baird/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_base.dat",
# #                  "NewExistentialFileName":"/home/baird/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_dev.dat"}
#
# # Resurrect - either from lisabeta data or saved source file
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
# # Merger = SYU.GetSourceFromLisabetaData(FileName,**Merger_kwargs)
#
# # Make some test telescopes
# ###### TO DO -- check when calculating orbit that it is intersecting properly between orbit start and obs times...
# t0 = '2034-01-01T00:00:00.00' # YYYY-MM-DDTHH:mm:SS.MS 01/01/2034
# t = Time(t0, format='isot', scale='utc').gps
# # Athena_kwargs={"FOV":1.,"exposuretime":60.,"slew_rate":1., "telescope":"Athena_base",
# #             "ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_base.dat"}
# Athena_kwargs={"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_base.dat",
#                 "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_dev.dat",
#                 "orbitFile":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/orbit_files/Athena_20340601_728d_inc60_R750Mkm_ecc4_ArgPeri20_AscNode-10_phi020_P90_frozenFalse_base.dat",
#                 "NeworbitFile":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/orbit_files/Athena_20340601_728d_inc60_R750Mkm_ecc4_ArgPeri20_AscNode-10_phi020_P90_frozenFalse_dev.dat",
#                 "telescope":"Athena",
#                 "tilesType" : "moc", # "greedy", # "hierarchical", # "ranked", # moc/greedy/hierarchical/ranked/galaxy.
#                 "timeallocationType" : "powerlaw", # "absmag" / "powerlaw" / "waw" / "manual" / "pem"
#                 "scheduleType" : "greedy",
#                 "doCalcTiles" : False, # Calculates the no. of "hierarchical" tiles based on sky area prob threashold and FoV
#                 "Ntiles" : None, # 50, # Speific to tilesType=hierarchical and greedy. Needs to be set if "doCalcTiles"=False
#                 "frozenAthena" : False, # False,
#                 "exposuretime" : 10000., # 6*60*60, # 60., #
#                 "min_observability_duration" : 10./3600., # in HOURS
#                 "inc" : 60., # 60., # In DEGREES, incline of orbital plane normal to Sun-Earth axis.
#                 "MeanRadius" : 750000000., # 750000000., # meters (from earth-orbit normal axis)
#                 "semi_maj" : 750000000., # 750000000., # equivalent to MeanRadius axis ONLY IF we say we are really orbiting the centre and not the focal point
#                 "eccentricity" : 0.4, # 0.8
#                 "ArgPeriapsis" : 20., # 0., # In DEGREES, angle of point of closest approach to FOCAL POINT IN ORBIT PLANE
#                 "AscendingNode" : -10., # 0., # In DEGREES
#                 "phi_0" : 10., # 0., # in DEGREES, initial phase of Athena when measurments start
#                 "period" : 90., # 180., # In days, for one complete halo orbit about L2
#                 "gps_science_start" : t, # 1703721618.0, # 01/01/2034 00:00:00.000 UTC -- gps start time of science meaasurements
#                 "mission_duration" : 2., # In years
#                 "filt_change_time" : 0., # In seconds?
#                 "overhead_per_exposure" : 0., # In seconds?
#                 "latitude" : 0., # None, # 20.7204,       ### None if we want a telesscopic orbit?
#                 "longitude" : 0., # None, # -156.1552,    ### None if we want a telesscopic orbit?
#                 "elevation" : 0., # None, # 3055.0,       ### None if we want a telesscopic orbit? GWEMOPT uses these for airmass calcs... Ask to raise flag for this?
#                 "slew_rate" : 100., # 1., # in deg/s -- Idea paper has 1 deg/sec
#                 "horizon" : 0., # 30.,                    ### None if we want a telesscopic orbit?
#                 "doMinimalTiling" : True, #  True,
#                 "readout" : 0.0001, # 6
#                 "doSingleExposure" : True, # False
#                 "iterativeOverlap" : 1.0,
#                 "maximumOverlap" : 1.0,
#                } ### What about doPerturbativeTiling? ### "doPerturbativeTiling" : True
# # Athena_kwargs["ExistentialFileName"]="/home/baird/SYNEX/Saved_Telescope_Dicts/Athena_base.dat"
# # Athena_kwargs["NewExistentialFileName"]="/home/baird/SYNEX/Saved_Telescope_Dicts/Athena_dev.dat"
# # Athena_kwargs["orbitFile"]="/home/baird/SYNEX/orbit_files/Athena_20340601_728d_inc60_R750Mkm_ecc4_ArgPeri20_AscNode-10_phi020_P90_frozenFalse_base.dat"
# # Athena_kwargs["NeworbitFile"]="/home/baird/SYNEX/orbit_files/Athena_20340601_728d_inc60_R750Mkm_ecc4_ArgPeri20_AscNode-10_phi020_P90_frozenFalse_dev.dat"
# Athena=SYDs.Athena(**Athena_kwargs)
#
# # Test tiling with detector cloning
# ex_times=np.logspace(1,4,num=2,endpoint=True,base=10.)
# cloning_params={"exposuretime":ex_times} # {"inc":np.linspace(0., 90., 5)}
# tiling_t0=time.time()
# detectors = SYU.TileSkyArea(Merger,detectors=Athena,base_telescope_params=None,cloning_params=cloning_params)
# # detectors = SYU.TileSkyArea(Merger,detectors=None,base_telescope_params=Athena_kwargs,cloning_params=cloning_params)
# tiling_t1=time.time()
# print("Total time for",len(detectors),"detectors:",tiling_t1-tiling_t0, "s") # 1216 for 5 detectors...
#
# for detector in detectors:
#     T0_mjd=detector.detector_source_coverage["Start time (mjd)"]
#     Xs,Ys=[],[]
#     for ExT0s,ExTs in zip(detector.detector_source_coverage["Source tile start times (s)"],detector.detector_source_coverage["Source tile exposuretimes (s)"]): Xs+=[ExT0s/86400.,(ExT0s+ExTs)/86400.]
#     for CumCounts in np.cumsum(detector.detector_source_coverage["Source photon counts"]): Ys+=[CumCounts,CumCounts]
#     print(detector,Xs,Ys,detector.detector_source_coverage["Source tile exposuretimes (s)"])
#     plt.plot(Xs,Ys,label=detector.detector_config_struct["telescope"])
#
# plt.xlabel(r"Time from "+str(T0_mjd)+" (mjd)")
# plt.ylabel(r"Source photons")
# plt.legend()
# plt.grid()
# plt.show()















########################### Example - Init source and telescope and tile using gwemopt ###########################

# # Source file we want to resurract - note we only specify everything afer .../inference_data or .../inference_param_files and don't need to give suffix either
# # FileName = "Randomized_SYNEX/RemakePaperPlots/Randomized_angles_spins_MRat_10_1wk" # "IdeaPaperSystem_9d"
#
# # Merger args
# Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_base.dat",
#                  "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_dev.dat",}
# # Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/TestSystem_9d_base.dat"}
#
# # Resurrect - either from lisabeta data or saved source file
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
# # Merger = SYU.GetSourceFromLisabetaData(FileName,**Merger_kwargs)
#
# # # Plot the CL bounds
# # SYU.PlotInferenceLambdaBeta(Merger.H5File, bins=50, SkyProjection=True, SaveFig=False, return_data=False)
# #
# # # Plot
# # SYU.PlotSkyMapData(Merger,SaveFig=False,plotName=None)
#
# # Make some test telescopes
# ###### TO DO -- check when calculating orbit that it is intersecting properly between orbit start and obs times...
# t0 = '2034-01-01T00:00:00.00' # YYYY-MM-DDTHH:mm:SS.MS 01/01/2034
# t = Time(t0, format='isot', scale='utc').gps
# # Athena_kwargs={"FOV":1.,"exposuretime":60.,"slew_rate":1., "telescope":"Athena_base",
# #             "ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_base.dat"}
# Athena_kwargs={"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_base.dat",
#                 "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_dev_list.dat",
#                 "orbitFile":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/orbit_files/Athena_20340601_728d_inc60_R750Mkm_ecc4_ArgPeri20_AscNode-10_phi020_P90_frozenFalse_base.dat",
#                 "NeworbitFile":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/orbit_files/Athena_20340601_728d_inc60_R750Mkm_ecc4_ArgPeri20_AscNode-10_phi020_P90_frozenFalse_dev.dat",
#                 "telescope":"Athena",
#                 "tilesType" : "moc", # "greedy", # "hierarchical", # "ranked", # moc/greedy/hierarchical/ranked/galaxy.
#                 "timeallocationType" : "powerlaw", # "absmag" / "powerlaw" / "waw" / "manual" / "pem"
#                 "scheduleType" : "greedy",
#                 "doCalcTiles" : False, # Calculates the no. of "hierarchical" tiles based on sky area prob threashold and FoV
#                 "Ntiles" : None, # 50, # Speific to tilesType=hierarchical and greedy. Needs to be set if "doCalcTiles"=False
#                 "frozenAthena" : False, # False,
#                 "exposuretime" : 10000., # 6*60*60, # 60., #
#                 "min_observability_duration" : 10./3600., # in HOURS
#                 "inc" : 60., # 60., # In DEGREES, incline of orbital plane normal to Sun-Earth axis.
#                 "MeanRadius" : 750000000., # 750000000., # meters (from earth-orbit normal axis)
#                 "semi_maj" : 750000000., # 750000000., # equivalent to MeanRadius axis ONLY IF we say we are really orbiting the centre and not the focal point
#                 "eccentricity" : 0.4, # 0.8
#                 "ArgPeriapsis" : 20., # 0., # In DEGREES, angle of point of closest approach to FOCAL POINT IN ORBIT PLANE
#                 "AscendingNode" : -10., # 0., # In DEGREES
#                 "phi_0" : 10., # 0., # in DEGREES, initial phase of Athena when measurments start
#                 "period" : 90., # 180., # In days, for one complete halo orbit about L2
#                 "gps_science_start" : t, # 1703721618.0, # 01/01/2034 00:00:00.000 UTC -- gps start time of science meaasurements
#                 "mission_duration" : 2., # In years
#                 "filt_change_time" : 0., # In seconds?
#                 "overhead_per_exposure" : 0., # In seconds?
#                 "latitude" : 0., # None, # 20.7204,       ### None if we want a telesscopic orbit?
#                 "longitude" : 0., # None, # -156.1552,    ### None if we want a telesscopic orbit?
#                 "elevation" : 0., # None, # 3055.0,       ### None if we want a telesscopic orbit? GWEMOPT uses these for airmass calcs... Ask to raise flag for this?
#                 "slew_rate" : 100., # 1., # in deg/s -- Idea paper has 1 deg/sec
#                 "horizon" : 0., # 30.,                    ### None if we want a telesscopic orbit?
#                 "doMinimalTiling" : True, #  True,
#                 "readout" : 0.0001, # 6
#                 "doSingleExposure" : False, # False
#                 "iterativeOverlap" : 1.0,
#                 "maximumOverlap" : 1.0,
#                } ### What about doPerturbativeTiling? ### "doPerturbativeTiling" : True
# Athena=SYDs.Athena(**Athena_kwargs)
#
# # # Check orbit
# import SYNEX.segments_athena as segs_a
# # config_struct = segs_a.calc_telescope_orbit(Athena.detector_config_struct,SAVETOFILE=False)
# # SYU.PlotOrbit(config_struct,MollProj=False,SaveFig=False)
# # SYU.AnimateOrbit(config_struct,include_sun=False,SaveAnim=False)
# # SYU.AnimateSkyProjOrbit(Athena.detector_config_struct,MollProj=True,SaveAnim=True)
#
#
#
# # Test tiling with detector cloning
# ex_times=np.logspace(1,4,num=5,endpoint=True,base=10.)
# cloning_params={"exposuretime":ex_times} # {"inc":np.linspace(0., 90., 5)}
# t0=time.time()
# go_params, map_struct, tile_structs, coverage_struct, DetectorInfo, detectors = SYU.TileSkyArea(Merger,detectors=Athena,base_telescope_params=None,cloning_params=cloning_params)
# print("Final check on detector:", [d.detector_config_struct["exposuretime"] for d in detectors])
# t1=time.time()
# print("Time for list of detectors:",t1-t0) # About 410 for 15 detectors
#
# for telescope in tile_structs.keys():
#     tile_struct=tile_structs[telescope]
#     TelTileProperties=np.array([(telescope, tile_struct[i]["ra"], tile_struct[i]["dec"], tile_struct[i]["exposureTime"], tile_struct[i]["nexposures"]) for i in tile_struct.keys()])
#
# # # Data is N*9 with columns ra,dec,mjd,mag,exposureTime,field,tile_probs,airmass,program_id
# # TelCovProperties=[list(zip(list(d.detector_coverage_struct["telescope"]), list(d.detector_coverage_struct["data"][:,0]), list(d.detector_coverage_struct["data"][:,1]), list(d.detector_coverage_struct["data"][:,2]), list(d.detector_coverage_struct["data"][:,4]))) for d in detectors]
# #
# # f = open('List_vs_Obo_coverage_info.txt', 'w')
# # f.write("Detector list\n")
# # for d in detectors:
# #     TelProp = list(zip(list(d.detector_coverage_struct["telescope"]), list(d.detector_coverage_struct["data"][:,0]), list(d.detector_coverage_struct["data"][:,1]), list(d.detector_coverage_struct["data"][:,2]), list(d.detector_coverage_struct["data"][:,4])))
# #     print("TelProp checks:",np.shape(TelProp),type(TelProp))
# #     for i in range(np.shape(TelProp)[0]): f.write(" ".join([str(el) for el in TelProp[i][:]])+"\n")
# #
# #
# # f.write("\n\n")
# # f.write("One by one")
# # f.write("\n")
#
#
#
# # Test input detectors one by one
# Athena_kwargs["NewExistentialFileName"]="/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_dev_ObO.dat"
# DetectorInfo2=[]
# t2=time.time()
# for detii,ex_t in enumerate(ex_times):
#     cloning_params={"exposuretime":ex_t}
#     go_params, map_struct, tile_structs, coverage_struct, DetInf, detector = SYU.TileSkyArea(Merger,detectors=Athena,base_telescope_params=None,cloning_params=cloning_params)
#     DetectorInfo2+=list(DetInf)
#     d=detector[0]
#     TelProp = list(zip(list(d.detector_coverage_struct["telescope"]), [detii]*len(d.detector_coverage_struct["telescope"]), list(d.detector_coverage_struct["data"][:,0]), list(d.detector_coverage_struct["data"][:,1]), list(d.detector_coverage_struct["data"][:,2]), list(d.detector_coverage_struct["data"][:,4])))
#     print("TelProp checks (obo):",np.shape(TelProp),type(TelProp))
#     # for i in range(np.shape(TelProp)[0]): f.write(" ".join([str(el) for el in TelProp[i][:]])+"\n")
# DetectorInfo2=np.array(DetectorInfo2)
# t3=time.time()
#
# # f.close()
#
# print("Time for one by one:",t3-t2) # About 1300 for 15 detectors
#
# # Make plots of results
# plt.plot(ex_times,DetectorInfo[:,1],label=r"List of dets")
# plt.plot(ex_times,DetectorInfo2[:,1],label=r"One-by-one dets")
# plt.xlabel(r"Exposure time [s]")
# plt.ylabel(r"Accumulated sky probability")
# ax=plt.gca()
# ax.set_xscale('log')
# plt.legend(fontsize="x-small")
# plt.show()
#
# plt.plot(ex_times,DetectorInfo[:,3],"-",label=r"Source captures (list)")
# plt.plot(ex_times,DetectorInfo[:,4],":",label=r"Tot. coverage (list)")
# plt.plot(ex_times,DetectorInfo[:,5],".-",label=r"Unique coverage (list)")
# plt.plot(ex_times,DetectorInfo2[:,3],"-",label=r"Source captures (Obo)")
# plt.plot(ex_times,DetectorInfo2[:,4],":",label=r"Tot. coverage (Obo)")
# plt.plot(ex_times,DetectorInfo2[:,5],".-",label=r"Unique coverage (Obo)")
# plt.xlabel(r"Exposure time [s]")
# plt.ylabel(r"No. of tiles")
# ax=plt.gca()
# ax.set_xscale("log")
# plt.legend(fontsize="x-small")
# plt.show()
#
# # Print out own summary...
# # SYU.GetCoverageInfo(go_params, map_struct, tile_structs, coverage_struct)
#
# # See some stats and plots plots
# # gwemopt.plotting.tiles(go_params, map_struct, tile_structs)
# # print("Start of summary call: -----------------------------------------------")
# # gwemopt.scheduler.summary(go_params, map_struct, coverage_struct)
# # SYU.PlotSkyMapData(Merger,SaveFig=False,plotName=None)
# # print("Start of plotting call: -----------------------------------------------")
# # gwemopt.plotting.coverage(go_params, map_struct, coverage_struct)



















# ########################### Example - Plotting Athena orbit ###########################
#
# import SYNEX.segments_athena as segs_a
#
# t0 = '2034-06-01T00:00:00.00' # YYYY-MM-DDTHH:mm:SS.MS
# t = Time(t0, format='isot', scale='utc').gps
#
# # Athena_kwargs={"FOV":1.,"exposuretime":60.,"slew_rate":1., "telescope":"Athena_1",
# #             "ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_1_base.dat"}
# Athena_kwargs={
#                 "ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_1_base.dat",
#                 "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_1_dev.dat",
#                 "frozenAthena" : False, # False,
#                 "exposuretime" : 6*60*60,
#                 "inc" : 60., # 60., # In DEGREES, incline of orbital plane normal to Sun-Earth axis.
#                 "MeanRadius" : 750000000, # 750000000., # meters (from earth-orbit normal axis)
#                 "semi_maj" : 750000000., # 750000000., # equivalent to MeanRadius axis ONLY IF we say we are really orbiting the centre and not the focal point
#                 "eccentricity" : 0.4, # 0.8
#                 "ArgPeriapsis" : 20., # 0., # In DEGREES, angle of point of closest approach to FOCAL POINT IN ORBIT PLANE
#                 "AscendingNode" : -10., # 0., # In DEGREES
#                 "phi_0" : 10., # 0., # in DEGREES, initial phase of Athena when measurments start
#                 "period" : 90., # 180., # In days, for one complete halo orbit about L2
#                 "gps_science_start" : t, # 1703721618.0, # 01/01/2034 00:00:00.000 UTC -- gps start time of science meaasurements
#                 "mission_duration" : 2.
#                }
#
# Athena_1=SYDs.Athena(**Athena_kwargs)
#
# config_struct = segs_a.calc_telescope_orbit(Athena_1.detector_config_struct,SAVETOFILE=False)
#
# SYU.PlotOrbit(config_struct, SaveFig=False)
# SYU.AnimateOrbit(config_struct,include_sun=False,SaveAnim=False)











# ########################### Example - GWEMOPT tiling strategies from within SYNEX ###########################
#
# # Initiate Athena
# Athena_kwargs = {"FoView":1.*(np.pi/180.)**2,"T_lat":10000.,"T_init":0.,"slew_rate":None} # slew_rate=None for no gaps for slew
# Athena = SYDs.Athena(**Athena_kwargs)
#
# # Define datafile
# dataFile = "IdeaPaperSystem_9d.h5"
#
# # Start tiling
# Athena.TileSkyArea(dataFile,TileStrat="moc",SAVE_SOURCE_EM_PROPERTIES_IN_TILE_JSON=True,go_params=None)
#
# # Plot tiles using gwemopt directly...
# # gwemopt.plotting.tiles(params, map_struct, tile_structs)
#
# # Plot tiles using SYNEX... this should break right now.
# SYU.PlotTilesArea_old(Athena.TilePickleName,n_tiles=50)
#
# # Tiling file name
# TilePickleName = "/Users/baird/Documents/LabEx_PostDoc/SYNEX/Tile_files/IdeaPaperSystem_9d_tiled.dat"
#
# # Get kuiper detection pvals
# Athena.GetKuiper(Athena.TilePickleName,source=Merger)
# # Athena.GetKuiper(TileJsonName)








# ########################### Example - GWEMOPT tiling strategies ###########################
#
# # Data for sky map
# FileLocAndName = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/IdeaPaperSystem_9d.h5"
#
# # Create the params dictionary and run the built-in gwemopt check function
# go_params={"nside":16}
# go_params = gou.params_checker(go_params)
#
# DO_PIXEL_LOC_PLOT = False
# if DO_PIXEL_LOC_PLOT:
#     npix = hp.nside2npix(go_params["nside"])
#     pix_thetas, pix_phis = hp.pix2ang(go_params["nside"], np.arange(npix))
#     fig = plt.figure()
#     ax = plt.gca()
#     pix_lambdas = pix_phis-np.pi
#     pix_betas = np.pi/2.-pix_thetas
#     plt.subplot(111, projection="aitoff")
#     plt.scatter(pix_lambdas,pix_betas,s=1)
#     plt.show()
#
# # Create the necessary sky_map dictionary
# go_params,map_struct = SYU.CreateSkyMapStruct(go_params,FileName=FileLocAndName)
#
# # Get useful stuff for checks
# pix_ras = map_struct["ra"]
# pix_decs = map_struct["dec"]
# pix_probs = map_struct["prob"]
#
# # Do a plot to check everything is ok so far
# DO_PROB_PLOT = True
# if DO_PROB_PLOT:
#     import matplotlib.tri as tri
#     triang = tri.Triangulation(pix_ras, pix_decs)
#     fig1, ax1 = plt.subplots()
#     fig1.colorbar(tcf)
#     # ax1.tricontour(triang, pix_probs, colors='k')
#     plt.show()





# ########################### Example - Calculate an example EM flux for a source ###########################
#
# # Create options to Initialize the detector object
# LISA_kwargs = {"TDI":'TDIAET'}
#
# # NB Dec = beta and Ra = lambda
# Merger_kwargs = {"q": 10., "M": 5e6, "dist": 13600, "chi1": 0.9,
#         "chi2": 0.95, "beta" : 0., "lambda" : 0.,
#         "inc": 40.3*np.pi/180., "psi": 0.,  "approximant" : 'IMRPhenomHM'} # ,
#         # "DeltatL_cut":T_obs_end_to_merger} # ,
#         # "timetomerger_max":2.5, "minf":1e-6}
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_kwargs)
#
# # Initialize the source object
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
# # Compute the xray flux
# Merger.GenerateEMFlux(LISA) # fstart22)
#
# # Plot the results to see if we get what they have in their paper
# Xs = Merger.xray_time/(60.*60.*24.)
# Ys = Merger.xray_flux
# plt.plot(Xs,Ys)
# ax = plt.gca()
# ax.set_yscale('log')
# plt.ylabel(r'Energy Flux [erg s$^{-1}$ cm$^{-2}$]')
# plt.show()
# plt.grid()























# ########################### Example - view and manipulate a fits file ###########################
# # image_file = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/SIXTE_examples/mcrab.fits"
# image_file = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/SIXTE_examples/img_mcrab.fits"
# hdu_list = fits.open(image_file)
# hdu_list.info()
# image_data = fits.getdata(image_file)
# print(type(image_data))
# print(image_data.shape)
# hdu_list.close()
# plt.imshow(image_data, cmap='gray')
# plt.colorbar()
# plt.show()
# print('Min:', np.min(image_data))
# print('Max:', np.max(image_data))
# print('Mean:', np.mean(image_data))
# print('Stdev:', np.std(image_data))
















# ########################### Example - Plot posterior probability as function of like ratio for octants ###########################
#
# # Get the jsons for a particular set of data
# # JsonLoc = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/Randomized_SYNEX/0cut/*.json"
# JsonLoc = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/Randomized_SYNEX/RemakePaperPlots/Randomized_angles_spins_MRat_"
# GridJsonFiles = glob.glob(JsonLoc+"*_1d.json")
#
# # Tell you how many data files there are
# print(len(GridJsonFiles), "json files found.")
#
# # JsonFileAndPath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/MassGrid_Mpi_0cut_6e7.json"
# fulllnLikeRatios = {}
# fullOctantPostProbs = {}
# for JsonFileAndPath in GridJsonFiles:
#     lnLikeRatios,OctantPostProbs = SYU.GetOctantLikeRatioAndPostProb(JsonFileAndPath,source=None,detector=None)
#     skymodes = lnLikeRatios.keys()
#     for skymode in skymodes:
#         if skymode not in fulllnLikeRatios:
#             fulllnLikeRatios[skymode] = list([lnLikeRatios[skymode]])
#         else:
#             fulllnLikeRatios[skymode].append(lnLikeRatios[skymode])
#         if skymode not in fullOctantPostProbs:
#             fullOctantPostProbs[skymode] = list([OctantPostProbs[skymode]])
#         else:
#             fullOctantPostProbs[skymode].append(OctantPostProbs[skymode])
#
# # useful lists for plotting
# skymode_markers = ["o","x","+","d","s","v","^","<"]
#
# # Plot results
# imarker = 0
# for skymode in skymodes:
#     # Create lists to plot
#     Xs = np.log10(np.abs(fulllnLikeRatios[skymode])) # fulllnLikeRatios[skymode] #
#     Ys = np.log(fullOctantPostProbs[skymode]) # fullOctantPostProbs[skymode] #
#     plt.plot(Xs,Ys,skymode_markers[imarker],label=skymode)
#     plt.xlabel(r"$\log(\mid\ln(\mathcal{L}_{(a,b)})\mid)$")
#     # plt.xlabel(r"$\ln(\mathcal{L}_{(a,b)})$")
#     # plt.ylabel(r"P$_{(a,b)}$")
#     plt.ylabel(r"$\log(\mathrm{P}_{(a,b)})$")
#     plt.legend(ncol=2)
#     imarker+=1
# plt.show()
# plt.grid()



















# ########################### Example - Calculate kuiper for Idea Paper example system ###########################
#
# # Create options to Initialize the detector object
# LISA_kwargs = {"TDI":'TDIAET'}
#
# # NB Dec = beta and Ra = lambda.... 13600Mpc == z=1.7639249603176368
# Merger_kwargs = {"q": 10., "M": 5e6, "dist": 13600., "chi1": 0.,
#         "chi2": 0., "beta" : 0., "lambda" : 0., "Lframe":True, "DeltatL_cut":-9.*24.*60.*60.,
#         "inc": 40.3*np.pi/180., "psi": 0.,  "approximant" : 'IMRPhenomHM'}
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_kwargs)
#
# # Initialize the source object
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
# # Compute the xray flux
# Merger.GenerateEMFlux(LISA)
#
# # Now compute the CTR
# ARF_file_loc_name = 'XIFU_CC_BASELINECONF_2018_10_10.arf'
# Merger.GenerateCTR(ARF_file_loc_name=ARF_file_loc_name,gamma=1.7)
#
# # Plot the results to see if we get what they have in their paper
# Xs = [time/(60.*60.*24.) for time in Merger.xray_time]
# Ys = Merger.CTR
# import matplotlib.pylab as pylab
# params = {'legend.fontsize': 8, # 'x-large',
#          'axes.labelsize': 8, # 'x-large',
#          'xtick.labelsize': 4, # 'x-large',
#          'ytick.labelsize': 4, # 'x-large'}
#          'lines.markersize': 2}
# pylab.rcParams.update(params)
# plt.plot(Xs,Ys)
# ax = plt.gca()
# ax.set_yscale('log')
# plt.xlabel(r"Time [d]")
# plt.ylabel(r"CTR [ph s$^{-1}$]")
# plt.show()
# plt.grid()
#
# # Initiate Athena
# Athena_kwargs = {"FoView":1.*(np.pi/180.)**2,"T_lat":10000.,"T_init":0.,"slew_rate":None} # slew_rate=None for no gaps for slew
# Athena = SYDs.Athena(**Athena_kwargs)
#
# # Start tiling
# dataFile = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/IdeaPaperSystem_9d.h5"
# Athena.TileSkyArea(dataFile,overlap=0.5,SAVE_SOURCE_EM_PROPERTIES_IN_TILE_JSON=True)
#
# # Plot tiles
# SYU.PlotTilesArea(Athena.TileJsonName,n_tiles=50)
#
# # Tiling file name
# TileJsonName = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/Tile_files/IdeaPaperSystem_9d_tiled.json"
#
# # Get kuiper detection pvals
# Athena.GetKuiper(Athena.TileJsonName,source=Merger)
# # Athena.GetKuiper(TileJsonName)










# ########################### Example - Run PSO with SYNEX / lisabeta ###########################
# # Create options to Initialize the detector object
# LISA_kwargs = {"TDI":'TDIAET'}
#
# # Create options to Initialize the source object
# # NB Dec = beta and Ra = lambda
# Merger_kwargs = {"q": 1.1, "M": 1.e6, "z": 1.8, "chi1": 0.9,
#    "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3., "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM'} # 'IMRPhenomD' }
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_kwargs)
#
# # Initialize the source object
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
# # Example prior object for ptemcee in lisabeta
# # inference_params = {
# #   "infer_params": ["chi1", "chi2", "dist", "inc", "phi", "lambda", "beta", "psi"],
# #   "params_range": [[-0.5, 1.], [-0.5, 1.], [5.0e3, 2.0e5], [0.,np.pi], [-np.pi, np.pi], [-np.pi, np.pi], [-np.pi/2.,np.pi/2.], [0.,np.pi]],
# #   "prior_type": ["uniform", "uniform", "uniform", "sin", "uniform", "uniform", "cos", "uniform"],
# #   "wrap_params": [False, False, False, False, True, True, False, True]
# # }
# inference_params = {
#   "infer_params": ["x", "y"],
#   "params_range": [[-5., 5.], [-5., 5.]],
#   "prior_type": ["uniform", "uniform"], # Not sure how to deal with these details in PSO yet... Maybe some kind of heating?
#   "wrap_params": [False, False]
# }
#
# # Define cost function for optimization routine
# def fitness_function(ps): # (source, detector):
#     # kwargs = {} # extra parameters for changing source / detector / whatever parameters at runtime
#     # fishercov = GetFisher_smbh(source, detector, **kwargs)
#     # SkyArea = lisatools.sky_area_cov(fishercov, sq_deg=True, n_sigma=None, prob=0.90)
#     f = []
#     for p in ps:
#         x,y = p[0],p[1]
#         if np.sqrt(x**2+y**2)<=1.:
#             f.append(10.)
#         else:
#             f.append(x**2 + (y+1.)**2 - 5.*np.cos(1.5*x+1.5) - 3.*np.cos(2.*y-1.5))
#     return f
#
# # Run with empty kwargs (for now)
# kwargs={}
# SYU.RunPSO(Merger, LISA, fitness_function=fitness_function, max_iter=200,
#                 NumSwarms=10, priors=inference_params, **kwargs)
#
# # Add now competition term to interact between swarms
#
#
#
# # Run PSO
# while PSO_Class.is_running:
#     PSO_Class.next()
#     # if PSO_Class.iter % 10 == 0:
#     #     PlotPSOSnap(PSO_Class)
#     PlotPSOSnap(PSO_Class)



























# ########################### Example - Plot inference results for a range of grid tests ###########################
#
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/Randomized_SYNEX/1d/Randomized_SYNEX_1d_" # /BetaGridTest_NoMpi_"
# jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/Randomized_SYNEX/1d/Randomized_SYNEX_1d_" # /BetaGridTest_NoMpi_"
# import glob
# GridJsons = glob.glob(jsonFilePath+"*.json")
# # Gridh5pys = glob.glob(DataFilePath+"*.h5")
#
# PostVals = {}
# for GridJson in GridJsons:
#     # find the coresponding data file
#     x = GridJson.split("_")[-1]
#     Gridh5py = DataFilePath + x[:-5] + ".h5"
#     print(Gridh5py)
#     posteriors_full = plotutils.load_params_posterior_lisa_smbh(GridJson, Gridh5py, format='multiemcee', load_fisher=True)
#     fig = plotutils.corner_plot(posteriors_full['injparams_Lframe'], posteriors_full['post'], add_posteriors=None, output=False, output_dir=None, output_file=None, histograms=True, fisher=False, fishercov=posteriors_full['fishercov'], Lframe=True, params=['chi1', 'chi2', 'dist', 'inc', 'phi', 'lambda', 'beta', 'psi'], color=plotutils.plotpalette[0], cov_color='k', add_colors=[plotutils.plotpalette[1]], show_truths=True, truth_color='k', bins=25, show_histograms=False, plot_datapoints=False);
#     plt.show()
#     # SYU.PlotInferenceData(Gridh5py, SaveFig=False)
#     # SYU.PlotInferenceLambdaBeta(Gridh5py, SaveFig=False)
#     # PostVals[MtotStr] = SYU.GetPosteriorVal(DataFile)









# ########################## Example - Simple plot data file inference results ###########################
# bins = 50
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/Randomized_SYNEX/4hr/Randomized_SYNEX_4hr_22.h5"
# JsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/Randomized_SYNEX/4hr/Randomized_SYNEX_4hr_22.json"
# posteriors_full = plotutils.load_params_posterior_lisa_smbh(JsonFilePath, DataFilePath, format='multiemcee', load_fisher=True)
# fig = plotutils.corner_plot(posteriors_full['injparams_Lframe'], posteriors_full['post'], add_posteriors=None, output=False, output_dir=None, output_file=None, histograms=True, fisher=False, fishercov=posteriors_full['fishercov'], Lframe=True, params=['chi1', 'chi2', 'dist', 'inc', 'phi', 'lambda', 'beta', 'psi'], color=plotutils.plotpalette[0], cov_color='k', add_colors=[plotutils.plotpalette[1]], show_truths=True, truth_color='k', bins=25, show_histograms=False, plot_datapoints=False);
# plt.show()
# SYU.PlotInferenceData(DataFilePath, SaveFig=False)
# # SYU.PlotInferenceLambdaBeta(DataFilePath, bins=bins, SkyProjection=False, SaveFig=False)













# ########################## Example - Plot sky loc area histograms for each cut time on one graph ###########################
#
# JsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/Randomized_SYNEX/"
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/Randomized_SYNEX/"
#
# CutsToPlot = ["0cut","4hr","1d","1wk"]
#
# for iCut in range(len(CutsToPlot)):
#     JsonFilesToSearch = JsonFilePath + CutsToPlot[iCut] + "/"
#     JsonFiles=glob.glob(JsonFilesToSearch+"*.json")
#     # JsonFiles = [JsonFilePath + "0cut/Randomized_SYNEX_0cut_7.json"]
#     TotAreas = []
#     for JsonFile in JsonFiles:
#         DataFileLocAndName = DataFilePath + CutsToPlot[iCut] + "/" + (JsonFile.split('/')[-1]).split('.')[0] + ".h5"
#         try:
#             TotArea = SYU.GetTotSkyAreaFromPostData(DataFileLocAndName,ConfLevel=0.9,bins=100)
#             TotAreas.append(TotArea)
#         except:
#             print("File load failed - doesn't exist.")
#     bins = np.logspace(np.log10(6.),np.log10(6e4), 20)
#     TotAreas = [(TotArea*(180./np.pi)**2) for TotArea in TotAreas]
#     plt.hist(TotAreas, bins=bins, label=CutsToPlot[iCut])
# plt.legend()
# plt.xlabel(r"$\log_{10}(\Omega)\,$[$\log_{10}(\mathrm{deg}^2)$]")
# plt.ylabel(r"Count")
# ax = plt.gca()
# ax.set_xscale('log')
# plt.show()
















# ########################## Example - Plot sky loc area vs. parameter scatters for each cut time on one graph ###########################
#
# JsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/Randomized_SYNEX/"
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/Randomized_SYNEX/"
#
# CutsToPlot = ["0cut","4hr","1d","1wk"]
#
# ParamsToPlot = ["M", "q", "chi1", "chi2", "lambda", "beta"] # , "phi", "psi"]
#
# for ParamToPlot in ParamsToPlot:
#     for iCut in range(len(CutsToPlot)):
#         Xs = []
#         JsonFilesToSearch = JsonFilePath + CutsToPlot[iCut] + "/"
#         JsonFiles=glob.glob(JsonFilesToSearch+"*.json")
#         # JsonFiles = [JsonFilePath + "0cut/Randomized_SYNEX_0cut_7.json"]
#         TotAreas = []
#         for JsonFile in JsonFiles:
#             DataFileLocAndName = DataFilePath + CutsToPlot[iCut] + "/" + (JsonFile.split('/')[-1]).split('.')[0] + ".h5"
#             try:
#                 TotArea = SYU.GetTotSkyAreaFromPostData(DataFileLocAndName,ConfLevel=0.9,bins=100)
#                 TotAreas.append(TotArea)
#                 with open(JsonFile) as f:
#                     JsonData = json.load(f)
#                 f.close()
#                 if ParamToPlot == "M":
#                     Xs.append(JsonData["source_params"]["m1"]+JsonData["source_params"]["m2"])
#                 elif ParamToPlot == "q":
#                     Xs.append(JsonData["source_params"]["m1"]/JsonData["source_params"]["m2"])
#                 else:
#                     Xs.append(JsonData["source_params"][ParamToPlot])
#             except:
#                 print("File load failed - doesn't exist.")
#         TotAreas = [(TotArea*(180./np.pi)**2) for TotArea in TotAreas]
#         plt.scatter(Xs, TotAreas, label=CutsToPlot[iCut])
#     plt.legend()
#     plt.ylabel(r"$\Omega\,[\mathrm{deg}^2$]")
#     plt.xlabel(ParamToPlot)
#     ax = plt.gca()
#     ax.set_yscale('log')
#     if ParamToPlot == "M" or ParamToPlot == "q":
#         ax.set_xscale('log')
#     plt.show()























# ########################## Example - Reproduce Athena doc plots ###########################
#
# DO_FULL_TIM_PLOTS = False
# DO_1d_PLOT = True
#
# JsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/Randomized_SYNEX/RemakePaperPlots/"
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/Randomized_SYNEX/RemakePaperPlots/"
#
# CutsToPlot = ["1d"] # ["0cut","1min","1hr","5hr", "10hr", "1d", "3d", "1wk", "1mon"]
# CutsToPlot_loc = [-24.*60.*60.] # [0.,-60.,-60.*60.,-5.*60.*60., -10.*60.*60., -24.*60.*60., -3.*24.*60.*60., -7.*24.*60.*60, -30.*24.*60.*60.]
# SystemNumbers = [str(ii) for ii in range(1,13)] #
#
# if DO_FULL_TIM_PLOTS:
#     for SystemNumber in SystemNumbers:
#         Xs = []
#         JsonFilesToSearch = JsonFilePath + "Randomized_angles_spins_MRat_" + str(SystemNumber)
#         JsonFiles=glob.glob(JsonFilesToSearch+"*.json")
#         # JsonFiles = [JsonFilePath + "0cut/Randomized_SYNEX_0cut_7.json"]
#         TotAreas = []
#         SNRs = []
#         for JsonFile in JsonFiles:
#             DataFileLocAndName = DataFilePath + (JsonFile.split('/')[-1]).split('.')[0] + ".h5"
#             with open(JsonFile) as f:
#                 JsonData = json.load(f)
#             f.close()
#             if JsonData["waveform_params"]["DeltatL_cut"] == None:
#                 Xs.append(0.)
#             else:
#                 Xs.append(JsonData["waveform_params"]["DeltatL_cut"])
#             JsonData["source_params"]["DeltatL_cut"] = JsonData["waveform_params"]["DeltatL_cut"]
#             Merger = SYSs.SMBH_Merger(**JsonData["source_params"])
#             LISA = SYDs.LISA(**JsonData["waveform_params"])
#             SNR, freqs, h_full, Snvals = SYU.ComputeSNR(Merger, LISA, freqs=None, Lframe=False, ReturnAllVariable=True)
#             SNRs.append(SNR)
#             # TotArea = SYU.GetTotSkyAreaFromPostData(DataFileLocAndName,ConfLevel=0.90)
#             # TotAreas.append(TotArea)
#         # TotAreas = [TotArea*(180./np.pi)**2 for TotArea in TotAreas]
#         # plt.scatter(Xs, TotAreas, label=str(SystemNumber))
#         plt.scatter(Xs, SNRs, label=str(SystemNumber))
#     plt.xlabel(r"Pre-merger Cut Time [s (to merger)]")
#     # plt.ylabel(r"$\Omega\,[\mathrm{deg}^2$]")
#     plt.ylabel(r"S/N")
#     ax = plt.gca()
#     ax.set_yscale('log')
#     plt.xscale('symlog')
#     plt.xticks(CutsToPlot_loc, CutsToPlot, rotation=60)
#     plt.xlim([np.min(CutsToPlot_loc), np.max(CutsToPlot_loc)])
#     plt.show()
#
# if DO_1d_PLOT:
#     Xs = []
#     JsonFilesToSearch = JsonFilePath + "Randomized_angles_spins_MRat_"
#     JsonFiles=glob.glob(JsonFilesToSearch+"*_1d.json")
#     print(len(JsonFiles),"files found")
#     # JsonFiles = [JsonFilePath + "0cut/Randomized_SYNEX_0cut_7.json"]
#     TotAreas = []
#     SNRs = []
#     for JsonFile in JsonFiles:
#         DataFileLocAndName = DataFilePath + (JsonFile.split('/')[-1]).split('.')[0] + ".h5"
#         with open(JsonFile) as f:
#             JsonData = json.load(f)
#         f.close()
#         if JsonData["waveform_params"]["DeltatL_cut"] == None:
#             Xs.append(0.)
#         else:
#             Xs.append(JsonData["waveform_params"]["DeltatL_cut"])
#         JsonData["source_params"]["DeltatL_cut"] = JsonData["waveform_params"]["DeltatL_cut"]
#         Merger = SYSs.SMBH_Merger(**JsonData["source_params"])
#         LISA = SYDs.LISA(**JsonData["waveform_params"])
#         SNR, freqs, h_full, Snvals = SYU.ComputeSNR(Merger, LISA, freqs=None, Lframe=False, ReturnAllVariable=True)
#         SNRs.append(SNR)
#         TotArea = SYU.GetTotSkyAreaFromPostData(DataFileLocAndName,ConfLevel=0.90)
#         TotAreas.append(TotArea)
#     print(len(TotAreas), "points to histogram")
#     TotAreas = [TotArea*(180./np.pi)**2 for TotArea in TotAreas]
#     # hist_plot_bins = np.logspace(np.log10(min(SNRs)),np.log10(max(SNRs)),25)
#     hist_plot_bins = np.logspace(np.log10(min(TotAreas)),np.log10(max(TotAreas)),25)
#     plt.hist(TotAreas,hist_plot_bins)
#     # plt.hist(SNRs,hist_plot_bins)
#     plt.xlabel(r"$\Omega\,[\mathrm{deg}^2$]")
#     # plt.xlabel(r"S/N")
#     ax = plt.gca()
#     ax.set_xscale('log')
#     plt.show()


























########################## Example - Plot sky loc areas for four cut times ###########################
#
# FilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/MassGrid_Mpi_"
# JsonFileEnd = "_2e6.json"
#
# CutsToPlot = ["0cut","4hr","1d","1wk"]
#
# for iCut in range(len(CutsToPlot)):
#     bins = bins_array[iCut]
#     JsonFilePath = FilePath + CutsToPlot[iCut] + JsonFileEnd
#     DataFilePath = JsonFilePath[:-5] + ".h5"
#     print(DataFilePath)
#     # posteriors_full = plotutils.load_params_posterior_lisa_smbh(JsonFilePath, DataFilePath, format='multiemcee', load_fisher=True)
#     # fig = plotutils.corner_plot(posteriors_full['injparams_Lframe'], posteriors_full['post'], add_posteriors=None, output=False, output_dir=None, output_file=None, histograms=True, fisher=False, fishercov=posteriors_full['fishercov'], Lframe=True, params=['chi1', 'chi2', 'dist', 'inc', 'phi', 'lambda', 'beta', 'psi'], color=plotutils.plotpalette[0], cov_color='k', add_colors=[plotutils.plotpalette[1]], show_truths=True, truth_color='k', bins=25, show_histograms=False, plot_datapoints=False);
#     # plt.show()
#     # SYU.PlotInferenceData(DataFilePath, SaveFig=False)
#     SYU.PlotInferenceLambdaBeta(DataFilePath, bins=bins, SkyProjection=False, SaveFig=False)












# ######################### Example - Plot sky areas using 2D histograms versus No. of bins ###########################
# # there is a trade off between resoltuion and scarcity of data that gives a plot
# # from this section with a min around 110 bins. Therefore we use this in future histograms for sky area from 2D hists.
# DataFile = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/MassGrid_Mpi_0cut_3e7.h5"
# jsonFile = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/MassGrid_Mpi_0cut_3e7.json"
#
# Bins = [ii for ii in range(25,135,5)]
# out = []
# for bins in Bins:
#     res = SYU.GetTotSkyAreaFromPostData(DataFile,ConfLevel=0.95,bins=bins)*180.**2/(np.pi*np.pi)
#     # print(res)
#     out.append(res)
#
# plt.plot(Bins, out)
# plt.xlabel(r"No. Bins")
# plt.ylabel(r"$\Delta \Omega \, $[$\mathrm{deg}^2$]")
# plt.show()















# ######################### Example - Get sky areas using 2D histograms ###########################
#
# FilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/MassGrid_Mpi_"
#
# CutsToPlot = ["0cut","4hr","1d","1wk"]
#
# markers = ["+","x","<","^"]
#
# for iCut in range(len(CutsToPlot)):
#     # set bins
#     # bins = bins_array[iCut]
#
#     # Get all json file names (json not h5 directly to avoid including '_raw' files too)
#     GridJsonFiles=glob.glob(FilePath+CutsToPlot[iCut]+"*.json")
#
#     # Loop over each file
#     TotAreas=[]
#     TotMasses=[]
#     for GridJsonFile in GridJsonFiles:
#         # Create h5 file name
#         InfDataPath = '/'.join(GridJsonFile.split('/')[:-2]) + "/inference_data/"
#         GridDataFile = str(InfDataPath + '_'.join((GridJsonFile.split('/')[-1]).split('_')[:-1]) + "_" + (GridJsonFile.split('_')[-1]).split('.')[0] + ".h5")
#         print(GridDataFile)
#
#         # Grab the total area for each file
#         # TotAreas.append(bins**(-0.5)*SYU.GetTotSkyAreaFromPostData(GridDataFile,ConfLevel=0.9,bins=bins)*180.**2/(np.pi*np.pi))
#         TotAreas.append(SYU.GetTotSkyAreaFromPostData(GridDataFile,ConfLevel=0.9)*180.**2/(np.pi*np.pi))
#         print(SYU.GetTotSkyAreaFromPostData(GridDataFile,ConfLevel=0.9)*180.**2/(np.pi*np.pi))
#
#         # Grab the total mass data
#         [_, inj_param_vals, _, _] = SYU.read_h5py_file(GridDataFile)
#         TotMasses.append(inj_param_vals["source_params_Lframe"]["m1"][0]+inj_param_vals["source_params_Lframe"]["m2"][0])
#
#     # Plot
#     # plt.scatter(np.log10(TotMasses),np.log10(TotAreas),marker=markers[iCut],label=CutsToPlot[iCut])
#     plt.scatter(TotMasses,TotAreas,marker=markers[iCut],label=CutsToPlot[iCut])
#
# # Finish up presentation of plot
# ax = plt.gca()
# # plt.xlabel(r"log$_{10}(\mathrm{M}_{\mathrm{tot}}) \, [\mathrm{M}_{\odot}]$")
# # plt.ylabel(r"log$_{10}(\Omega) \, [\mathrm{rad}^2]$")
# plt.xlabel(r"$\mathrm{M}_{\mathrm{tot}} \, [\mathrm{M}_{\odot}]$")
# plt.ylabel(r"$\Omega \, [\mathrm{deg}^2]$")
# ax.set_xscale('log')
# ax.set_yscale('log')
# plt.legend()
# plt.show()
# plt.grid()























# ########################### Example - Plot list of grid runs ###########################
#
# # DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/SylvainTest_"
# # jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/SylvainTest_"
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/MassGrid_Mpi_0cut_"
# jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/MassGrid_Mpi_0cut_"
# import glob
# GridJsons = glob.glob(jsonFilePath+"*.json")
#
# # # find the coresponding data file
# # GridJson = jsonFilePath + "1.json"
# # x = GridJson.split("_")[-1]
# # Gridh5py = DataFilePath + x[:-5] + ".h5"
# # SYU.PlotInferenceData(Gridh5py, SaveFig=False)
# # SYU.PlotInferenceLambdaBeta(Gridh5py, SkyProjection=False, SaveFig=False)
# #
# for GridJson in GridJsons:
#     # find the coresponding data file
#     print("GridJson:", GridJson)
#     x = GridJson.split("_")[-1]
#     Gridh5py = DataFilePath + x[:-5] + ".h5"
#     print(Gridh5py)
#     # print(Gridh5pys)
#     # posteriors_full = plotutils.load_params_posterior_lisa_smbh(GridJson, Gridh5py, format='multiemcee', load_fisher=True)
#     # fig = plotutils.corner_plot(posteriors_full['injparams_Lframe'], posteriors_full['post'], add_posteriors=None, output=False, output_dir=None, output_file=None, histograms=True, fisher=False, fishercov=posteriors_full['fishercov'], Lframe=True, params=['chi1', 'chi2', 'dist', 'inc', 'phi', 'lambda', 'beta', 'psi'], color=plotutils.plotpalette[0], cov_color='k', add_colors=[plotutils.plotpalette[1]], show_truths=True, truth_color='k', bins=25, show_histograms=False, plot_datapoints=False);
#     # plt.show()
#     # SYU.PlotInferenceData(Gridh5py, SaveFig=False)
#     SYU.PlotInferenceLambdaBeta(Gridh5py, SkyProjection=True, SaveFig=False)





















# ########################## Example - Run FoM over range ###########################
# # Create options to Initialize the detector object
# LISA_kwargs = {"TDI":'TDIAET'}
#
# # Set the cut time for the premerger portion only
# T_str = "4hr"
# T_obs_end_to_merger = -4.*60.*60.
#
# # NB Dec = beta and Ra = lambda
# Merger_kwargs = {"q": 1.1, "M": 6000000., "z": 3., "chi1": 0.9,
#         "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3.,
#         "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM',
#         "DeltatL_cut":T_obs_end_to_merger,
#         "minf":6e-6}
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_kwargs)
#
# # Initialize the source object
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
# # FoM Variables
# # FoM_loop_kwargs = {"inc":[0.000000001,np.pi-0.000000001,10], "maxf": [1.e-4,1.,9]} # "beta": [-np.pi/2., np.pi/2., 70]} # {"M": [1.e6,5.e8,70],"z": [1.,9.,70]} # {"z": [1.,9.,50]} # {"M": [1.e6,5.e8,50]}
# FoM_loop_kwargs = {"M": [1.e6,5.e8,50], "z": [1.,9.,50]}
#
# FigureOfMerit = "SkyModeLikelihoodsBETA" # "SkyModeLikelihoodsLAMBDA" # "SNR" # "SkyArea" #
# RunGrid = True
#
# # Save the results?
# SaveOutput = True
#
# t1 = time.time()
# FoMOverRange = SYU.RunFoMOverRange(Merger,LISA,FoM_loop_kwargs,FigureOfMerit=FigureOfMerit,RunGrid=RunGrid)
# t2 = time.time()
# print("Time for SkyLoc Calc = " + str(t2-t1) + "s")
#
# if FoMOverRange["IsGrid"]:
#     # Get the variables
#     Vars = FoMOverRange["LoopVariables"]
#
#     # Save the data
#     if SaveOutput:
#         import json
#         json_file = "FoMOverRange_" + FoMOverRange["FoM"] + "_" + Vars[0] + "_" + Vars[1] + "_" + T_str + ".json"
#         with open(json_file, 'w') as f:
#             json.dump(FoMOverRange, f, indent=2)
#         f.close()
#
#     # Change the sign of things if needed
#     if Vars[0]=="DeltatL_cut":
#         FoMOverRange[Vars[0]+"_xs"] = [FoM/(60.*60.) for FoM in FoMOverRange[Vars[0]+"_xs"]]
#     if Vars[1]=="DeltatL_cut":
#         FoMOverRange[Vars[1]+"_xs"] = [FoM/(60.*60.) for FoM in FoMOverRange[Vars[1]+"_xs"]]
#
#     # Plot the data
#     X,Y = np.meshgrid(FoMOverRange[Vars[0]+"_xs"], FoMOverRange[Vars[1]+"_xs"])
#     if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA" or FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
#         Z = np.log10(abs(np.array(FoMOverRange["grid_ys"])+20.))
#     else:
#         Z = np.log10(FoMOverRange["grid_ys"])
#
#     fig, ax = plt.subplots(constrained_layout=True)
#     im = ax.pcolormesh(X, Y, Z, shading='gouraud', vmin=Z.min(), vmax=Z.max())
#     cbar = fig.colorbar(im, ax=ax)
#     if Vars[0] == "M":
#         ax.set_xscale('log')
#         plt.xlabel(r'M$_{\mathrm{tot}} \; (M_{\odot})$')
#     elif Vars[0] == "inc":
#         plt.xlabel(r'$\iota \; (\mathrm{deg.})$')
#     elif Vars[0] == "maxf":
#         ax.set_xscale('log')
#         plt.xlabel(r'$\mathrm{f}_{\mathrm{max}} \; (\mathrm{Hz})$')
#     elif Vars[0] == "beta":
#         plt.xlabel(r'$\beta_{SSB} \; (\mathrm{deg.})$')
#     elif Vars[0] == "z":
#         plt.xlabel(r'z')
#     elif Vars[0] == "DeltatL_cut":
#         ax.set_xscale('log')
#         plt.xlabel(r'$\Delta T_{L,cut} \; $(hr)')
#
#     if Vars[1] == "M":
#         ax.set_yscale('log')
#         plt.ylabel(r'M$_{\mathrm{tot}} \; (M_{\odot})$')
#     elif Vars[1] == "inc":
#         plt.ylabel(r'$\iota \; (\mathrm{deg.})$')
#     elif Vars[1] == "maxf":
#         ax.set_yscale('log')
#         plt.ylabel(r'$\mathrm{f}_{\mathrm{max}} \; (\mathrm{Hz})$')
#     elif Vars[1] == "beta":
#         plt.ylabel(r'$\beta_{SSB} \; (\mathrm{deg.})$')
#     elif Vars[1] == "z":
#         plt.ylabel(r'z')
#     elif Vars[1] == "DeltatL_cut":
#         ax.set_yscale('log')
#         plt.ylabel(r'$\Delta T_{L,cut} \; $(hr)')
#
#     if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA":
#         cbar.set_label(r'$\log_{10}(|\mathcal{B}_{\beta}|)$')
#     elif FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
#         cbar.set_label(r'$\log_{10}(|\mathcal{B}_{\lambda}|)$')
#     else:
#         cbar.set_label(r'$\log_{10}(\Delta \Omega \; (\mathrm{sq. deg.}))$')
#
#     plt.show()
# else:
#     for var in FoMOverRange["LoopVariables"]:
#         fig = plt.figure()
#         ax = plt.gca()
#
#         if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA" or FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
#             Z = np.log10(abs(np.array(FoMOverRange[var+"_ys"])+20.))
#         else:
#             Z = np.log10(FoMOverRange[var+"_ys"])
#
#         plt.plot(FoMOverRange[var+"_xs"], Z)
#         if var == "M":
#             ax.set_xscale('log')
#             plt.xlabel(r'M$_{\mathrm{tot}} \; (M_{\odot})$')
#         elif var == "inc":
#             plt.xlabel(r'$\iota \; (\mathrm{deg.})$')
#         elif var == "maxf":
#             ax.set_xscale('log')
#             plt.xlabel(r'$\mathrm{f}_{\mathrm{max}} \; (\mathrm{Hz})$')
#         elif var == "beta":
#             plt.xlabel(r'$\beta_{SSB} \; (\mathrm{deg.})$')
#         elif var == "z":
#             plt.xlabel(r'z')
#         elif var == "T_obs_end_to_merger":
#             ax.set_xscale('log')
#             plt.xlabel(r'T$_{ontm} \; $(s)')
#
#         if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA":
#             plt.ylabel(r'$\log_{10}(|\mathcal{B}_{\beta}|)$')
#         elif FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
#             plt.ylabel(r'$\log_{10}(|\mathcal{B}_{\lambda}|)$')
#         else:
#             plt.ylabel(r'$\log_{10}(\Delta \Omega \; (\mathrm{sq. deg.}))$')
#
#         plt.show()
#         plt.grid()
#
#         if SaveOutput:
#             json_file = "FoMOverRange_" + FoMOverRange["FoM"] + "_" + var + "_" + T_str + ".json"
#             with open(json_file, 'w') as f:
#                 json.dump(FoMOverRange, f, indent=2)
#             f.close()













# ########################## Example - Load FoM over range and plot with inlay infers ###########################
#
# json_file = "FoMOverRange_SkyModeLikelihoodsLambda_M_z_0cut.json"
# BF_lims = np.linspace(10.,30., 5)
# SaveFig = False
# for BF_lim in BF_lims:
#     SYU.PlotLikeRatioFoMFromJsonWithInferenceInlays(json_file, BF_lim, SaveFig)










# ########################## Example - Load FoM as func of Mtot and plot lines of z with inlay waveforms ###########################
# # Fags to normaize by Fisher elements?
# NormaliseByFisher = False
#
# # Create options to Initialize the source object
# LISA_kwargs = {"TDI":'TDIAET'}
#
# # NB Dec = beta and Ra = lambda
# # Merger_kwargs = {"q": 1.1, "M": 5130000., "z": 3.11, "chi1": 0.9,
# Merger_kwargs = {"q": 1.1, "M": 5130000., "z": 3.11, "chi1": 0.9,
#         "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3., "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM'}
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_kwargs)
#
# # Initialize the source object
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
# # Load data
# json_file = "FoMOverRange_z_M.json"
# with open(json_file) as f:
#     FoMOverRange = json.load(f)
# f.close()
# # print(data)
#
# # Check that you loaded a grided FoMOverRange file
# if not FoMOverRange["IsGrid"]:
#     raise ValueError("Loaded json file does not contain grid data. Check filename.")
#
# # Extract x data
# M_tot_xs = FoMOverRange["M_xs"]
# Red_xs = FoMOverRange["z_xs"]
#
# # Extract y data
# SkyArea_ys = np.array(FoMOverRange["grid_ys"]) # * (180./np.pi)**2
#
# # Get dimensions
# dims = np.shape(SkyArea_ys)
#
# # Loop over a handful of redshifts to plot lines
# ChosenReds = []
# for iRed in range(0,dims[1],5):
#     ys = [y for y in SkyArea_ys[:,iRed]]
#     label = r'z = %0.2f' % (Red_xs[iRed])
#     if NormaliseByFisher:
#         FisherElement = []
#         for mass in M_tot_xs:
#             from lisabeta.lisa.lisa_fisher import default_steps
#             import lisabeta.lisa.lisa_fisher as lisa_fisher
#             Par1 = "M"
#             Par2 = "M"
#             wftdi = SYU.GetSMBHGWDetection(Merger, LISA)
#             freqs = wftdi[(2,2)]["freq"]
#             BaseParamDict = {"M":mass}
#             param_dict, _, _ = SYU.ClassesToParams(Merger, LISA, CollectionMethod="Fisher",**BaseParamDict)
#             FisherElement.append(lisa_fisher.fisher_element(param_dict, Par1, Par2, default_steps[Par1], default_steps[Par2], freqs))
#         plt.plot(M_tot_xs,ys/np.array(FisherElement),label=label)
#     else:
#         plt.plot(M_tot_xs,ys,label=label)
#
#     # store the redshift plotted for getting inlay positions later
#     ChosenReds.append(Red_xs[iRed])
#
# # Find the ys positions of each mass inference plot to be inlay
# i0 = np.where(abs(np.array(M_tot_xs)-3.e6)==min(abs(np.array(M_tot_xs)-3.e6)))[0][0]
# i1 = np.where(abs(np.array(M_tot_xs)-1.e7)==min(abs(np.array(M_tot_xs)-1.e7)))[0][0]
# i2 = np.where(abs(np.array(M_tot_xs)-8.e7)==min(abs(np.array(M_tot_xs)-8.e7)))[0][0]
# iRed = np.where(abs(np.array(Red_xs)-ChosenReds[1])==min(abs(np.array(Red_xs)-ChosenReds[1])))[0][0]
# MassYs = [0.067, 0.062, 0.42]# [SkyArea_ys[i0,iRed],
#          #  SkyArea_ys[i1,iRed],
#          #  SkyArea_ys[i2,iRed]]
#
# # Log scale for mass and axis labels
# ax1 = plt.gca()
# ax1.set_xscale('log')
# plt.xlabel(r'$M_{tot} \; (M_{\odot})$')
# plt.ylabel(r'$\Delta \Omega \; (\mathrm{sq. deg.})$')
# if not NormaliseByFisher:
#     plt.ylim([0.,3.])
#     plt.xlim([8.e5,6.e8])
#
# # Add legend, grid and show
# plt.legend(loc='upper right')
# plt.grid()
#
# # Overlay inference plots
# MtotTestedStrs = ["3e6", "1e7", "8e7"]
# MtotTestedVals = [3.e6, 1.e7, 8.e7]
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/MassGridTest_NoMpi_"
# jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/MassGridTest_NoMpi_"
#
# # Positioning arguments for inlay inference plots
# pos_lefts = np.array([0.18, 0.45, 0.69])
# pos_bottoms = np.array([0.6, 0.6, 0.4])
# pos_widths = np.array([0.2, 0.2, 0.2])
# pos_heights = np.array([0.2, 0.2, 0.2])
#
# # Arrow coords and widths
# Arrow_xs = [2.e6, 2.e7, 1.01e8]
# Arrow_ys = [1.9, 1.9, 1.1]
# Arrow_dxs = [3.e6-Arrow_xs[0], 1.e7-Arrow_xs[1], 8.e7-Arrow_xs[2]]
# Arrow_dys = [MassYs[0]-Arrow_ys[0], MassYs[1]-Arrow_ys[1], MassYs[2]-Arrow_ys[2]]
#
# # Add arrows
# for iArrow in range(len(Arrow_xs)):
#     mpl.pyplot.arrow(Arrow_xs[iArrow], Arrow_ys[iArrow], Arrow_dxs[iArrow], Arrow_dys[iArrow])
#
# for iMtot in range(len(MtotTestedStrs)):
#     # Change base total mass
#     kwargs = {"M":MtotTestedVals[iMtot]}
#
#     # Get LISA Noise and full waveforms
#     SNR, freqs, h_full, Snvals = SYU.ComputeSNR(Merger, LISA, freqs=None, Lframe=False, ReturnAllVariable=True, **kwargs)
#
#     # Create new axes
#     ax = plt.axes([pos_lefts[iMtot], pos_bottoms[iMtot], pos_widths[iMtot], pos_heights[iMtot]], facecolor='y')
#
#     h_full_tot = 0.
#     Snvals_tot = 0.
#     for chan in Snvals.keys():
#         # # Plot waveform
#         # plt.loglog(freqs, np.sqrt(h_full[chan]*np.conj(h_full[chan])))
#         #
#         # # Plot noise
#         # plt.loglog(freqs, np.sqrt(Snvals[chan]))
#
#         # Totals
#         h_full_tot += h_full[chan]
#         Snvals_tot += Snvals[chan]
#
#     # Plot waveform
#     plt.loglog(freqs, np.sqrt(h_full_tot*np.conj(h_full_tot)))
#
#     # Plot noise
#     plt.loglog(freqs, np.sqrt(Snvals_tot))
#
#     # Labels and grid
#     plt.grid()
#     # plt.xlabel(labels[0])
#     # plt.ylabel(labels[5])
#
#     # Ranges
#     plt.ylim([1.e-23, 1.e-18])
#
#     # Move to position
#     plt.xticks([])
#     plt.yticks([])
#
# plt.show()
























# ########################## Example - Load FoM as func of z and plot lines of Mtot with inlay waveforms ###########################
# # Fags to normaize by Fisher elements?
# NormaliseByFisher = False
#
# # Create options to Initialize the source object
# LISA_kwargs = {"TDI":'TDIAET'}
#
# # NB Dec = beta and Ra = lambda
# # Merger_kwargs = {"q": 1.1, "M": 5130000., "z": 3.11, "chi1": 0.9,
# Merger_kwargs = {"q": 1.1, "M": 5130000., "z": 3.11, "chi1": 0.9,
#         "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3., "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM'}
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_kwargs)
#
# # Initialize the source object
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
# # Load data
# json_file = "FoMOverRange_z_M.json"
# with open(json_file) as f:
#     FoMOverRange = json.load(f)
# f.close()
# # print(data)
#
# # Check that you loaded a grided FoMOverRange file
# if not FoMOverRange["IsGrid"]:
#     raise ValueError("Loaded json file does not contain grid data. Check filename.")
#
# # Extract x data
# M_tot_xs = FoMOverRange["M_xs"]
# Red_xs = FoMOverRange["z_xs"]
#
# # Extract y data
# SkyArea_ys = np.array(FoMOverRange["grid_ys"]) # * (180./np.pi)**2
#
# # Get dimensions
# dims = np.shape(SkyArea_ys)
#
# # Loop over a handful of redshifts to plot lines
# ChosenMasses = []
# for iMass in range(0,dims[1],5):
#     ys = [y for y in SkyArea_ys[iMass,:]]
#     label = r'$M_{tot} = %d\times10^{6} \, M_{\odot}$' % (M_tot_xs[iMass]/1000000.)
#     if NormaliseByFisher:
#         FisherElement = []
#         for z in Red_xs:
#             from lisabeta.lisa.lisa_fisher import default_steps
#             import lisabeta.lisa.lisa_fisher as lisa_fisher
#             Par1 = "dist"
#             Par2 = "dist"
#             wftdi = SYU.GetSMBHGWDetection(Merger, LISA)
#             freqs = wftdi[(2,2)]["freq"]
#             BaseParamDict = {"z":z}
#             param_dict, _, _ = SYU.ClassesToParams(Merger, LISA, CollectionMethod="Fisher",**BaseParamDict)
#             FisherElement.append(lisa_fisher.fisher_element(param_dict, Par1, Par2, default_steps[Par1], default_steps[Par2], freqs))
#         if Par1 == "dist":
#             z_to_dist = cosmo.luminosity_distance(Merger.z).to("Mpc").value/Merger.z
#             dist_to_z = 1./z_to_dist
#             FisherElement = FisherElement/z_to_dist
#         if Par2 == "dist":
#             z_to_dist = cosmo.luminosity_distance(Merger.z).to("Mpc").value/Merger.z
#             dist_to_z = 1./z_to_dist
#             FisherElement = FisherElement/z_to_dist
#         plt.plot(Red_xs,ys*np.array(FisherElement),label=label)
#         plt.plot(Red_xs,ys,label=label)
#     else:
#         plt.plot(Red_xs,ys,label=label)
#
#     # store the redshift plotted for getting inlay positions later
#     ChosenMasses.append(M_tot_xs[iMass])
#
# # Find the ys positions of each mass inference plot to be inlay
# i0 = np.where(abs(np.array(Red_xs)-1.)==min(abs(np.array(Red_xs)-1.)))[0][0]
# i1 = np.where(abs(np.array(Red_xs)-5.)==min(abs(np.array(Red_xs)-5.)))[0][0]
# i2 = np.where(abs(np.array(Red_xs)-9.)==min(abs(np.array(Red_xs)-9.)))[0][0]
# iRed = np.where(abs(np.array(M_tot_xs)-ChosenMasses[1])==min(abs(np.array(M_tot_xs)-ChosenMasses[1])))[0][0]
# MassYs = [0.02, 0.1, 0.42]# [SkyArea_ys[i0,iRed],
#          #  SkyArea_ys[i1,iRed],
#          #  SkyArea_ys[i2,iRed]]
#
# # Log scale for mass and axis labels
# ax1 = plt.gca()
# plt.xlabel(r'z$ \; ( )$')
# plt.ylabel(r'$\Delta \Omega \; (\mathrm{sq. deg.})$')
# plt.ylim([0.,3.])
# plt.xlim([0.5,9.5])
#
# # Add legend, grid and show
# plt.legend(loc='upper right')
# plt.grid()
#
# # Overlay inference plots
# RedTestedStrs = ["1", "5", "9"]
# RedTestedVals = [1., 5., 9.]
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/RedGridTest_NoMpi_"
# jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/RedGridTest_NoMpi_"
#
# # Positioning arguments for inlay inference plots
# pos_lefts = np.array([0.16, 0.2, 0.54])
# pos_bottoms = np.array([0.3, 0.6, 0.46])
# pos_widths = np.array([0.2, 0.2, 0.2])
# pos_heights = np.array([0.2, 0.2, 0.2])
#
# # Arrow coords and widths
# Arrow_xs = [1.5, 3., 6.5]
# Arrow_ys = [0.75, 1.9, 1.4]
# Arrow_dxs = [1.-Arrow_xs[0], 5.-Arrow_xs[1], 9.-Arrow_xs[2]]
# Arrow_dys = [MassYs[0]-Arrow_ys[0], MassYs[1]-Arrow_ys[1], MassYs[2]-Arrow_ys[2]]
#
# # Add arrows
# for iArrow in range(len(Arrow_xs)):
#     mpl.pyplot.arrow(Arrow_xs[iArrow], Arrow_ys[iArrow], Arrow_dxs[iArrow], Arrow_dys[iArrow])
#
# for iRed in range(len(RedTestedStrs)):
#     # Change base total mass
#     kwargs = {"z":RedTestedVals[iRed]}
#
#     # Get LISA Noise and full waveforms
#     SNR, freqs, h_full, Snvals = SYU.ComputeSNR(Merger, LISA, freqs=None, Lframe=False, ReturnAllVariable=True, **kwargs)
#
#     # Create new axes
#     ax = plt.axes([pos_lefts[iRed], pos_bottoms[iRed], pos_widths[iRed], pos_heights[iRed]], facecolor='y')
#
#     h_full_tot = 0.
#     Snvals_tot = 0.
#     for chan in Snvals.keys():
#         # # Plot waveform
#         # plt.loglog(freqs, np.sqrt(h_full[chan]*np.conj(h_full[chan])))
#         #
#         # # Plot noise
#         # plt.loglog(freqs, np.sqrt(Snvals[chan]))
#
#         # Totals
#         h_full_tot += h_full[chan]
#         Snvals_tot += Snvals[chan]
#
#     # Plot waveform
#     plt.loglog(freqs, np.sqrt(h_full_tot*np.conj(h_full_tot)))
#
#     # Plot noise
#     plt.loglog(freqs, np.sqrt(Snvals_tot))
#
#     # Labels and grid
#     plt.grid()
#     # plt.xlabel(labels[0])
#     # plt.ylabel(labels[5])
#
#     # Ranges
#     plt.ylim([1.e-23, 1.e-18])
#
#     # Move to position
#     plt.xticks([])
#     plt.yticks([])
#
# plt.show()





























# ########################## Example - Load FoM as func of Mtot and plot lines of z with inlay infers ###########################
# # Fags to normaize by Fisher elements?
# NormaliseByFisher = True
#
# # Create options to Initialize the source object
# LISA_kwargs = {"TDI":'TDIAET'}
#
# # NB Dec = beta and Ra = lambda
# # Merger_kwargs = {"q": 1.1, "M": 5130000., "z": 3.11, "chi1": 0.9,
# Merger_kwargs = {"q": 1.1, "M": 5130000., "z": 3.11, "chi1": 0.9,
#         "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3., "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM'}
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_kwargs)
#
# # Initialize the source object
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
# # Load data
# json_file = "FoMOverRange_z_M.json"
# with open(json_file) as f:
#     FoMOverRange = json.load(f)
# f.close()
# # print(data)
#
# # Check that you loaded a grided FoMOverRange file
# if not FoMOverRange["IsGrid"]:
#     raise ValueError("Loaded json file does not contain grid data. Check filename.")
#
# # Extract x data
# M_tot_xs = FoMOverRange["M_xs"]
# Red_xs = FoMOverRange["z_xs"]
#
# # Extract y data
# SkyArea_ys = np.array(FoMOverRange["grid_ys"]) # * (180./np.pi)**2
#
# # Get dimensions
# dims = np.shape(SkyArea_ys)
#
# # Loop over a handful of redshifts to plot lines
# ChosenReds = []
# for iRed in range(0,dims[1],5):
#     ys = [y for y in SkyArea_ys[:,iRed]]
#     label = r'z = %0.2f' % (Red_xs[iRed])
#     if NormaliseByFisher:
#         FisherElement = []
#         for mass in M_tot_xs:
#             from lisabeta.lisa.lisa_fisher import default_steps
#             import lisabeta.lisa.lisa_fisher as lisa_fisher
#             Par1 = "M"
#             Par2 = "M"
#             wftdi = SYU.GetSMBHGWDetection(Merger, LISA)
#             freqs = wftdi[(2,2)]["freq"]
#             BaseParamDict = {"M":mass}
#             param_dict, _, _ = SYU.ClassesToParams(Merger, LISA, CollectionMethod="Fisher",**BaseParamDict)
#             FisherElement.append(lisa_fisher.fisher_element(param_dict, Par1, Par2, default_steps[Par1], default_steps[Par2], freqs))
#         plt.plot(M_tot_xs,ys/np.array(FisherElement),label=label)
#     # plt.plot(M_tot_xs,ys,label=label)
#
#     # store the redshift plotted for getting inlay positions later
#     ChosenReds.append(Red_xs[iRed])
#
# # Find the ys positions of each mass inference plot to be inlay
# i0 = np.where(abs(np.array(M_tot_xs)-3.e6)==min(abs(np.array(M_tot_xs)-3.e6)))[0][0]
# i1 = np.where(abs(np.array(M_tot_xs)-1.e7)==min(abs(np.array(M_tot_xs)-1.e7)))[0][0]
# i2 = np.where(abs(np.array(M_tot_xs)-8.e7)==min(abs(np.array(M_tot_xs)-8.e7)))[0][0]
# iRed = np.where(abs(np.array(Red_xs)-ChosenReds[1])==min(abs(np.array(Red_xs)-ChosenReds[1])))[0][0]
# MassYs = [0.067, 0.062, 0.42]# [SkyArea_ys[i0,iRed],
#          #  SkyArea_ys[i1,iRed],
#          #  SkyArea_ys[i2,iRed]]
#
# # Log scale for mass and axis labels
# ax1 = plt.gca()
# ax1.set_xscale('log')
# plt.xlabel(r'$M_{tot} \; (M_{\odot})$')
# plt.ylabel(r'$\Delta \Omega \; (\mathrm{sq. deg.})$')
# if not NormaliseByFisher:
#     plt.ylim([0.,3.])
#     plt.xlim([8.e5,6.e8])
#
# # Add legend, grid and show
# plt.legend(loc='upper right')
# plt.grid()
#
# # Overlay inference plots
# MtotTestedStrs = ["3e6", "1e7", "8e7"]
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/MassGridTest_NoMpi_"
# jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/MassGridTest_NoMpi_"
#
# # Positioning arguments for inlay inference plots
# pos_lefts = np.array([0.18, 0.45, 0.69])
# pos_bottoms = np.array([0.6, 0.6, 0.4])
# pos_widths = np.array([0.2, 0.2, 0.2])
# pos_heights = np.array([0.2, 0.2, 0.2])
#
# # Arrow coords and widths
# Arrow_xs = [2.e6, 2.e7, 1.01e8]
# Arrow_ys = [1.9, 1.9, 1.1]
# Arrow_dxs = [3.e6-Arrow_xs[0], 1.e7-Arrow_xs[1], 8.e7-Arrow_xs[2]]
# Arrow_dys = [MassYs[0]-Arrow_ys[0], MassYs[1]-Arrow_ys[1], MassYs[2]-Arrow_ys[2]]
#
# # Add arrows
# for iArrow in range(len(Arrow_xs)):
#     mpl.pyplot.arrow(Arrow_xs[iArrow], Arrow_ys[iArrow], Arrow_dxs[iArrow], Arrow_dys[iArrow])
#
# PostVals = {}
# for iMtot in range(len(MtotTestedStrs)):
#     DataFile = DataFilePath + MtotTestedStrs[iMtot] + ".h5"
#     jsonFile = jsonFilePath + MtotTestedStrs[iMtot] + ".json"
#     DataFileRaw = DataFile[:-3]+"_raw.h5"
#
#     # Make the inference plots from scratch since we don't pass back axes handles in the utils functions
#     DataFileLocAndName = DataFile
#     hist_n_bins=1000
#     [infer_params, inj_param_vals, static_params, meta_data] = SYU.read_h5py_file(DataFileLocAndName)
#     if not inj_param_vals: # Raw files at first didn't record this, so make sure it's there...
#         # Get the inj values from the processed data file instead
#         DataFileLocAndName_NotRaw = DataFileLocAndName[:-7] + ".h5"
#         [_,inj_param_vals,_,_] = read_h5py_file(DataFileLocAndName_NotRaw)
#     ndim = len(infer_params.keys())
#     labels = list(infer_params.keys())
#     if np.size(infer_params[labels[0]][0])>1:
#         nsamples = len(infer_params[labels[0]][0])
#     else:
#         nsamples = len(infer_params[labels[0]])
#     print("Posterior sample length: " + str(nsamples) + ", number of infered parameters: " + str(ndim))
#     # Grab data for infered parameters
#     data = np.empty([ndim,nsamples])
#     SampleModes = []
#     for ii in range(ndim):
#         if np.size(infer_params[labels[0]][0])>1:
#             data[ii][:] = infer_params[labels[ii]][0]
#         else:
#             data[ii][:] = infer_params[labels[ii]]
#         histn,histbins = np.histogram(data[ii,:], bins=hist_n_bins)
#         histbins = histbins[:-1] + 0.5*(histbins[2]-histbins[1]) # change from left bins edges to middle of bins
#         mode = histbins[histn==max(histn)]
#         if len(mode)>1:
#             mode = mode[1] # Doe now take the first on the list, but if there are several we need to work out what to do there...
#         SampleModes.append(mode)
#     data = np.transpose(np.array(data)) # should have shape [nsamples, ndim]
#
#     # Get injected values
#     InjParam_InjVals = []
#     for key in infer_params.keys():
#         InjParam_InjVals.append(inj_param_vals["source_params_Lframe"][key][0]) # Lframe is right.
#
#     # Sctter plot of labmda and beta posterior chains, marginalized over all other params
#     ax = plt.axes([pos_lefts[iMtot], pos_bottoms[iMtot], pos_widths[iMtot], pos_heights[iMtot]], facecolor='y')
#     plt.scatter(data[:,0], data[:,5], marker=".", linestyle="None")
#
#     # Add injected and mode vertical and horizontal lines
#     ax.axvline(InjParam_InjVals[0], color="g", linestyle=":")
#     ax.axvline(SampleModes[0], color="r", linestyle=":")
#     ax.axhline(InjParam_InjVals[5], color="g", linestyle=":")
#     ax.axhline(SampleModes[5], color="r", linestyle=":")
#
#     # Add points at injected and mode values
#     ax.plot(InjParam_InjVals[0], InjParam_InjVals[5], "sg")
#     ax.plot(SampleModes[0], SampleModes[5], "sr")
#
#     # Labels and grid
#     plt.grid()
#     plt.xlabel(labels[0])
#     plt.ylabel(labels[5])
#
#     # Move to position
#     plt.xticks([])
#     plt.yticks([])
#
# plt.show()











# ########################## Example - Load FoM as func of z and plot lines of Mtot with inlay infers ###########################
# # Fags to normaize by Fisher elements?
# NormaliseByFisher = True
#
# # Create options to Initialize the source object
# LISA_kwargs = {"TDI":'TDIAET'}
#
# # NB Dec = beta and Ra = lambda
# # Merger_kwargs = {"q": 1.1, "M": 5130000., "z": 3.11, "chi1": 0.9,
# Merger_kwargs = {"q": 1.1, "M": 5130000., "z": 3.11, "chi1": 0.9,
#         "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3., "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM'}
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_kwargs)
#
# # Initialize the source object
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
# # Load data
# json_file = "FoMOverRange_z_M.json"
# with open(json_file) as f:
#     FoMOverRange = json.load(f)
# f.close()
# # print(data)
#
# # Check that you loaded a grided FoMOverRange file
# if not FoMOverRange["IsGrid"]:
#     raise ValueError("Loaded json file does not contain grid data. Check filename.")
#
# # Extract x data
# M_tot_xs = FoMOverRange["M_xs"]
# Red_xs = FoMOverRange["z_xs"]
#
# # Extract y data
# SkyArea_ys = np.array(FoMOverRange["grid_ys"]) # * (180./np.pi)**2
#
# # Get dimensions
# dims = np.shape(SkyArea_ys)
#
# # Loop over a handful of redshifts to plot lines
# ChosenMasses = []
# for iMass in range(0,dims[1],5):
#     ys = [y for y in SkyArea_ys[iMass,:]]
#     label = r'$M_{tot} = %d\times10^{6} \, M_{\odot}$' % (M_tot_xs[iMass]/1000000.)
#     if NormaliseByFisher:
#         FisherElement = []
#         for z in Red_xs:
#             from lisabeta.lisa.lisa_fisher import default_steps
#             import lisabeta.lisa.lisa_fisher as lisa_fisher
#             Par1 = "dist"
#             Par2 = "dist"
#             wftdi = SYU.GetSMBHGWDetection(Merger, LISA)
#             freqs = wftdi[(2,2)]["freq"]
#             BaseParamDict = {"z":z}
#             param_dict, _, _ = SYU.ClassesToParams(Merger, LISA, CollectionMethod="Fisher",**BaseParamDict)
#             FisherElement.append(lisa_fisher.fisher_element(param_dict, Par1, Par2, default_steps[Par1], default_steps[Par2], freqs))
#         if Par1 == "dist":
#             z_to_dist = cosmo.luminosity_distance(Merger.z).to("Mpc").value/Merger.z
#             dist_to_z = 1./z_to_dist
#             FisherElement = FisherElement/z_to_dist
#         if Par2 == "dist":
#             z_to_dist = cosmo.luminosity_distance(Merger.z).to("Mpc").value/Merger.z
#             dist_to_z = 1./z_to_dist
#             FisherElement = FisherElement/z_to_dist
#         plt.plot(Red_xs,ys*np.array(FisherElement),label=label)
#         plt.plot(Red_xs,ys,label=label)
#     else:
#         plt.plot(Red_xs,ys,label=label)
#
#     # store the redshift plotted for getting inlay positions later
#     ChosenMasses.append(M_tot_xs[iMass])
#
# # Find the ys positions of each mass inference plot to be inlay
# i0 = np.where(abs(np.array(Red_xs)-1.)==min(abs(np.array(Red_xs)-1.)))[0][0]
# i1 = np.where(abs(np.array(Red_xs)-5.)==min(abs(np.array(Red_xs)-5.)))[0][0]
# i2 = np.where(abs(np.array(Red_xs)-9.)==min(abs(np.array(Red_xs)-9.)))[0][0]
# iRed = np.where(abs(np.array(M_tot_xs)-ChosenMasses[1])==min(abs(np.array(M_tot_xs)-ChosenMasses[1])))[0][0]
# MassYs = [0.02, 0.1, 0.42]# [SkyArea_ys[i0,iRed],
#          #  SkyArea_ys[i1,iRed],
#          #  SkyArea_ys[i2,iRed]]
#
# # Log scale for mass and axis labels
# ax1 = plt.gca()
# plt.xlabel(r'z$ \; ( )$')
# plt.ylabel(r'$\Delta \Omega \; (\mathrm{sq. deg.})$')
# plt.ylim([0.,3.])
# plt.xlim([0.5,9.5])
#
# # Add legend, grid and show
# plt.legend(loc='upper right')
# plt.grid()
#
# # Overlay inference plots
# RedTestedStrs = ["1", "5", "9"]
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/RedGridTest_NoMpi_"
# jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/RedGridTest_NoMpi_"
#
# # Positioning arguments for inlay inference plots
# pos_lefts = np.array([0.16, 0.2, 0.54])
# pos_bottoms = np.array([0.3, 0.6, 0.46])
# pos_widths = np.array([0.2, 0.2, 0.2])
# pos_heights = np.array([0.2, 0.2, 0.2])
#
# # Arrow coords and widths
# Arrow_xs = [1.5, 3., 6.5]
# Arrow_ys = [0.75, 1.9, 1.4]
# Arrow_dxs = [1.-Arrow_xs[0], 5.-Arrow_xs[1], 9.-Arrow_xs[2]]
# Arrow_dys = [MassYs[0]-Arrow_ys[0], MassYs[1]-Arrow_ys[1], MassYs[2]-Arrow_ys[2]]
#
# # Add arrows
# for iArrow in range(len(Arrow_xs)):
#     mpl.pyplot.arrow(Arrow_xs[iArrow], Arrow_ys[iArrow], Arrow_dxs[iArrow], Arrow_dys[iArrow])
#
# PostVals = {}
# for iRed in range(len(RedTestedStrs)):
#     DataFile = DataFilePath + RedTestedStrs[iRed] + ".h5"
#     jsonFile = jsonFilePath + RedTestedStrs[iRed] + ".json"
#     DataFileRaw = DataFile[:-3]+"_raw.h5"
#
#     # Make the inference plots from scratch since we don't pass back axes handles in the utils functions
#     DataFileLocAndName = DataFile
#     hist_n_bins=1000
#     [infer_params, inj_param_vals, static_params, meta_data] = SYU.read_h5py_file(DataFileLocAndName)
#     if not inj_param_vals: # Raw files at first didn't record this, so make sure it's there...
#         # Get the inj values from the processed data file instead
#         DataFileLocAndName_NotRaw = DataFileLocAndName[:-7] + ".h5"
#         [_,inj_param_vals,_,_] = read_h5py_file(DataFileLocAndName_NotRaw)
#     ndim = len(infer_params.keys())
#     labels = list(infer_params.keys())
#     if np.size(infer_params[labels[0]][0])>1:
#         nsamples = len(infer_params[labels[0]][0])
#     else:
#         nsamples = len(infer_params[labels[0]])
#     print("Posterior sample length: " + str(nsamples) + ", number of infered parameters: " + str(ndim))
#     # Grab data for infered parameters
#     data = np.empty([ndim,nsamples])
#     SampleModes = []
#     for ii in range(ndim):
#         if np.size(infer_params[labels[0]][0])>1:
#             data[ii][:] = infer_params[labels[ii]][0]
#         else:
#             data[ii][:] = infer_params[labels[ii]]
#         histn,histbins = np.histogram(data[ii,:], bins=hist_n_bins)
#         histbins = histbins[:-1] + 0.5*(histbins[2]-histbins[1]) # change from left bins edges to middle of bins
#         mode = histbins[histn==max(histn)]
#         if len(mode)>1:
#             mode = mode[1] # Doe now take the first on the list, but if there are several we need to work out what to do there...
#         SampleModes.append(mode)
#     data = np.transpose(np.array(data)) # should have shape [nsamples, ndim]
#
#     # Get injected values
#     InjParam_InjVals = []
#     for key in infer_params.keys():
#         InjParam_InjVals.append(inj_param_vals["source_params_Lframe"][key][0]) # Lframe is right.
#
#     # Sctter plot of labmda and beta posterior chains, marginalized over all other params
#     ax = plt.axes([pos_lefts[iRed], pos_bottoms[iRed], pos_widths[iRed], pos_heights[iRed]], facecolor='y')
#     plt.scatter(data[:,0], data[:,5], marker=".", linestyle="None")
#
#     # Add injected and mode vertical and horizontal lines
#     ax.axvline(InjParam_InjVals[0], color="g", linestyle=":")
#     ax.axvline(SampleModes[0], color="r", linestyle=":")
#     ax.axhline(InjParam_InjVals[5], color="g", linestyle=":")
#     ax.axhline(SampleModes[5], color="r", linestyle=":")
#
#     # Add points at injected and mode values
#     ax.plot(InjParam_InjVals[0], InjParam_InjVals[5], "sg")
#     ax.plot(SampleModes[0], SampleModes[5], "sr")
#
#     # Labels and grid
#     plt.grid()
#     plt.xlabel(labels[0])
#     plt.ylabel(labels[5])
#
#     # Move to position
#     plt.xticks([])
#     plt.yticks([])
#
# plt.show()











# ########################## Example - Plot FoM over range saved file ###########################
# JsonFileAndPath = "FoMOverRange_SkyModeLikelihoodsBeta_M_z_1wk.json"
# SaveFig = False
# BF_lim = 20.
# fig, ax = SYU.PlotLikeRatioFoMFromJson(JsonFileAndPath, BF_lim=BF_lim, SaveFig=SaveFig)
# plt.show()











# ########################### Example - Run ptmcee inference randomized angles and spins ###########################
#
# # Functions to randomize spins and angles
# def draw_random_angles(size=1):
#     inc = np.arccos(np.random.uniform(low=-1., high=1., size=size))
#     phi = np.random.uniform(low=-np.pi, high=np.pi, size=size)
#     lambd = np.random.uniform(low=-np.pi, high=np.pi, size=size)
#     beta = np.arcsin(np.random.uniform(low=-1., high=1., size=size))
#     psi = np.random.uniform(low=0., high=np.pi, size=size)
#     return np.array([inc, phi, lambd, beta, psi]).T
# def draw_random_spins(low=-1., high=1., size=1):
#     chi1 = np.random.uniform(low=low, high=high, size=size)
#     chi2 = np.random.uniform(low=low, high=high, size=size)
#     return np.array([chi1, chi2]).T
# def draw_random_masses(low=np.log10(7e5), high=np.log10(1e8), size=1):
#     m1 = 10.**(np.random.uniform(low=low, high=high, size=size))
#     m2 = 10.**(np.random.uniform(low=low, high=high, size=size))
#     return np.array([m1, m2]).T
#
# # Draw the random values
# n = 20
# rand_masses = draw_random_masses(size=n)
# rand_spins = draw_random_spins(size=n)
# rand_angles = draw_random_angles(size=n)
#
# # Create options to Initialize the detector object
# LISA_base_kwargs = {"TDI":'TDIAET'}
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_base_kwargs)
#
# # Set the time to merger cut if looking at pre-merger only - 'None' to ignore
# T_obs_end_to_merger = None # -4.*60.*60.
#
# # Create options to Initialize the source object
# # NB Dec = beta and Ra = lambda
# Merger_base_kwargs = {"q": 1.1, "M": 6e6, "z": 1., "chi1": 0.9,
#         "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3.,
#         "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM',
#         "Lframe":True, "DeltatL_cut":T_obs_end_to_merger}
#
# # Create dictionary of inference params along with their prior ranges and prior types
# # Note that there are default prior ranges for the anglular parameters: chi1, chi2, inc, phi, lambda, beta, psi
# inference_params = {
#   "infer_params": ["chi1", "chi2", "dist", "inc", "phi", "lambda", "beta", "psi"],
#   "params_range": [[-1., 1.], [-1., 1.], [5e3, 2e5], [0.,np.pi], [-np.pi, np.pi], [-np.pi, np.pi], [-np.pi/2.,np.pi/2.], [0.,np.pi]],
#   "prior_type": ["uniform", "uniform", "uniform", "sin", "uniform", "uniform", "cos", "uniform"],
#   "wrap_params": [False, False, False, False, True, True, False, True]
# }
#
# # Create a dictionary of additional parameters to run the inference with (the 'run_params' and 'waveform_params' options)
# # See SYNEX_PTMC.py for default dictionary that you can use as a base for what parameters are modifiable. Can also add a key for
# # any kwarg in source or detector classes and this will be modified in the run, but not updated in the class parameter list.
# RunTimekwargs = {"print_info": True, ### Run param options
#                 "n_walkers": 96, # must be greater than or equal to twice the inference cube dimension
#                 "n_iter": 8000,
#                 "burn_in": 5000, # Throw away at least 2/3 of it
#                 "autocor_method": "acor",
#                 "thin_samples": True, # for speed set this to False
#                 "TDI": "TDIAET", ### waveform param options. These are taken from the source and detector classes first, and then overridden here if specified
#                 "multimodal": True,
#                 "multimodal_pattern": "8modes",
#                 "p_jump": 0.5,
#                 "init_method": "fisher",
#                 "skip_fisher": False,
#                 "n_temps": 10}
#
# # Loop over the variable values. We can either run with the same class of detector and source
# # but using the new variables in the last optional dict, or we can recreate the classes on each iteration
# # of the loop. Since Alexis mentioned he wants to create a class on every iteration of the final product,
# # here we will do the same.
# for iiLoop in range(n):
#     # Copy then update the source class kwargs
#     Merger_kwargs = copy.deepcopy(Merger_base_kwargs)
#     Merger_kwargs["M"] = np.max(rand_masses[iiLoop]) + np.min(rand_masses[iiLoop])
#     Merger_kwargs["q"] = np.max(rand_masses[iiLoop])/np.min(rand_masses[iiLoop])
#     Merger_kwargs["chi1"] = rand_spins[iiLoop][0]
#     Merger_kwargs["chi2"] = rand_spins[iiLoop][1]
#     Merger_kwargs["inc"] = rand_angles[iiLoop][0]
#     Merger_kwargs["phi"] = rand_angles[iiLoop][1]
#     Merger_kwargs["lambda"] = rand_angles[iiLoop][2]
#     Merger_kwargs["beta"] = rand_angles[iiLoop][3]
#     Merger_kwargs["psi"] = rand_angles[iiLoop][4]
#
#     # Initialize the source object
#     Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
#     # Can call everything at once by handing arguments to RunInference function directly
#     DoPlots = False
#     OutFileName = "Randomized_angles_spins_0cut_"+str(iiLoop+1)
#     try:
#         print("Trying inference on system", str(iiLoop+1))
#         print("M:", Merger_kwargs["M"], "q:", Merger_kwargs["q"], "chi1:", Merger_kwargs["chi1"], "chi2:", Merger_kwargs["chi2"], ", inc:", Merger_kwargs["inc"], ", phi:", Merger_kwargs["phi"], ", lambda:", Merger_kwargs["lambda"], ", beta:", Merger_kwargs["beta"], ", psi:", Merger_kwargs["psi"])
#         # SYU.RunInference(Merger, LISA, inference_params=inference_params, Plots=DoPlots, OutFileName=OutFileName,JsonFileAndPath=None,**RunTimekwargs)
#     except:
#         print("Error in inference on system", str(iiLoop+1), "... Params:")
#         print("M:", Merger_kwargs["M"], "q:", Merger_kwargs["q"], "chi1:", Merger_kwargs["chi1"], "chi2:", Merger_kwargs["chi2"], ", inc:", Merger_kwargs["inc"], ", phi:", Merger_kwargs["phi"], ", lambda:", Merger_kwargs["lambda"], ", beta:", Merger_kwargs["beta"], ", psi:", Merger_kwargs["psi"])
#         print("Skipping to next system...")




















# ########################### Example - Run ptmcee randomized angles, spins, mass ratio ###########################
#
# # Functions to randomize spins and angles
# def draw_random_angles(size=1):
#     inc = np.arccos(np.random.uniform(low=-1., high=1., size=size))
#     phi = np.random.uniform(low=-np.pi, high=np.pi, size=size)
#     lambd = np.random.uniform(low=-np.pi, high=np.pi, size=size)
#     beta = np.arcsin(np.random.uniform(low=-1., high=1., size=size))
#     psi = np.random.uniform(low=0., high=np.pi, size=size)
#     return np.array([inc, phi, lambd, beta, psi]).T
# def draw_random_spins(low=0., high=1., size=1):
#     chi1 = np.random.uniform(low=low, high=high, size=size)
#     chi2 = np.random.uniform(low=low, high=high, size=size)
#     return np.array([chi1, chi2]).T
# def draw_random_massratio(low=np.log10(0.1), high=np.log10(1.), size=1):
#     q = 10.**(np.random.uniform(low=low, high=high, size=size))
#     return np.array([q]).T
#
# # Draw the random values
# n = 20
# rand_spins = draw_random_spins(size=n)
# rand_angles = draw_random_angles(size=n)
# rand_massratios = draw_random_massratio(size=n)
#
# # Create options to Initialize the detector object
# LISA_base_kwargs = {"TDI":'TDIAET'}
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_base_kwargs)
#
# # Set the time to merger cut if looking at pre-merger only - 'None' to ignore
# OneHour = 60.*60.
# T_obs_end_to_mergers = [-30.*24.*OneHour, -7.*24.*OneHour, -3.*24.*OneHour, -24.*OneHour, -10.*OneHour, -5.*OneHour, -OneHour, -60., None]
# T_obs_labels = ["1mon", "1wk", "3d", "1d", "10hr", "5hr", "1hr", "1min", "0cut"]
#
# # Create options to Initialize the source object
# # NB Dec = beta and Ra = lambda
# Merger_base_kwargs = {"q": 1.1, "M": 3e6, "z": 1., "chi1": 0.9,
#         "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3.,
#         "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM',
#         "Lframe":True, "DeltatL_cut":None}
#
# # Create dictionary of inference params along with their prior ranges and prior types
# # Note that there are default prior ranges for the anglular parameters: chi1, chi2, inc, phi, lambda, beta, psi
# inference_params = {
#   "infer_params": ["chi1", "chi2", "dist", "inc", "phi", "lambda", "beta", "psi"],
#   "params_range": [[-1., 1.], [-1., 1.], [5e3, 2e5], [0.,np.pi], [-np.pi, np.pi], [-np.pi, np.pi], [-np.pi/2.,np.pi/2.], [0.,np.pi]],
#   "prior_type": ["uniform", "uniform", "uniform", "sin", "uniform", "uniform", "cos", "uniform"],
#   "wrap_params": [False, False, False, False, True, True, False, True]
# }
#
# # Create a dictionary of additional parameters to run the inference with (the 'run_params' and 'waveform_params' options)
# # See SYNEX_PTMC.py for default dictionary that you can use as a base for what parameters are modifiable. Can also add a key for
# # any kwarg in source or detector classes and this will be modified in the run, but not updated in the class parameter list.
# RunTimekwargs = {"print_info": True, ### Run param options
#                 "n_walkers": 96, # must be greater than or equal to twice the inference cube dimension
#                 "n_iter": 8000,
#                 "burn_in": 5000, # Throw away at least 2/3 of it
#                 "autocor_method": "acor",
#                 "thin_samples": True, # for speed set this to False
#                 "TDI": "TDIAET", ### waveform param options. These are taken from the source and detector classes first, and then overridden here if specified
#                 "multimodal": True,
#                 "multimodal_pattern": "8modes",
#                 "p_jump": 0.5,
#                 "init_method": "fisher",
#                 "skip_fisher": False,
#                 "n_temps": 10}
#
# # Loop over the variable values. We can either run with the same class of detector and source
# # but using the new variables in the last optional dict, or we can recreate the classes on each iteration
# # of the loop. Since Alexis mentioned he wants to create a class on every iteration of the final product,
# # here we will do the same.
# for iiLoop in range(n):
#     for iiCut in range(len(T_obs_end_to_mergers)):
#     # Copy then update the source class kwargs
#         Merger_kwargs = copy.deepcopy(Merger_base_kwargs)
#         Merger_kwargs["DeltatL_cut"] = T_obs_end_to_mergers[iiCut]
#         Merger_kwargs["q"] = rand_massratios[iiLoop][0]
#         Merger_kwargs["chi1"] = rand_spins[iiLoop][0]
#         Merger_kwargs["chi2"] = rand_spins[iiLoop][1]
#         Merger_kwargs["inc"] = rand_angles[iiLoop][0]
#         Merger_kwargs["phi"] = rand_angles[iiLoop][1]
#         Merger_kwargs["lambda"] = rand_angles[iiLoop][2]
#         Merger_kwargs["beta"] = rand_angles[iiLoop][3]
#         Merger_kwargs["psi"] = rand_angles[iiLoop][4]
#
#         # Initialize the source object
#         Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
#         # Can call everything at once by handing arguments to RunInference function directly
#         DoPlots = False
#         OutFileName = "Randomized_angles_spins_MRat_" + str(iiLoop+1) + "_" + T_obs_labels[iiCut]
#         try:
#             print("Trying inference on system", str(iiLoop+1))
#             print("M:", Merger_kwargs["M"], "q:", Merger_kwargs["q"], "chi1:", Merger_kwargs["chi1"], "chi2:", Merger_kwargs["chi2"], ", inc:", Merger_kwargs["inc"], ", phi:", Merger_kwargs["phi"], ", lambda:", Merger_kwargs["lambda"], ", beta:", Merger_kwargs["beta"], ", psi:", Merger_kwargs["psi"])
#             # SYU.RunInference(Merger, LISA, inference_params=inference_params, Plots=DoPlots, OutFileName=OutFileName,JsonFileAndPath=None,**RunTimekwargs)
#         except:
#             print("Error in inference on system", str(iiLoop+1), "... Params:")
#             print("M:", Merger_kwargs["M"], "q:", Merger_kwargs["q"], "chi1:", Merger_kwargs["chi1"], "chi2:", Merger_kwargs["chi2"], ", inc:", Merger_kwargs["inc"], ", phi:", Merger_kwargs["phi"], ", lambda:", Merger_kwargs["lambda"], ", beta:", Merger_kwargs["beta"], ", psi:", Merger_kwargs["psi"])
#             print("Skipping to next system...")
















# ########################### Example - Run ptmcee inference with lisabeta ###########################
#
# # Create options to Initialize the detector object
# LISA_kwargs = {"TDI":'TDIAET'}
#
# # Set the cut time for time to merger
# T_obs_end_to_merger
#
# # Create options to Initialize the source object
# # NB Dec = beta and Ra = lambda
# Merger_kwargs = {"q": 10., "M": 6e6, "z": 1., "chi1": 0.,
#         "chi2": 0., "beta" : 0., "lambda" : 0.,
#         "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM',
#         "Lframe":True, "DeltatL_cut":T_obs_end_to_merger}
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_kwargs)
#
# # Initialize the source object
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
# # Create dictionary of inference params along with their prior ranges and prior types
# # Note that there are default prior ranges for the anglular parameters: chi1, chi2, inc, phi, lambda, beta, psi
# inference_params = {
#   "infer_params": ["chi1", "chi2", "dist", "inc", "phi", "lambda", "beta", "psi"],
#   "params_range": [[-1., 1.], [-1., 1.], [5e3, 2e5], [0.,np.pi], [-np.pi, np.pi], [-np.pi, np.pi], [-np.pi/2.,np.pi/2.], [0.,np.pi]],
#   "prior_type": ["uniform", "uniform", "uniform", "sin", "uniform", "uniform", "cos", "uniform"],
#   "wrap_params": [False, False, False, False, True, True, False, True]
# }
#
# # Create a dictionary of additional parameters to run the inference with (the 'run_params' and 'waveform_params' options)
# # See SYNEX_PTMC.py for default dictionary that you can use as a base for what parameters are modifiable. Can also add a key for
# # any kwarg in source or detector classes and this will be modified in the run, but not updated in the class parameter list.
# RunTimekwargs = {"print_info": True, ### Run param options
#                 "n_walkers": 96, # must be greater than or equal to twice the inference cube dimension
#                 "n_iter": 8000,
#                 "burn_in": 5000, # Throw away at least 2/3 of it
#                 "autocor_method": "acor",
#                 "thin_samples": True, # for speed set this to False
#                 "TDI": "TDIAET", ### waveform param options. These are taken from the source and detector classes first, and then overridden here if specified
#                 "multimodal": True,
#                 "multimodal_pattern": "8modes",
#                 "p_jump": 0.5,
#                 "init_method": "fisher",
#                 "skip_fisher": False,
#                 "n_temps": 10}
#
# # Can call everything at once by handing arguments to RunInference function directly
# DoPlots = False
# OutFileName = "RedGridTest_Mpi_0cut_4"
# SYU.RunInference(Merger, LISA, inference_params=inference_params, Plots=DoPlots, OutFileName=OutFileName,JsonFileAndPath=None,**RunTimekwargs)
#








# # ########################### Example - Run ptmcee inference with lisabeta ###########################
#
# # Create options to Initialize the detector object
# LISA_kwargs = {"TDI":'TDIAET'}
#
# # Create options to Initialize the source object
# # NB Dec = beta and Ra = lambda
# Merger_kwargs = {"q": 1.1, "M": 5.e7, "z": 3.1, "chi1": 0.9,
#    "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3., "inc": np.pi/10., "psi": -0.4,  "approximant" : 'IMRPhenomD' }
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_kwargs)
#
# # Initialize the source object
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
# # Create dictionary of inference params along with their prior ranges and prior types
# # Note that there are default prior ranges for the anglular parameters: chi1, chi2, inc, phi, lambda, beta, psi
# inference_params = {
#   "infer_params": ["chi1", "chi2", "dist", "inc", "phi", "lambda", "beta", "psi"],
#   "params_range": [[-0.5, 1.], [-0.5, 1.], [5.0e3, 2.0e5], [0.,np.pi], [0.,1.], [-np.pi, np.pi], [-np.pi/2.,np.pi/2.], [0.,1.]],
#   "prior_type": ["uniform", "uniform", "uniform", "sin", "uniform", "uniform", "cos", "uniform"],
#   "wrap_params": [False, False, False, False, True, True, False, True]
# }
#
# # Create a dictionary of additional parameters to run the inference with (the 'run_params' and 'waveform_params' options)
# # See SYNEX_PTMC.py for default dictionary that you can use as a base for what parameters are modifiable. Can also add a key for
# # any kwarg in source or detector classes and this will be modified in the run, but not updated in the class parameter list.
# InferenceTechkwargs = {"print_info": False, ### Run param options
#                 "n_walkers": 16, # must be greater than or equal to twice the inference cube dimension
#                 "n_iter": 40000,
#                 "burn_in": 20000, # This is taken off n_iter, NOT dont in addition to n_iter.
#                 "autocor_method": "acor",
#                 "thin_samples": True, # for speed set this to False
#                 "TDI": "TDIAET", ### waveform param options. These are taken from the source and detector classes first, and then overridden here if specified
#                 }
#
# # Param dict of which parameters we cant to calculate our FoC over
# FoM_loop_kwargs = {"inc":[0.,np.pi,20]}
#
# # Can call everything at once by handing arguments to RunInference function directly
# RunGrid = False
# FigureOfMerit="SkyLocInfer"
# t_start = time.time()
# FoM_array = SYU.RunFoMOverRange(Merger,LISA,FoM_loop_kwargs,FigureOfMerit,RunGrid,inference_params,**InferenceTechkwargs)
# t_end = time.time()
# print("Time to complete RunFoMOverRange: ", round(10.*(t_end-t_start))/10., "s")
#
# # Save the data incase things break later...
# FoMFileName = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/RangeInference_inc/FoM_array.json"
# a_file = open(FoMFileName, "w")
# with open(FoMFileName, 'w') as f:
#     json.dump(FoM_array, f, indent=2)
# f.close()
# print("Finished data-dump to file: ", FoMFileName)
#
# # Open example data file and plot
# FilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/RangeInference_inc/InjValues_inc_0.0_m1_26190476.2_m2_23809523.8_chi1_0.9_chi2_1.0_Deltat_0.0_dist_27064.8_phi_0.0_lambda_1.0_beta_-1.2_psi_-0.4.h5"
# SYU.PlotInferenceData(FilePath, False)
# FilePath = FilePath[:-3]+"_raw.h5"
# SYU.PlotInferenceData(FilePath, False)
#
# # Open data file of FoM info
# a_file = open(FoMFileName, "r")
# FoM_array = json.load(a_file)
# a_file.close()
#
# # Each entry to the FoM array is a dictionary of useful stuff
# # let's plot a difference between injected and most probable values:
# beta_mode_inj_diff = []
# beta_mode_inj_LowError = []
# beta_mode_inj_HighError = []
# lamda_mode_inj_diff = []
# lamda_mode_inj_LowError = []
# lamda_mode_inj_HighError = []
#
# # for PosteriorValDict in FoM_array["inc_ys"]:
# # Mean error is stdv, median error are Q1 and Q3, mode erros is FWHM.
# for ii in range(len(FoM_array["inc_ys"])):
#     for key in FoM_array["inc_ys"][ii]:
#         mean = FoM_array["inc_ys"][ii][key]["mean"]
#         mode = FoM_array["inc_ys"][ii][key]["mode"]
#         median = FoM_array["inc_ys"][ii][key]["median"]
#         StDv = FoM_array["inc_ys"][ii][key]["StDv"]
#         # FWHMs = FoM_array["inc_ys"][ii][key]["FWHMs"]
#         Q1 = FoM_array["inc_ys"][ii][key]["LowerQuart"]
#         Q3 = FoM_array["inc_ys"][ii][key]["UpperQuart"]
#         x_inj = FoM_array["inc_ys"][ii][key]["injected_value_Lframe"] # Lframe is right.
#         if key == "beta": # Some sort of percentage error thing
#             beta_mode_inj_diff.append(100.*(x_inj-mode)/x_inj)
#             beta_mode_inj_LowError.append(100.*abs(x_inj-Q1)/x_inj)
#             beta_mode_inj_HighError.append(100.*abs(x_inj-Q3)/x_inj)
#         elif key == "lambda":
#             lamda_mode_inj_diff.append(100.*(x_inj-mode)/x_inj)
#             lamda_mode_inj_LowError.append(100.*abs(x_inj-Q1)/x_inj)
#             lamda_mode_inj_HighError.append(100.*abs(x_inj-Q3)/x_inj)
#
# LoopedParam = FoM_array["inc_xs"]
#
# plt.figure()
# plt.errorbar(LoopedParam, beta_mode_inj_diff, [beta_mode_inj_LowError, beta_mode_inj_HighError], linestyle="None", label="beta")
# plt.errorbar(LoopedParam, lamda_mode_inj_diff, [lamda_mode_inj_LowError, lamda_mode_inj_HighError], linestyle="None", label="lambda")
# plt.grid()
# plt.legend()
# plt.show()








# ########################### Example - Simple Get waveform ###########################
#
# # Set the total time of the signal
# # T = 31558149.763545576 # 1 year in seconds
# # toffset = 0. # units of seconds ### Didnt change anything when modified
# # timetomerger_max = 1. # Time to merger (yr) to set an fstart cutting the waveform ### Changes fstart (larger for longer times)
# # DeltatL_cut = None # Total time of observations in seconds from t0 ### Didnt change anything when modified
# # t0 = (4.*60.*60)/T # 0. # Start time of observations in years
# # tmin = 0. # Starting time of observation window (yr)
# # tmax = 0.1 # (T-4.*60.*60)/T # Ending time of observation window (yr)
# T_obs_end_to_merger = -7.*24.*60.*60 # Cut waveform four hours before merger
#
# # Create options to Initialize the detector object
# LISA_kwargs = {"TDI":'TDIAET'}
#
# # NB Dec = beta and Ra = lambda
# # Merger_kwargs = {"q": 1.1, "M": 6e6, "z": 3., "chi1": 0.9,
# #         "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3.,
# #         "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM',
# #         "DeltatL_cut":T_obs_end_to_merger}
# Merger_kwargs = {"q": 1.1, "M": 2e8, "z": 3., "chi1": 0.9,
#         "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3.,
#         "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM',
#         "DeltatL_cut":T_obs_end_to_merger} # ,
#         # "timetomerger_max":2.5, "minf":1e-6}
#
# # Initialize the detector object
# LISA = SYDs.LISA(**LISA_kwargs)
#
# # Initialize the source object
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
# # Set the maxf property based on what time is needed to orient telescopes
# # SYU.fmaxFromTimeToMerger(Merger, LISA, T_obs_end_to_merger=T_obs_end_to_merger)
#
# # convert to params dictionaries
# [param_dict, waveform_params, _ ] = SYU.ClassesToParams(Merger,LISA,"Inference")
#
# # Get the detection data
# tdifreqseries_base = SYU.GetSMBHGWDetection_FreqSeries(Merger, LISA)
#
# # Output stuff into useful forms - 22 mode`
# modes = tdifreqseries_base['modes']
#
# # Compute noises
# SYU.ComputeDetectorNoise(Merger, LISA)
#
# # Extract the noise for plotting
# Snvals = LISA.Snvals
#
# # Create the container for the full waveform
# freqs = ['linear', None]
# freqs = SYU.GenerateFreqs(freqs, param_dict, **waveform_params)
# h_full = {}
# for chan in Snvals.keys(): # channels:
#     h_full[chan] = np.zeros_like(freqs, dtype=complex)
#     for lm in modes:
#         h_full[chan] += tdifreqseries_base[lm][chan]
#
# # Add higher harmonics to get the full waveform
# h_full_tot = 0.
# Snvals_tot = 0.
# for chan in Snvals.keys():
#     # Totals
#     h_full_tot += h_full[chan]
#     Snvals_tot += Snvals[chan]
#
# # Find the amp
# h_full_tot_amp = np.sqrt(h_full_tot*np.conj(h_full_tot))*np.sqrt(freqs)
# h_full_tot_amp = abs(h_full_tot)*np.sqrt(freqs)
#
# # Plot waveform
# plt.loglog(freqs, h_full_tot_amp)
#
# # Plot noise
# plt.loglog(freqs, np.sqrt(Snvals_tot))
#
# # labels
# plt.xlabel(r"Frequency [Hz]")
# plt.ylabel(r"$\tilde{\mathrm{h}}_{\mathrm{full}}$ [Hz$^{-1/2}$]")
#
# # Ranges
# # plt.ylim([1.e-23, 1.e-18])
#
# # Show and grid
# plt.show()
# plt.grid()
