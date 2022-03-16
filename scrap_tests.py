import SYNEX.SYNEX_Utils as SYU
import SYNEX.SYNEX_Detectors as SYDs
import SYNEX.SYNEX_Sources as SYSs
import numpy as np
import os
import lisabeta.utils.plotutils as plotutils
import time
import glob
from astropy.cosmology import Planck13, z_at_value # needed only to convert a given distance to redshift at initialization
from astropy.cosmology import WMAP9 as cosmo
from astropy.time import Time
import json
import healpy as hp
import gwemopt
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
Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_base.dat"} # ,
                 # "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_dev.dat"}

# Base telescope args
t0 = '2034-01-01T00:00:00.00' # YYYY-MM-DDTHH:mm:SS.MS 01/01/2034
t = Time(t0, format='isot', scale='utc').gps
Athena_kwargs={"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_base.dat", # "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_dev.dat",
                "orbitFile":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/orbit_files/Athena_20340601_728d_inc60_R750Mkm_ecc4_ArgPeri20_AscNode-10_phi020_P90_frozenFalse_base.dat",
                "NeworbitFile":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/orbit_files/Athena_20340601_728d_inc60_R750Mkm_ecc4_ArgPeri20_AscNode-10_phi020_P90_frozenFalse_dev.dat",
                "telescope":"Athena",
                "tilesType" : "moc", # "greedy", # "hierarchical", # "ranked", # moc/greedy/hierarchical/ranked/galaxy.
                "timeallocationType" : "powerlaw", # "absmag" / "powerlaw" / "waw" / "manual" / "pem"
                "scheduleType" : "greedy",
                "doCalcTiles" : False, # Calculates the no. of "hierarchical" tiles based on sky area prob threashold and FoV
                "Ntiles" : None, # 50, # Speific to tilesType=hierarchical and greedy. Needs to be set if "doCalcTiles"=False
                "frozenAthena" : False, # False,
                "exposuretime" : 10000., # 6*60*60, # 60., #
                "min_observability_duration" : 10./3600., # in HOURS
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
                "slew_rate" : 100., # 1., # in deg/s -- Idea paper has 1 deg/sec
                "horizon" : 0., # 30.,                    ### None if we want a telesscopic orbit?
                "doMinimalTiling" : True, #  True,
                "readout" : 0.0001, # 6
                "doSingleExposure" : True, # False
                "iterativeOverlap" : 1.0,
                "maximumOverlap" : 1.0,
               } ### What about doPerturbativeTiling? ### "doPerturbativeTiling" : True
# Athena_kwargs["ExistentialFileName"]="/home/baird/SYNEX/Saved_Telescope_Dicts/Athena_base.dat"
# Athena_kwargs["NewExistentialFileName"]="/home/baird/SYNEX/Saved_Telescope_Dicts/Athena_dev.dat"
# Athena_kwargs["orbitFile"]="/home/baird/SYNEX/orbit_files/Athena_20340601_728d_inc60_R750Mkm_ecc4_ArgPeri20_AscNode-10_phi020_P90_frozenFalse_base.dat"
# Athena_kwargs["NeworbitFile"]="/home/baird/SYNEX/orbit_files/Athena_20340601_728d_inc60_R750Mkm_ecc4_ArgPeri20_AscNode-10_phi020_P90_frozenFalse_dev.dat"

# Test tiling with detector cloning
ex_times=np.logspace(1,4,num=2,endpoint=True,base=10.)
cloning_params={"exposuretime":ex_times}
tiling_t0=time.time()
detectors = SYU.TileSkyArea(Merger_kwargs,detectors=None,base_telescope_params=Athena_kwargs,cloning_params=cloning_params)
tiling_t1=time.time()
print("Total time for",len(detectors),"detectors:",tiling_t1-tiling_t0, "s")

for detector in detectors:
    T0_mjd=detector.detector_source_coverage["Start time (mjd)"]
    Xs,Ys=[],[]
    for ExT0s,ExTs in zip(detector.detector_source_coverage["Source tile start times (s)"],detector.detector_source_coverage["Source tile exposuretimes (s)"]): Xs+=[ExT0s/86400.,(ExT0s+ExTs)/86400.]
    for CumCounts in np.cumsum(detector.detector_source_coverage["Source photon counts"]): Ys+=[CumCounts,CumCounts]
    print(detector,Xs,Ys,detector.detector_source_coverage["Source tile exposuretimes (s)"])






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











# ######################### Example - Plot sky loc area histogram for 1d data for paper plot ###########################
#
# JsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/Randomized_SYNEX/RemakePaperPlots/"
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/Randomized_SYNEX/RemakePaperPlots/"
#
# CutsToPlot = ["1d"]
# TotAreas_json_contents = {}
#
# # for iCut in range(len(CutsToPlot)):
# #     JsonFilesToSearch = JsonFilePath + "Randomized_angles_spins_MRat_*_" + CutsToPlot[iCut]
# #     JsonFiles=glob.glob(JsonFilesToSearch+".json")
# #     # JsonFiles = [JsonFilePath + "0cut/Randomized_SYNEX_0cut_7.json"]
# #     TotAreas = []
# #     loop_count = 0
# #     print(len(JsonFiles), "systems found...")
# #     for JsonFile in JsonFiles:
# #         file_ID = JsonFile.split("_")[-2]
# #         DataFileLocAndName = DataFilePath + (JsonFile.split('/')[-1]).split('.')[0] + ".h5"
# #         try:
# #             TotArea = SYU.GetTotSkyAreaFromPostData(DataFileLocAndName,ConfLevel=0.9)
# #             TotAreas.append(TotArea)
# #             loop_count += 1
# #             TotAreas_json_contents[file_ID] = {"ID": file_ID, "Tot Sky Area": TotArea,
# #                                                "CL":0.9, "json file": JsonFile, # Maybe CL can be the real CL for each file... Each is a little above 90%.
# #                                                "data file": DataFileLocAndName}
# #         except:
# #             print("File load failed - doesn't exist.")
# #     bins = np.logspace(np.log10(6.),np.log10(6e4), 40) # 20
# #     TotAreas = [(TotArea*(180./np.pi)**2) for TotArea in TotAreas]
# #     plt.hist(TotAreas, bins=bins, label=CutsToPlot[iCut])
# # plt.legend()
# # plt.xlabel(r"$\log_{10}(\Omega)\,$[$\log_{10}(\mathrm{deg}^2)$]")
# # plt.ylabel(r"Count")
# # ax = plt.gca()
# # ax.set_xscale('log')
# # plt.show()
#
# tile_hierarchy_file = DataFilePath.split("SYNEX")[0] + "SYNEX/DataBase_1d_90CL_TotalSkyAreas.json"
# f = open(tile_hierarchy_file)
# TotAreas_json_contents = json.load(f)
# f.close()
#
# test = TotAreas_json_contents["1"].get("Tot Sky Area")
# print(test)
# with open(tile_hierarchy_file, 'w') as f:
#     json.dump(TotAreas_json_contents, f, indent=2)
# f.close()













# ########################## Example - Plot sky loc area histograms for each cut time on one graph ###########################
#
# JsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/Randomized_SYNEX/RemakePaperPlots/"
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/Randomized_SYNEX/RemakePaperPlots/"
#
# CutsToPlot = "0cut" # ["0cut","4hr","1d","1wk"]
# SystemsToPlot = [ii for ii in range(1,7)]
#
# # JsonFilesToSearch = JsonFilePath + "Randomized_angles_spins_MRat_"
# # JsonFiles=glob.glob(JsonFilesToSearch+"*" + CutsToPlot + ".json")
#
# JsonFilesToSearch = JsonFilePath + "Randomized_angles_spins_MRat_1_"
# JsonFiles=glob.glob(JsonFilesToSearch+"*"+".json")
#
# TotAreas = []
# for JsonFile in JsonFiles:
#     DataFileLocAndName = DataFilePath + (JsonFile.split('/')[-1]).split('.')[0] + ".h5"
#     TotArea = SYU.GetTotSkyAreaFromPostData(DataFileLocAndName,ConfLevel=0.90,bins=150)
#     print(TotArea*(180./np.pi)**2)
#     SYU.PlotInferenceLambdaBeta(DataFileLocAndName, bins=150, SkyProjection=False, SaveFig=False)
#
#












# ########################### Example - Plot inference results points over grid of Modified Bayes Factors ###########################
# json_file = "FoMOverRange_SkyModeLikelihoodsLambda_Inc_Beta.json" # "FoMOverRange_SkyModeLikelihoodsLambda_inc_maxf.json" # "FoMOverRange_SkyModeLikelihoodsBeta_M_z.json" #
# BF_lim = 20. # 4. # 100. # for Lambda
# SaveFig = True
# InlayType="histogram"
# ModeJumpLimit = 0.1
# SYU.PlotBayesFactorFoMFromJsonWithAllInferencePoints(json_file, BF_lim=BF_lim, ModeJumpLimit=ModeJumpLimit, SaveFig=SaveFig)
# SYU.PlotBayesFactorFoMFromJsonWithInferenceInlays(json_file, BF_lim=BF_lim, SaveFig=SaveFig, InlayType=InlayType)
# DataFileLocAndName = "inference_data/IncGridTest_NoMpi_3PiBy6.h5"
# SYU.GetPosteriorStats(DataFileLocAndName, ModeJumpLimit=ModeJumpLimit)

# GridDataFile = "inference_data/MassGridTest_NoMpi_3e7.h5" # "inference_data/maxfGridTest_NoMpi_4eM4.h5" # "inference_data/MassGridTest_NoMpi_5e7.h5"
# SYU.PlotInferenceLambdaBeta(GridDataFile, SkyProjection=False, SaveFig=False)
# SYU.PlotHistsLambdaBeta(GridDataFile, SaveFig=False,ParamToPlot=["beta","lambda"])
# ParamToPlot=["beta", "lambda"]
#
# import glob
# GridDataFiles = glob.glob("inference_data/RedGridTest_NoMpi_*.h5")
#
# for GridDataFile in GridDataFiles:
#     print(GridDataFile[-6:-3])
#     # SYU.PlotHistsLambdaBeta(GridDataFile, ParamToPlot=ParamToPlot)
#     SYU.PlotInferenceLambdaBeta(GridDataFile)




# t1 = time.time()
# # Otions to Initialize detector
# LISA_kwargs = {"TDI":'TDIAET'}
#
# # Otions to Initialize source
# Merger_kwargs = {"q": 1.1, "M": 5130000., "z": 3.11, "chi1": 0.9,
#         "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3.,
#          "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM'}
#
# # Initialize detector object
# LISA = SYD.LISA(**LISA_kwargs)
#
# # Initialize source object
# Merger = SYS.SMBH_Merger(**Merger_kwargs)
#
# # Get likelihood dictionary for sky modes
# lnL_skymodes,params_skymode = SYU.GetSkyMultiModeProbFromClasses(Merger, LISA)
# t2 = time.time()

# import glob
# MassGridJsons = glob.glob("/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/MassGridTest_NoMpi_*.json")
# RedGridJsons = glob.glob("/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/RedGridTest_NoMpi_*.json")
#
# t1 = time.time()
#
# Colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'] # ['b', 'g', 'r', 'c', 'm', 'y', 'k']
#
# plt.figure()
# IncLabels = True
# for MassGridFile in MassGridJsons:
#     lnL_skymodes,params_skymode = SYU.GetSkyMultiModeProbFromJson(MassGridFile)
#     Colori = 0
#     for key,value in lnL_skymodes.items():
#         if IncLabels:
#             if key == (1,0):
#                 label = str(key) + ", inj"
#             else:
#                 label = key
#             plt.plot(params_skymode[key]["M"], value, 'x', c=Colors[Colori], label=label)
#         else:
#             plt.plot(params_skymode[key]["M"], value, 'x', c=Colors[Colori])
#         Colori += 1
#     IncLabels = False
# ax = plt.gca()
# ax.set_xscale('log')
# plt.xlabel(r"M$_{\mathrm{tot}} \; (M_{\odot})$")
# plt.ylabel(r"$\ln(\mathcal{L})$")
# plt.legend(ncol=2)
# plt.ylim([-50., 1.])
# plt.show()
# plt.grid()

# plt.figure()
# IncLabels = True
# for RedGridFile in RedGridJsons:
#     lnL_skymodes,params_skymode = SYU.GetSkyMultiModeProbFromJson(RedGridFile)
#     Colori = 0
#     for key,value in lnL_skymodes.items():
#         if IncLabels:
#             if key == (1,0):
#                 label = str(key) + ", inj"
#             else:
#                 label = key
#             plt.plot(params_skymode[key]["dist"], value, 'x', c=Colors[Colori], label=label)
#             # plt.plot(z_at_value(Planck13.distmod, params_skymode[key]["dist"]), value, 'x', c=Colors[Colori], label=label)
#         else:
#             plt.plot(params_skymode[key]["dist"], value, 'x', c=Colors[Colori])
#             # plt.plot(z_at_value(Planck13.distmod, params_skymode[key]["dist"]), value, 'x', c=Colors[Colori])
#         Colori += 1
#     IncLabels = False
# ax = plt.gca()
# plt.xlabel(r"$\mathcal{D} \; $(Mpc)")
# plt.ylabel(r"$\ln(\mathcal{L})$")
# plt.legend(ncol=2)
# plt.ylim([-50., 1.])
# plt.show()
# plt.grid()
#
# t2 = time.time()
#
# print("Time to execute:", round((t2-t1)*1000.)/1000., "s")
#
# for key,value in lnL_skymodes.items():
#     print(key, value, params_skymode[key]["inc"], params_skymode[key]["lambda"], params_skymode[key]["beta"], params_skymode[key]["psi"])


# MtotTestedStrs = ["2e6", "3e6", "4e6", "5e6", "6e6", "7e6", "8e6", "9e6", "1e7", "2e7", "3e7", "4e7", "5e7", "6e7", "7e7", "8e7", "9e7", "1e8", "2e8"]
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/MassGridTest_NoMpi_"
# jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/MassGridTest_NoMpi_"
#
# PostVals = {}
# for MtotStr in MtotTestedStrs:
#     print(MtotStr)
#     DataFile = DataFilePath + MtotStr + ".h5" # "Mass.h5"
#     jsonFile = jsonFilePath + MtotStr + ".json" # "Mass.json"
#     DataFileRaw = DataFile[:-3]+"_raw.h5"
#     # posteriors_full = plotutils.load_params_posterior_lisa_smbh(jsonFile, DataFile, format='multiemcee', load_fisher=True)
#     # fig = plotutils.corner_plot(posteriors_full['injparams_Lframe'], posteriors_full['post'], add_posteriors=None, output=False, output_dir=None, output_file=None, histograms=True, fisher=False, fishercov=posteriors_full['fishercov'], Lframe=True, params=['chi1', 'chi2', 'dist', 'inc', 'phi', 'lambda', 'beta', 'psi'], color=plotutils.plotpalette[0], cov_color='k', add_colors=[plotutils.plotpalette[1]], show_truths=True, truth_color='k', bins=25, show_histograms=False, plot_datapoints=False);
#     # plt.show()
#     # SYU.PlotInferenceData(DataFile, SaveFig=False)
#     SYU.PlotInferenceLambdaBeta(DataFile, SaveFig=False)
#     PostVals[MtotStr] = SYU.GetPosteriorVal(DataFile)

# # Same again for redshift
# RedTestedStrs = ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
# DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/RedGridTest_NoMpi_"
# jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/RedGridTest_NoMpi_"
#
# PostVals = {}
# for red in RedTestedStrs:
#     print(red)
#     DataFile = DataFilePath + red + ".h5" # "Mass.h5"
#     jsonFile = jsonFilePath + red + ".json" # "Mass.json"
#     DataFileRaw = DataFile[:-3]+"_raw.h5"
#     # posteriors_full = plotutils.load_params_posterior_lisa_smbh(jsonFile, DataFile, format='multiemcee', load_fisher=True)
#     # fig = plotutils.corner_plot(posteriors_full['injparams_Lframe'], posteriors_full['post'], add_posteriors=None, output=False, output_dir=None, output_file=None, histograms=True, fisher=False, fishercov=posteriors_full['fishercov'], Lframe=True, params=['chi1', 'chi2', 'dist', 'inc', 'phi', 'lambda', 'beta', 'psi'], color=plotutils.plotpalette[0], cov_color='k', add_colors=[plotutils.plotpalette[1]], show_truths=True, truth_color='k', bins=25, show_histograms=False, plot_datapoints=False);
#     # plt.show()
#     # SYU.PlotInferenceData(DataFile, SaveFig=False)
#     SYU.PlotInferenceLambdaBeta(DataFile, SaveFig=False)
#     # PostVals[red] = SYU.GetPosteriorVal(DataFile)

# # Create y entries for plotting
# Lambdas = []
# dLambdas_lower = []
# dLambdas_Upper = []
# betas = []
# dbetas_lower = []
# dbetas_Upper = []
# InjLambda = []
# InjBeta = []
# plt.figure()
# ax = plt.gca()
# InjValLam = PostVals[MtotStr]["lambda"]["injected_value_Lframe"]
# InjValBeta = PostVals[MtotStr]["beta"]["injected_value_Lframe"]
# ax.axhline(InjValLam, linestyle=":", color="black", label=r"$\lambda_{L,inj}$")
# ax.axhline(InjValBeta, linestyle=":", color="blue", label=r"$\beta_{L,inj}$")
# LegendCount = 0
# for MtotStr in MtotTestedStrs:
#     # Lambdas
#     ys = PostVals[MtotStr]["lambda"]["mode"]
#     LowErr = [Low1-Low2 for Low1,Low2 in zip(ys,PostVals[MtotStr]["lambda"]["LowerFWHMs"])]
#     HighErr = [Low2-Low1 for Low1,Low2 in zip(ys,PostVals[MtotStr]["lambda"]["UpperFWHMs"])]
#     yErrors = [LowErr, HighErr]
#     print(ys)
#     print(yErrors)
#     # InjVal = PostVals[MtotStr]["lambda"]["injected_value_Lframe"]
#     xs = [str(MtotStr+r"$\,M_{t}:\,\lambda_"+str(ii+1))+"$" for ii in range(len(ys))]
#     plt.errorbar(xs, ys, yErrors, linestyle="None", marker='+', color='black', label=r"$\lambda_L$")
#     # plt.plot(xs, InjVal*np.ones(np.shape(ys)), linestyle="None", marker='x', markeredgecolor='black', label=r"$lambda_{L,inj}$")
#
#     # Betas
#     ys = PostVals[MtotStr]["beta"]["mode"]
#     LowErr = [Low1-Low2 for Low1,Low2 in zip(ys,PostVals[MtotStr]["beta"]["LowerFWHMs"])]
#     HighErr = [Low2-Low1 for Low1,Low2 in zip(ys,PostVals[MtotStr]["beta"]["UpperFWHMs"])]
#     yErrors = [LowErr, HighErr]
#     # InjVal = PostVals[MtotStr]["beta"]["injected_value_Lframe"]
#     xs = [str(MtotStr+r"$\,M_{t}:\,\beta_"+str(ii+1))+"$" for ii in range(len(ys))]
#     plt.errorbar(xs, ys, yErrors, linestyle="None", marker='+', color='blue', label=r"$\beta_L$")
#     # plt.plot(xs, InjVal*np.ones(np.shape(ys)), linestyle="None", marker='x', markeredgecolor='blue', label=r"$beta_{L,inj}$")
#
#     # Push legend only after first set of points plotted, otherwise too many legend entries
#     if LegendCount == 0:
#         plt.legend()
#         LegendCount += 1
#
# plt.ylabel(r'$\lambda_L \, \mathrm{or} \, \beta_L \; (\mathrm{deg.})$')
# plt.show()
# plt.grid()





### LISABETA has equivalent plotting tools... Why don't you just use this?
# posteriors_full = plotutils.load_params_posterior_lisa_smbh(jsonFilePath, FilePath, format='multiemcee', load_fisher=True)
#
# fig = plotutils.corner_plot(posteriors_full['injparams_Lframe'], posteriors_full['post'], add_posteriors=None, output=False, output_dir=None, output_file=None, histograms=True, fisher=False, fishercov=posteriors_full['fishercov'], Lframe=True, params=['chi1', 'chi2', 'dist', 'inc', 'phi', 'lambda', 'beta', 'psi'], color=plotutils.plotpalette[0], cov_color='k', add_colors=[plotutils.plotpalette[1]], show_truths=True, truth_color='k', bins=25, show_histograms=False, plot_datapoints=False);
# plt.show()
