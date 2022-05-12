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
import SYNEX.SYNEX_Utils as SYU
import SYNEX.SYNEX_Detectors as SYDs
import SYNEX.SYNEX_Sources as SYSs
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pylab as pylab
from SYNEX.SYNEX_Utils import pylab_params
pylab.rcParams.update(pylab_params)







### Check we can load a tiled detector correctly ###

FileName = "/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Randomized_SYNEX2/Randomized_angles_spins_MRat_10_2wk.dat"
Athena_kwargs={"ExistentialFileName":FileName}
Athena = SYDs.Athena(**Athena_kwargs)
print("Check 1:", Athena_kwargs.keys(), Athena.detector_source_coverage)

t0 = '2034-01-01T00:00:00.00' # YYYY-MM-DDTHH:mm:SS.MS 01/01/2034
t = Time(t0, format='isot', scale='utc').gps
Athena_kwargs={"ExistentialFileName":FileName,
                "verbose":False,
                "telescope":"Athena",
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
Athena = SYDs.Athena(**Athena_kwargs)
print("Check 2:", Athena_kwargs.keys(), Athena.detector_source_coverage)

##### Both checks now passed.


# FileName = "Randomized_SYNEX2/Randomized_angles_spins_MRat_1_5hr"
# JsonFileLocAndName,H5FileLocAndName=SYU.CompleteLisabetaDataAndJsonFileNames(FileName)
# Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/Randomized_angles_spins_MRat_1_5hr.dat"}
# Merger = SYU.GetSourceFromLisabetaData(FileName,**Merger_kwargs)
# print(JsonFileLocAndName, Merger.JsonFile, -Merger.DeltatL_cut/(24*60*60))
# Merger.CreateSkyMapStruct() #### Everything until here is right... Check this function!
# print(len(Merger.map_struct["prob"]))
# SYU.PlotSkyMapData(Merger)
# posteriors_full = plotutils.load_params_posterior_lisa_smbh(JsonFileLocAndName, H5FileLocAndName, format='multiemcee', load_fisher=True)
# fig = plotutils.corner_plot(posteriors_full['injparams_Lframe'], posteriors_full['post'], add_posteriors=None, output=False, output_dir=None, output_file=None, histograms=True, fisher=False, fishercov=posteriors_full['fishercov'], Lframe=True, params=['chi1', 'chi2', 'dist', 'inc', 'phi', 'lambda', 'beta', 'psi'], color=plotutils.plotpalette[0], cov_color='k', add_colors=[plotutils.plotpalette[1]], show_truths=True, truth_color='k', bins=25, show_histograms=False, plot_datapoints=False);
# ax_list = fig.axes
# for ii in range(len(ax_list)):
#     if "ylabel" in list([ax_list[ii].tick_params()]): ax_list[ii].set_ylabel(ax_list[ii].get_ylabel(),fontsize="x-small")
#     if "xlabel" in list([ax_list[ii].tick_params()]): ax_list[ii].set_xlabel(ax_list[ii].get_xlabel(),fontsize="x-small")
# plt.show()









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
