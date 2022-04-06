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








FileName = "Randomized_SYNEX2/Randomized_angles_spins_MRat_9_10hr"
JsonFileLocAndName,H5FileLocAndName=SYU.CompleteLisabetaDataAndJsonFileNames(FileName)
# Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/Randomized_angles_spins_MRat_1_1mon.dat"}
# Merger = SYU.GetSourceFromLisabetaData(FileName,**Merger_kwargs)
# SYU.PlotSkyMapData(Merger)
posteriors_full = plotutils.load_params_posterior_lisa_smbh(JsonFileLocAndName, H5FileLocAndName, format='multiemcee', load_fisher=True)
fig = plotutils.corner_plot(posteriors_full['injparams_Lframe'], posteriors_full['post'], add_posteriors=None, output=False, output_dir=None, output_file=None, histograms=True, fisher=False, fishercov=posteriors_full['fishercov'], Lframe=True, params=['chi1', 'chi2', 'dist', 'inc', 'phi', 'lambda', 'beta', 'psi'], color=plotutils.plotpalette[0], cov_color='k', add_colors=[plotutils.plotpalette[1]], show_truths=True, truth_color='k', bins=25, show_histograms=False, plot_datapoints=False);
ax_list = fig.axes
for ii in range(len(ax_list)):
    if "ylabel" in list([ax_list[ii].tick_params()]): ax_list[ii].set_ylabel(ax_list[ii].get_ylabel(),fontsize="x-small")
    if "xlabel" in list([ax_list[ii].tick_params()]): ax_list[ii].set_xlabel(ax_list[ii].get_xlabel(),fontsize="x-small")
plt.show()









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
