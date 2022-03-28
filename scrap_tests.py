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




########################### Example - Load saved objects after tiling and plot useful results ###########################

# Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_dev.dat"}
# Merger=SYSs.SMBH_Merger(**Merger_kwargs)

DetSaveFiles1=sorted(glob.glob("/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/ExpTime_doSingleExposure/Athena_base_exposuretime_*.dat"))
DetSaveFiles2=sorted(glob.glob("/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/ExpTime_Not_doSingleExposure/Athena_base_exposuretime_*.dat"))
DetSaveFiles3=sorted(glob.glob("/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/ExpTime_Overhead_doSingleExposure/Athena_base_exposuretime_*.dat"))
DetSaveFiles4=sorted(glob.glob("/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/ExpTime_Overhead_Not_doSingleExposure/Athena_base_exposuretime_*.dat"))
detectors1=[SYDs.Athena(**{"ExistentialFileName":DetSaveFile}) for DetSaveFile in DetSaveFiles1]
detectors2=[SYDs.Athena(**{"ExistentialFileName":DetSaveFile}) for DetSaveFile in DetSaveFiles2]
detectors3=[SYDs.Athena(**{"ExistentialFileName":DetSaveFile}) for DetSaveFile in DetSaveFiles3]
detectors4=[SYDs.Athena(**{"ExistentialFileName":DetSaveFile}) for DetSaveFile in DetSaveFiles4]
detectors=[detectors1,detectors2,detectors3,detectors4]

labels=[r"No overhead, doSingleExposure=True",r"No overhead, doSingleExposure=False",
        r"Overhead, doSingleExposure=True",r"Overhead, doSingleExposure=False"]

# SYU.PlotPhotonAccumulation(detectors, SaveFig=False, SaveFileName=None)
SYU.PlotSourcePhotons(detectors, labels=labels, SaveFig=False, SaveFileName=None)


























#########################################################################################################################
