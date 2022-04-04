import numpy as np
try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
except ModuleNotFoundError:
    plt = None
    mpl = None

import astropy.constants as const
import time
import astropy.units as u

import lisabeta.lisa.lisatools as lisatools
# import lisabeta.utils.plotutils as plotutils

from SYNEX import SYNEX_Detectors as SYDs
from SYNEX import SYNEX_Sources as SYSs
from SYNEX import SYNEX_Utils as SYU

from numpy.random import rand

import time
import json
import glob
import copy

from astropy.cosmology import WMAP9 as cosmo



########################### Run ptmcee randomized angles, spins, mass ratio ###########################

# Functions to randomize spins and angles
def draw_random_angles(size=1):
    inc = np.arccos(np.random.uniform(low=-1., high=1., size=size))
    phi = np.random.uniform(low=-np.pi, high=np.pi, size=size)
    lambd = np.random.uniform(low=-np.pi, high=np.pi, size=size)
    beta = np.arcsin(np.random.uniform(low=-1., high=1., size=size))
    psi = np.random.uniform(low=0., high=np.pi, size=size)
    return np.array([inc, phi, lambd, beta, psi]).T
def draw_random_spins(low=0., high=1., size=1):
    chi1 = np.random.uniform(low=low, high=high, size=size)
    chi2 = np.random.uniform(low=low, high=high, size=size)
    return np.array([chi1, chi2]).T
def draw_random_massratio(low=np.log10(0.1), high=np.log10(1.), size=1):
    q = 10.**(np.random.uniform(low=low, high=high, size=size))
    return np.array([q]).T

# Draw the random values
n = 400
rand_spins = draw_random_spins(size=n)
rand_angles = draw_random_angles(size=n)
rand_massratios = draw_random_massratio(size=n)

# Create options to Initialize the detector object
LISA_base_kwargs = {"TDI":'TDIAET'}

# Initialize the detector object
LISA = SYDs.LISA(**LISA_base_kwargs)

# Set the time to merger cut if looking at pre-merger only - 'None' to ignore
OneHour = 60.*60.
T_obs_end_to_mergers = [-30.*24.*OneHour, -3.*7.*24.*OneHour, -2.*7.*24.*OneHour, -7.*24.*OneHour, -3.*24.*OneHour, -24.*OneHour, -10.*OneHour, -5.*OneHour, None] # [-24.*OneHour] # [None] #
T_obs_labels = ["1mon", "3wk", "2wk", "1wk", "3d", "1d", "10hr", "5hr", "0cut"] # ["1d"] # ["0cut"] #

# Create options to Initialize the source object
# NB Dec = beta and Ra = lambda
Merger_base_kwargs = {"q": 1.1, "M": 8e6, "z": 1.3, "chi1": 0.9,
        "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3.,
        "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM',
        "Lframe":True, "DeltatL_cut":None,
        "ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/TestSystem_9d_base.dat",
        "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/TestSystem_9d_tmp.dat"}

# Create dictionary of inference params along with their prior ranges and prior types
# Note that there are default prior ranges for the anglular parameters: chi1, chi2, inc, phi, lambda, beta, psi
inference_params = {
  "infer_params": ["chi1", "chi2", "dist", "inc", "phi", "lambda", "beta", "psi"],
  "params_range": [[-1., 1.], [-1., 1.], [5e3, 2e5], [0.,np.pi], [-np.pi, np.pi], [-np.pi, np.pi], [-np.pi/2.,np.pi/2.], [0.,np.pi]],
  "prior_type": ["uniform", "uniform", "uniform", "sin", "uniform", "uniform", "cos", "uniform"],
  "wrap_params": None
}

# Create a dictionary of additional parameters to run the inference with (the 'run_params' and 'waveform_params' options)
# See SYNEX_PTMC.py for default dictionary that you can use as a base for what parameters are modifiable. Can also add a key for
# any kwarg in source or detector classes and this will be modified in the run, but not updated in the class parameter list.
RunTimekwargs = {"print_info": True, ### Run param options
                "n_walkers": 96, # 16, # must be greater than or equal to twice the inference cube dimension
                "n_iter": 8000, # 400, #
                "burn_in": 5000, #  1, #  Throw away at least 2/3 of it
                "autocor_method": "autocor_new", # "acor",
                "thin_samples": True, # for speed set this to False
                "TDI": "TDIAET", ### waveform param options. These are taken from the source and detector classes first, and then overridden here if specified
                "multimodal": True,
                "multimodal_pattern": "8modes",
                "p_jump": 0.5,
                "init_method": "fisher",
                "skip_fisher": False,
                "n_temps": 10}

# Loop over the variable values. We can either run with the same class of detector and source
# but using the new variables in the last optional dict, or we can recreate the classes on each iteration
# of the loop. Since Alexis mentioned he wants to create a class on every iteration of the final product,
# here we will do the same.
for iiLoop in range(n):
    for iiCut in range(len(T_obs_end_to_mergers)):
        # Copy then update the source class kwargs
        Merger_kwargs = copy.deepcopy(Merger_base_kwargs)
        Merger_kwargs["DeltatL_cut"] = T_obs_end_to_mergers[iiCut]
        Merger_kwargs["q"] = 1./rand_massratios[iiLoop][0]
        Merger_kwargs["chi1"] = rand_spins[iiLoop][0]
        Merger_kwargs["chi2"] = rand_spins[iiLoop][1]
        Merger_kwargs["inc"] = rand_angles[iiLoop][0]
        Merger_kwargs["phi"] = rand_angles[iiLoop][1]
        Merger_kwargs["lambda"] = rand_angles[iiLoop][2]
        Merger_kwargs["beta"] = rand_angles[iiLoop][3]
        Merger_kwargs["psi"] = rand_angles[iiLoop][4]

        # Initialize the source object
        # Merger = SYSs.SMBH_Merger(**Merger_kwargs)

        # Set output file names
        OutFileName = "Randomized_angles_spins_MRat_" + str(iiLoop+1) + "_" + T_obs_labels[iiCut]
        RunTimekwargs["out_file"]=OutFileName

        # Can call everything at once by handing arguments to RunInference function directly
        DoPlots = False
        # try:
        print("Trying inference on system", str(iiLoop+1), T_obs_labels[iiCut])
        SYU.RunInference(Merger_kwargs, LISA, inference_params=inference_params, PlotInference=DoPlots,PlotSkyMap=DoPlots,**RunTimekwargs)
        # except:
        #     print("Error in inference on system", str(iiLoop+1), "... Params:")
        #     print("M:", Merger_kwargs["M"], "q:", Merger_kwargs["q"], "chi1:", Merger_kwargs["chi1"], "chi2:", Merger_kwargs["chi2"], ", inc:", Merger_kwargs["inc"], ", phi:", Merger_kwargs["phi"], ", lambda:", Merger_kwargs["lambda"], ", beta:", Merger_kwargs["beta"], ", psi:", Merger_kwargs["psi"])
        #     print("Skipping to next system...")
