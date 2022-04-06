"""
TO DO:
-Do we always want to match filenames?
-Do we want to restructure outputs so they are all in one output directory and labeled under event?
Might get tricky if we are e.g. running randomized sources... For the moment keep outputs seperated.
"""
# Import stuff
import sys
import os
import pathlib
# Set the path to root SYNEX directory so file reading and writing is easier -- NOTE: NO TRAILING SLASH
SYNEX_PATH=str(pathlib.Path(__file__).parent.resolve()).split("SYNEX")[0]+"SYNEX"
import astropy.constants as const
import astropy.units as u
from astropy.time import Time
import numpy as np
import healpy as hp
import gwemopt
from gwemopt import utils as gou
from gwemopt import tiles as gots
from gwemopt.plotting import add_edges
import astropy, astroplan
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import Planck13, z_at_value
import astropy.units as u
import copy
from ast import literal_eval
import statistics
from scipy import stats
import time
import json, h5py

# Import lisebeta stuff
import lisabeta
import lisabeta.lisa.lisa_fisher as lisa_fisher
import lisabeta.lisa.lisa as lisa
import lisabeta.lisa.lisatools as lisatools
import lisabeta.lisa.pyLISAnoise as pyLISAnoise
import lisabeta.lisa.pyresponse as pyresponse
import lisabeta.tools.pytools as pytools
import lisabeta.pyconstants as pyconstants
import lisabeta.inference.inference as inference

# Import ptemcee stuff
import SYNEX.SYNEX_PTMC as SYP
import ptemcee
import ptemcee.mpi_pool as mpi_pool

# import detector and source classes
import SYNEX.SYNEX_Sources as SYSs
import SYNEX.SYNEX_Detectors as SYDs
import SYNEX.segments_athena as segs_a

# Plotting stuff
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pylab as pylab
try:
    import ligo.skymap.plot
    cmap = "cylon"
except:
    cmap = 'PuBuGn'
pylab_params = {'legend.fontsize': 4, # 'x-large',
         'axes.labelsize': 4, # 'x-large',
         'axes.titlesize': 4,
         'xtick.labelsize': 4, # 'x-large',
         'ytick.labelsize': 4, # 'x-large'}
         'lines.markersize': 0.7,
         'lines.linewidth': 0.7,
         'font.size': 4} # ,
#          'font.family':'serif'} # 0.1} #### 'fontname'
pylab.rcParams.update(pylab_params)
try:
    mpl.use('MacOSX')
except:
    a=1

# Stop warnings about deprecated methods...
import warnings
warnings.filterwarnings("ignore")

# # Import lalsimulation stuff incase we decide to use it later - you have to download it from the website I think like acor etc.
# import lalsimulation

# mpi stuff
try:
    from mpi4py import MPI
except ModuleNotFoundError:
    MPI = None

# global variable for histogram plots and mode calculations
hist_n_bins=1000






#                     #################################
#                     #                               #
#                     #   General Utility functions   #
#                     #                               #
#                     #################################

def ClassesToParams(source, detector=None, CollectionMethod="Inference",**kwargs):
    """ Function to change parameters embedded in the classes and handed to functions
    in dictionaries, to dictionaries and parameter lists for function calls in lisabeta.

    Parameters
    ---------
    Source : SYNEX source object.

    Detector : SYNEX GW detector object

    NOTE:
    -------
    If source has JsonFile attribute then the "waveform_params" handed back
    will be loaded preferentially from this instead of from GW detector object
    """
    # First gather the parameters into the right format to pass to the lisabeta functions
    param_dict = {
              "m1": source.m1,                # Redshifted mass of body 1 (solar masses)
              "m2": source.m2,                # Redshifted mass of body 2 (solar masses)
              "chi1": source.chi1,            # Dimensionless spin of body 1 (in [-1, 1])
              "chi2": source.chi2,            # Dimensionless spin of body 2 (in [-1, 1])
              "Deltat": source.Deltat,        # Frequency at the start of the observations (Hz)
              "dist": source.dist,            # Luminosity distance (Mpc)
              "inc": source.inc,              # Inclination angle (rad)
              "phi": source.phi,              # Observer's azimuthal phase (rad)
              "lambda": source.lamda,         # Source longitude in SSB-frame (rad)
              "beta": source.beta,            # Source latitude in SSB-frame (rad)
              "psi": source.psi,              # Polarization angle (rad)
              "Lframe": source.Lframe,                 # Params always in Lframe - this is easier for the sampler.
              }

    # Take out any modification, to these parameters in tge kwargs dict
    for key, value in kwargs.items():
        if key in param_dict: param_dict[key] = value

    # Sort out case where z is specified in kwargs
    if "z" in kwargs:
        param_dict["dist"] = cosmo.luminosity_distance(kwargs["z"]).to("Mpc").value

    # Sort out case where q and/or M are specified in kwargs
    if "M" in kwargs and "q" in kwargs:
        param_dict["m1"] = kwargs["M"]*(kwargs["q"](1.+kwargs["q"]))
        param_dict["m2"] = kwargs["M"]/(1.+kwargs["q"])
    elif "M" in kwargs and "m1" in kwargs:
        param_dict["m2"] = kwargs["M"]-kwargs["m1"]
    elif "M" in kwargs and "m2" in kwargs:
        param_dict["m1"] = kwargs["M"]-kwargs["m2"]
    elif "q" in kwargs and "m1" in kwargs:
        param_dict["m2"] = kwargs["m1"]/kwargs["q"]
    elif "q" in kwargs and "m2" in kwargs:
        param_dict["m1"] = kwargs["q"]*kwargs["m2"]
    elif "q" in kwargs: # Now cases where only one is specified - take other needed values from source
        # vary q and keep M the same
        param_dict["m1"] = source.M*(kwargs["q"]/(1.+kwargs["q"]))
        param_dict["m2"] = source.M/(1.+kwargs["q"])
    elif "M" in kwargs:
        # vary M and keep q the same
        param_dict["m1"] = kwargs["M"]*source.q/(1.+source.q)
        param_dict["m2"] = kwargs["M"]/(1.+source.q)

    # Unwrap the positional arguments and override the source and detector values set
    # source params
    # elif source.JsonFile and os.path.isfile(source.JsonFile):
    #     with open(source.JsonFile, 'r') as f: input_params=json.load(f)
    #     waveform_params=input_params["waveform_params"]
    if hasattr(source,"waveform_params"):
        waveform_params=source.waveform_params
    else:
        waveform_params = {}
        if source.timetomerger_max!=None:
            waveform_params["timetomerger_max"] = source.timetomerger_max                        # Included
        if source.minf!=None:
            waveform_params["minf"] = source.minf                        # Included
        if source.maxf!=None and not source.DeltatL_cut:
            waveform_params["maxf"] = source.maxf                        # Included
        if source.t0!=None:
            waveform_params["t0"] = source.t0                        # Included
        if source.tref!=None:
            waveform_params["tref"] = source.tref                        # Included
        if source.phiref!=None:
            waveform_params["phiref"] = source.phiref                        # Included
        if source.fref_for_phiref!=None:
            waveform_params["fref_for_phiref"] = source.fref_for_phiref                        # Included
        if source.fref_for_tref!=None:
            waveform_params["fref_for_tref"] = source.fref_for_tref                        # Included
        if source.force_phiref_fref!=None:
            waveform_params["force_phiref_fref"] = source.force_phiref_fref                        # Included
        if source.toffset!=None:
            waveform_params["toffset"] = source.toffset                        # Included
        if source.modes!=None:
            waveform_params["modes"] = source.modes                        # Included
        if source.acc!=None:
            waveform_params["acc"] = source.acc                        # Included
        if source.approximant!=None:
            waveform_params["approximant"] = source.approximant                        # Included
        if hasattr(source, "DeltatL_cut"):
            waveform_params["DeltatL_cut"] = source.DeltatL_cut                        # Included
        # detector params
        if hasattr(detector, "tmin") and detector.tmin!=None:
            waveform_params["tmin"] = detector.tmin                        # Included
        if hasattr(detector, "tmax") and detector.tmax!=None:
            waveform_params["tmax"] = detector.tmax                        # Included
        if hasattr(detector, "TDI") and detector.TDI!=None:
            waveform_params["TDI"] = detector.TDI                        # Included
        if hasattr(detector, "order_fresnel_stencil") and detector.order_fresnel_stencil!=None:
            waveform_params["order_fresnel_stencil"] = detector.order_fresnel_stencil                        # Included
        if hasattr(detector, "LISAconst") and detector.LISAconst!=None:
            waveform_params["LISAconst"] = detector.LISAconst                        # Included
        if hasattr(detector, "responseapprox") and detector.responseapprox!=None:
            waveform_params["responseapprox"] = detector.responseapprox                        # Included
        if hasattr(detector, "frozenLISA") and detector.frozenLISA!=None:
            waveform_params["frozenLISA"] = detector.frozenLISA                        # Included
        if hasattr(detector, "TDIrescaled") and detector.TDIrescaled!=None:
            waveform_params["TDIrescaled"] = detector.TDIrescaled
        if hasattr(detector, "LISAnoise") and detector.LISAnoise!=None:
            waveform_params["LISAnoise"] = detector.LISAnoise                        # Included

    # update waveform param given in kwarks dict
    for key,val in kwargs.items():
        if key in waveform_params: waveform_params[key]=val

    # Optional parameters that depend on what method is being called
    # Need to understand if these are important differences...
    if CollectionMethod=="Inference":
        if "fend" in kwargs:
            waveform_params["fend"] = kwargs["fend"]
        elif source.fend!=None:
            waveform_params["fend"] = source.fend
    elif CollectionMethod=="Fisher":
        if "tf_method" in kwargs:
            waveform_params["tf_method"] = kwargs["tf_method"]
        elif hasattr(source, "tf_method") and source.tf_method!=None:
            waveform_params["tf_method"] = source.tf_method
        if "fstart22" in kwargs:
            waveform_params["fstart22"] = kwargs["fstart22"]
        elif hasattr(source, "fstart22") and source.fstart22!=None:
            waveform_params["fstart22"] = source.fstart22
        if "fend22" in kwargs:
            waveform_params["fend22"] = kwargs["fend22"]
        elif hasattr(source, "fend22") and source.fend22!=None:
            waveform_params["fend22"] = source.fend22
        if "gridfreq" in kwargs:
            waveform_params["gridfreq"] = kwargs["gridfreq"]
        elif hasattr(source, "gridfreq") and source.gridfreq!=None:
            waveform_params["gridfreq"] = source.gridfreq

    # Extra parameters important for fisher only (?)
    extra_params={'scale_mode_freq': True, 'use_buggy_LAL_tpeak': False}

    return param_dict, waveform_params, extra_params

def ParamsToClasses(input_params,CollectionMethod="Inference",**kwargs):
    """ Function to return SYNEX classes from dictionary of parameters.

    Parameters
    ----------
    input_params : Dictionary
            Params to give to lisabeta inference methods. Usually saved by SYNEX
            in Json format in file SYNEX/inference_param_files/...

    CollectionMethod : String
            Method of formatting dictionaries depending on what 'input_params' was
            used for. Options are 'Inference' or 'Fisher' depending on if the
            output from the lisabeta launch was a Fisher localization estimate
            or a full MCMC inference (with ptemcee inside lisabeta).

    kwargs : Dictionary
            Optional dictionary of values to adjust the input params. This can be
            useful if we want many copies of the same objects with a subset of
            parameters changed. In this case we would call this function many times
            with the same input parameters but with **kwargs containing the adjusted
            param value at each call iteration.


    Return
    ----------
    Source : SYNEX source class
            Source class object created wth attributed given in input_params and
            modified by anything in kwargs, using formatting according to 'CollectionMethod'.


    TO DO
    ----------
        -Include detector class in return
        -Include EM classes like Athena by adding options to 'CollectionMethod'
        -Include EM variables in kwargs and input_params
        -Create file with default source and detector params like with gwemopt but for lisabeta objects
    """
    # Form the three sub dictionaries
    run_params = input_params["run_params"]
    waveform_params = input_params["waveform_params"]
    prior_params = input_params["prior_params"]
    param_dict = input_params["source_params"]

    # Hand everything useful to one dict - anything not used in class creation is left alone.
    source_kwargs = {**param_dict, "lisabetaFile":run_params["out_dir"] + run_params["out_name"] + ".h5", **waveform_params}

    # Check if we will modify something
    MUTATED=False
    if len(kwargs.keys())>0:
        # more than just 'ExistentialFileName' specified
        ValueCheck = [value!=source_kwargs[key] for key,value in kwargs.items() if key in source_kwargs]
        if any([ValueCheck]):
            # Values are changed so tell source init it's mutated
            MUTATED = True
    source_kwargs.update(kwargs)
    source_kwargs["MUTATED"]=MUTATED

    # # Optional parameters that depend on what method is being called
    # # Need to understand if these are important differences...
    # if CollectionMethod=="Inference":
    #     if "fend" in kwargs:
    #         source_kwargs["fend"] = kwargs["fend"]
    # elif CollectionMethod=="Fisher":
    #     if "tf_method" in kwargs:
    #         source_kwargs["tf_method"] = kwargs["tf_method"]
    #     if "fstart22" in kwargs:
    #         source_kwargs["fstart22"] = kwargs["fstart22"]
    #     if "fend22" in kwargs:
    #         source_kwargs["fend22"] = kwargs["fend22"]
    #     if "gridfreq" in kwargs:
    #         source_kwargs["gridfreq"] = kwargs["gridfreq"]

    # Now create source object
    source = SYSs.SMBH_Merger(**source_kwargs)

    return source

def CompleteLisabetaDataAndJsonFileNames(FileName):
    """
    Function to give back the json and h5 filenames complete with paths
    based on a single input Filename that could be to either json or h5
    file, with or without full path to file.

    This could be done just by reading in the json file and returning it
    and the saved datafile path inside... Something feels faster about
    avoiding reading in jsons and just manipulating strings though, and
    since we want to apply this code to an optimization routine I am
    starting to avoid too much overhead per loop...
    """

    # Create paths to data and json folders
    LisabetaJsonPath = SYNEX_PATH + "/inference_param_files"
    LisabetaDataPath = SYNEX_PATH + "/inference_data"

    if FileName[0]!="/":
        LisabetaJsonPath+="/"
        LisabetaDataPath+="/"

    # Figure out if its a json or h5 filename
    if FileName[-3]==".":
        JsonFileLocAndName = FileName[:-3] + '.json'
        H5FileLocAndName = FileName
    elif FileName[-5]==".":
        JsonFileLocAndName = FileName
        H5FileLocAndName = FileName[:-5] + '.h5'
    else:
        JsonFileLocAndName = FileName + '.json'
        H5FileLocAndName = FileName + '.h5'

    # Add file path if only the filenames specified
    bool_vec = [len(FileName.split("inference_param_files"))==1,len(FileName.split("inference_data"))==1]
    if bool_vec[0] and bool_vec[1]:
        H5FileLocAndName = LisabetaDataPath + H5FileLocAndName
        JsonFileLocAndName = LisabetaJsonPath + JsonFileLocAndName
    elif bool_vec[0] and not bool_vec[1]:
        JsonFileLocAndName = LisabetaJsonPath + H5FileLocAndName .split("inference_data")[-1]
        JsonFileLocAndName = ".".join(JsonFileLocAndName.split(".")[:-1])+".json"
    elif not bool_vec[0] and bool_vec[1]:
        H5FileLocAndName = LisabetaDataPath + JsonFileLocAndName.split("inference_param_files")[-1]
        H5FileLocAndName = ".".join(H5FileLocAndName.split(".")[:-1])+".h5"

    # Check if the subdirectories exist for both data and json files
    JsonPathOnly="/".join(JsonFileLocAndName.split("/")[:-1])
    DataPathOnly="/".join(H5FileLocAndName.split("/")[:-1])
    try:
        # See if the directories exist in case we load ource from savefile or something
        pathlib.Path(JsonPathOnly).mkdir(parents=True, exist_ok=True)
        pathlib.Path(DataPathOnly).mkdir(parents=True, exist_ok=True)
    except:
        JsonPathOnly=LisabetaJsonPath
        DataPathOnly=LisabetaDataPath
        JsonFileLocAndName=JsonPathOnly+JsonFileLocAndName.split("inference_param_files")[-1]
        H5FileLocAndName=DataPathOnly+H5FileLocAndName.split("inference_data")[-1]
        pathlib.Path(JsonPathOnly).mkdir(parents=True, exist_ok=True)
        pathlib.Path(DataPathOnly).mkdir(parents=True, exist_ok=True)

    # Return full paths and names all harmonious and what not
    return JsonFileLocAndName,H5FileLocAndName

def GWEMOPTPathChecks(go_params, config_struct):
    """
    Function to check formatting of file paths specified at initiation of gwemopt
    params dictionary. This will check the locations as well as file names are
    coherent with telescope type... But haven't coded that yet...

    TO DO: Need to add checks that the proper file locations are specified (tile_files etc)...
    """

    PATH_VARS = ["outputDir", "tilingDir", "catalogDir"] # , "configDirectory"]
    FILE_VARS = ["coverageFiles", "lightcurveFiles"]
    CONFIG_FILE_VARS = ["tesselationFile"]

    # Check if paths are complete
    for PathVar in PATH_VARS:
        if len(go_params[PathVar].split("SYNEX"))==1:
            # Needs proper gwemopt_output or whatever added...
            go_params[PathVar] = SYNEX_PATH+go_params[PathVar]
        # Check if it exists
        try:
            pathlib.Path(go_params[PathVar]).mkdir(parents=True, exist_ok=True)
        except:
            go_params[PathVar] = SYNEX_PATH+"/"+go_params[PathVar].split("/SYNEX/")[-1]
            pathlib.Path(go_params[PathVar]).mkdir(parents=True, exist_ok=True)

    # Now specified files
    for FileVar in FILE_VARS:
        if len(go_params[FileVar].split("SYNEX"))==1:
            # Needs SYNEX_PATH added
            go_params[FileVar] = SYNEX_PATH+go_params[FileVar]
        # Check file endings to match gwemopt hardcoded stuff
        if len(go_params[FileVar].split("."))==1:
            go_params[FileVar] = go_params[FileVar] + ".dat"
        elif go_params[FileVar].split(".")[-1]!="dat":
            go_params[FileVar] = go_params[FileVar].split(".")[0] + ".dat"
        # Check if it exists
        try:
            PathOnly = "/".join(go_params[FileVar].split("/")[:-1])
            pathlib.Path(PathOnly).mkdir(parents=True, exist_ok=True)
        except:
            go_params[FileVar] = SYNEX_PATH+"/"+go_params[FileVar].split("/SYNEX/")[-1]
            PathOnly = "/".join(go_params[FileVar].split("/")[:-1])
            pathlib.Path(PathOnly).mkdir(parents=True, exist_ok=True)

    # Now config_struct paths
    for FileVar in CONFIG_FILE_VARS:
        if len(config_struct[FileVar].split("SYNEX"))==1:
            # Needs SYNEX_PATH added
            config_struct[FileVar] = SYNEX_PATH+config_struct[FileVar]
        # Check file endings to match gwemopt hardcoded stuff
        if len(config_struct[FileVar].split("."))==1:
            config_struct[FileVar] = config_struct[FileVar] + ".tess"
        elif config_struct[FileVar].split(".")[-1]!="tess":
            config_struct[FileVar] = ".".join(config_struct[FileVar].split(".")[:-1]) + ".tess"
        # Check if it exists
        try:
            PathOnly = "/".join(config_struct[FileVar].split("/")[:-1])
            pathlib.Path(PathOnly).mkdir(parents=True, exist_ok=True)
        except:
            config_struct[FileVar] = SYNEX_PATH+"/"+config_struct[FileVar].split("/SYNEX/")[-1]
            PathOnly = "/".join(config_struct[FileVar].split("/")[:-1])
            pathlib.Path(PathOnly).mkdir(parents=True, exist_ok=True)

    return go_params, config_struct






#                ###########################################
#                #                                         #
#                #   SYNEX -- Lisabeta Utility functions   #
#                #                                         #
#                ###########################################

def GetPosteriorStats(FileName, ModeJumpLimit=None, LogLikeJumpLimit=None):
    """
    Function to extract a posterior estimate, with errors, for a list of parameters infered
    using ptemcee.
    """

    if not ModeJumpLimit and not LogLikeJumpLimit:
        raise print("Specified neither ModeJumpLimit nor LogLikeJumpLimit- taking only likelihood limit of 20 for posterior peak significance.")
        ModeJumpLimit = None
        LogLikeJumpLimit = 20.
        LoopLimit = -LogLikeJumpLimit
    if not ModeJumpLimit and not LogLikeJumpLimit:
        raise print("Specified both ModeJumpLimit and LogLikeJumpLimit- taking only likelihood limit for posterior peak significance.")
        ModeJumpLimit = None
        LoopLimit = ModeJumpLimit

    if LogLikeJumpLimit:
        LoopLimit = -LogLikeJumpLimit
    elif ModeJumpLimit:
        LoopLimit = ModeJumpLimit

    from scipy.signal import chirp, find_peaks, peak_widths

    [infer_params, inj_param_vals, static_params, meta_data] = read_h5py_file(FileName)
    # print("metadata keys: ", meta_data.keys())
    ndim = len(infer_params.keys())
    infer_param_keys = infer_params.keys()

    # Import utils for mode and FWHM estimations
    import scipy.optimize as opt

    # Grab data for infered parameters
    PosteriorStats = {}
    for key,values in infer_params.items():

        # First histogram for mode locations
        histn,histbins = np.histogram(values, bins=hist_n_bins)
        dhistbins = histbins[2]-histbins[1]
        histbins = histbins[:-1] + 0.5*(dhistbins) # change from left bins edges to middle of bins

        # convert to lists so pop works later
        histn = list(histn)
        histbins = list(histbins)

        # mode list
        Modes = []
        LowerFWHMs = []
        UpperFWHMs = []

        # Initalize error
        LoopError = 1.

        # Bin population must be within a limit of primary mode bin population to be counted as a mode
        while LoopError>LoopLimit:

            # Check if the histogram is still populated - for low limits this can happen.
            if not histn:
                break

            # Get properties of mode
            ModeLocationLogical = histn==max(histn)
            ModeLocationIndex = np.where(ModeLocationLogical)[0][0] # take smallest if more than one - second [0]
            mode = histbins[ModeLocationIndex]
            mode_population = histn[ModeLocationIndex]

            if key == "lambda" and False:
                print("Mode:", mode, ", Mode pop.:", mode_population)
                plt.figure()
                plt.plot(histbins, histn)
                plt.show()

            # Check first if population is sufficient to count as a mode
            if len(Modes)>0:
                if ModeJumpLimit:
                    LoopError = mode_population/primary_mode_population
                    if LoopError < LoopLimit:
                        break
                elif LogLikeJumpLimit:
                    Upper_bool_lnlikes = values<(mode+dhistbins)
                    Lower_bool_lnlikes = values>(mode-dhistbins)
                    bool_lnlikes = [Upper_bool_lnlike and Lower_bool_lnlike for Upper_bool_lnlike, Lower_bool_lnlike in zip(Upper_bool_lnlikes,Lower_bool_lnlikes)]
                    mode_lnlike = statistics.mean(meta_data["lnlike"][bool_lnlikes])
                    LoopError = - primary_lnlike + mode_lnlike # this should be always negative ?
                    if LoopError < LoopLimit: # If less negative than the limit then the new mode has likelihood comparable with the primary mode
                    # Needed to put the extra negatives in so the while loop statement doesn't contradict things
                        break
            else:
                primary_mode_population = mode_population
                Upper_bool_lnlikes = values<(mode+dhistbins)
                Lower_bool_lnlikes = values>(mode-dhistbins)
                bool_lnlikes = [Upper_bool_lnlike and Lower_bool_lnlike for Upper_bool_lnlike, Lower_bool_lnlike in zip(Upper_bool_lnlikes,Lower_bool_lnlikes)]
                primary_lnlike = statistics.mean(meta_data["lnlike"][bool_lnlikes])

            # Scan to find where bin populations fall to 1% or less of mode population
            # Cutting out more than enough incase the data is noisy.
            LowerHalfBins = np.array(histbins[:ModeLocationIndex])
            LowerHalfn = histn[:ModeLocationIndex]
            LowerLim = LowerHalfBins[LowerHalfn<=0.002*mode_population]
            UpperHalfBins = np.array(histbins[ModeLocationIndex:])
            UpperHalfn = histn[ModeLocationIndex:]
            UpperLim = UpperHalfBins[UpperHalfn<=0.002*mode_population]

            # Take the closest values to the maw in case of local mins and other sus behaviour
            if len(LowerLim)==0:
                LowerLim = histbins[0]
            else:
                LowerLim = LowerLim[-1]
            if len(UpperLim)==0:
                UpperLim = histbins[-1]
            else:
                UpperLim = UpperLim[0]

            # pop data out around mode
            IndicesToPop = np.where(np.array(histbins<=UpperLim)&np.array(histbins>=LowerLim))[0]
            histn_popped = histn[IndicesToPop[0]:IndicesToPop[-1]]
            histbins_popped = histbins[IndicesToPop[0]:IndicesToPop[-1]]

            # Recalculate the histogram around the current mode
            values_popped = values[np.where(np.array(values<=UpperLim)&np.array(values>=LowerLim))[0][:]]
            histn_popped,histbins_popped = np.histogram(values_popped, bins=150)
            histbins_popped = histbins_popped[:-1] + 0.5*(histbins_popped[2]-histbins_popped[1])

            # Calculate an average mode using bins around mode within 5% of mode bin population
            mode = np.sum(histbins_popped[histn_popped>=0.95*max(histn_popped)])/(1.*len(histbins_popped[histn_popped>=0.95*max(histn_popped)]))
            mode_population = np.sum(histn_popped[histn_popped>=0.95*max(histn_popped)])/(1.*len(histn_popped[histn_popped>=0.95*max(histn_popped)]))

            # repeat this at 50% (pm 5%) of mode bin population to get FWHM
            LowerInds = np.where(np.array(histn_popped>=0.45*mode_population)&np.array(histn_popped<=0.55*mode_population)&np.array(histbins_popped<mode))[0]
            UpperInds = np.where(np.array(histn_popped>=0.45*mode_population)&np.array(histn_popped<=0.55*mode_population)&np.array(histbins_popped>mode))[0]
            histbins_LowerFWHM = [histbins_popped[LowerInd] for LowerInd in LowerInds]
            histn_LowerFWHM = [histn_popped[LowerInd] for LowerInd in LowerInds]
            if len(histbins_LowerFWHM)==0:
                LowerFWHM = histbins_popped[0]
            else:
                LowerFWHM = np.sum(histbins_LowerFWHM)/(1.*len(histbins_LowerFWHM))
            histbins_UpperFWHM = [histbins_popped[UpperInd] for UpperInd in UpperInds]
            if len(histbins_UpperFWHM)==0:
                UpperFWHM = histbins_popped[-1]
            else:
                UpperFWHM = np.sum(histbins_UpperFWHM)/(1.*len(histbins_UpperFWHM))

            if key == "lambda" and False:
                print(key, ", eval over popped range ", str(mode), str(mode_population))
                plt.figure()
                plt.plot(histbins_popped, histn_popped, 'r')
                plt.axvline(x=LowerFWHM)
                plt.axvline(x=mode)
                plt.axvline(x=UpperFWHM)
                plt.show()

            # continue appending after check for error
            Modes.append(mode)
            LowerFWHMs.append(LowerFWHM)
            UpperFWHMs.append(UpperFWHM)

            # pop out values around modes to update histogram bins and populations
            for x in reversed(IndicesToPop): # relative index changes on each pop so have to remove them in reverse order
                histn.pop(x)
                histbins.pop(x)

        PosteriorStats[key] = {
                         # Averages of data
                         "mean": statistics.mean(list(values)),
                         "median": statistics.median(list(values)),
                         "mode": Modes,
                         # Errors of averages
                         "StDv": np.sqrt(np.sum((statistics.mean(list(values))-list(values))**2)/(len(list(values))*1.-1.)), # fit_stdev,
                         "LowerQuart": np.quantile(list(values),0.25),
                         "UpperQuart": np.quantile(list(values),0.75),
                         "LowerFWHMs": LowerFWHMs,
                         "UpperFWHMs": UpperFWHMs,
                         # Injected values
                         "injected_value_SSBframe": inj_param_vals["source_params_SSBframe"][key][0],
                         "injected_value_Lframe": inj_param_vals["source_params_Lframe"][key][0],
                         } # Note inj value is npfloat64,like all other variables except FWHMs which is not list. This is json serializable.
    return PosteriorStats

def GetTotSkyAreaFromPostData(FileName,ConfLevel=0.9):
    """
    Counting function to give the total sky area for a confidence level given some
    posterior data
    """
    # Unpack data
    [infer_params, inj_param_vals, static_params, meta_data] = read_h5py_file(FileName)
    if not inj_param_vals: # Raw files at first didn't record this, so make sure it's there...
        # Get the inj values from the processed data file instead
        DataFileLocAndName_NotRaw = DataFileLocAndName[:-7] + ".h5"
        [_,inj_param_vals,_,_] = read_h5py_file(DataFileLocAndName_NotRaw)
    labels = ["lambda","beta"]
    if np.size(infer_params[labels[0]][0])>1:
        nsamples = len(infer_params["lambda"][0])
    else:
        nsamples = len(infer_params["lambda"])
    print("Posterior sample length: " + str(nsamples))
    # Grab data for sky location
    data = np.empty([2,nsamples])
    SampleModes = []
    for ii in range(2):
        data[ii][:] = infer_params[labels[ii]]
        histn,histbins = np.histogram(data[ii,:], bins=hist_n_bins)
        histbins = histbins[:-1] + 0.5*(histbins[2]-histbins[1]) # change from left bins edges to middle of bins
        mode = histbins[histn==max(histn)]
        if len(mode)>1:
            mode = mode[1] # Doe now take the first on the list, but if there are several we need to work out what to do there...
        SampleModes.append(mode)
    data = np.transpose(np.array(data)) # should have shape [nsamples, 2]

    # Get injected values
    InjParam_InjVals = []
    for key in labels:
        InjParam_InjVals.append(inj_param_vals["source_params_Lframe"][key][0]) # Lframe is right.

    # Manually do 2D histogram because I don't trust the ones I found online
    bin_max = max(max(data[:,0]),SampleModes[0]+np.pi/10000.)
    bin_min = min(min(data[:,0]),SampleModes[0]-np.pi/10000.)
    if bin_max>np.pi:
        bin_max=np.pi
    if bin_min<-np.pi:
        bin_max=-np.pi
    if isinstance(bin_min,np.ndarray):
        bin_min = bin_min[0]
    if isinstance(bin_max,np.ndarray):
        bin_max = bin_max[0]
    if bin_max-bin_min>3.*np.pi/2.:
        print("Changing lambda bins")
        nbins_lambda = 200
    else:
        nbins_lambda = 100
    lambda_bins = np.linspace(bin_min, bin_max, nbins_lambda+1) # np.linspace(np.min(data[:,5]), np.max(data[:,5]), bins+1)
    bin_max = max(max(data[:,1]),SampleModes[1]+np.pi/20000.)
    bin_min = min(min(data[:,1]),SampleModes[1]-np.pi/20000.)
    if bin_max>np.pi/2.:
        bin_max=np.pi/2.
    if bin_min<-np.pi/2.:
        bin_max=-np.pi/2.
    if isinstance(bin_min,np.ndarray):
        bin_min = bin_min[0]
    if isinstance(bin_max,np.ndarray):
        bin_max = bin_max[0]
    if bin_max-bin_min>3.*np.pi/4.:
        print("Changing beta bins")
        nbins_beta = 200
    else:
        nbins_beta = 100
    beta_bins = np.linspace(bin_min, bin_max, nbins_beta+1) # np.linspace(np.min(data[:,0]), np.max(data[:,0]), bins+1)
    hist2D_pops, areas = hist_lam_bet(data,lambda_bins,beta_bins)

    # update total area if hist value is above sig level value
    TotArea = 0.
    Count = 0.
    CutLim = ConfLevel*np.sum(hist2D_pops)
    iiloop=0
    while iiloop<=np.size(hist2D_pops) and Count<=CutLim: # Second term just there because when bins are large it takes forever as it counts bins with populations 1... Can we find a workaroundf or this? Just sum the singular stuff? Idk...
        TotArea += sum(areas[hist2D_pops==np.max(hist2D_pops)])
        Count += sum(hist2D_pops[hist2D_pops==np.max(hist2D_pops)])
        areas[hist2D_pops==np.max(hist2D_pops)] = 0.
        hist2D_pops[hist2D_pops==np.max(hist2D_pops)] = 0.
        iiloop+=1
    print("Counted bin total population = ", Count*100.*ConfLevel/CutLim, "% of total population.")
    if Count*100.*ConfLevel/CutLim>98.:
        print(DataFileLocAndName)

    return TotArea

def hist_lam_bet(data,lambda_bins,beta_bins):
    # data is n*2 numpy array :: 0 element is lambda and 1 element is beta
    nbins_lambda = len(lambda_bins)-1
    nbins_beta = len(beta_bins)-1
    lambda_BinW = np.diff(lambda_bins)
    beta_BinW = np.diff(beta_bins)
    hist2D_pops = np.empty([nbins_lambda,nbins_beta])
    areas = np.empty([nbins_lambda,nbins_beta])
    for xii in range(nbins_lambda):
        list_len = list(range(len(data[:,0])))
        if xii == 0:
            values_ii = [valii for valii in list_len if (lambda_bins[xii]<=data[valii,0]<=lambda_bins[xii+1])]
        else:
            values_ii = [valii for valii in list_len if (lambda_bins[xii]<=data[valii,0]<lambda_bins[xii+1])]
        beta_pops, beta_bins = np.histogram(data[values_ii,1], bins=beta_bins)
        for yii in range(nbins_beta):
            hist2D_pops[xii,yii] = beta_pops[yii]
            areas[xii,yii] = lambda_BinW[xii]*beta_BinW[yii]
    return hist2D_pops, areas

def GetOctantLikeRatioAndPostProb(FileName,source=None,detector=None):
    # To do
    # Change the betas to lisa frame
    # Add errors to the numbers (root(n) for probs etc)

    # Check filenames
    JsonFileLocAndName,H5FileLocAndName = CompleteLisabetaDataAndJsonFileNames(FileName)

    # Load the json file to grab source param values
    f = open(JsonFileLocAndName)
    data = json.load(f)
    f.close()
    param_dict = data["source_params"]
    waveform_params = data["waveform_params"]

    likelihoodClass = lisa.LikelihoodLISASMBH(param_dict, **waveform_params)
    lnL_skymodes,params_skymode = lisatools.func_loglikelihood_skymodes(likelihoodClass)

    # list skymodes
    skymodes = [(1,0), (1,1), (1,2), (1,3), (-1,0), (-1,1), (-1,2), (-1,3)]

    # Read h5 data
    [infer_params, inj_param_vals, static_params, meta_data] = read_h5py_file(H5FileLocAndName)
    labels = list(infer_params.keys())
    beta_posteriors = infer_params["beta"]
    lambda_posteriors = infer_params["lambda"]

    OctantPostProbs = {}
    lnLikeRatios = {}
    for skymode in skymodes:
        # Create dictionary of skymode likelihood ratios
        lnLikeRatios[skymode] = lnL_skymodes[skymode]# -lnL_skymodes[(1,0)] # this is always negative # lnL_skymodes[skymode] #

        # get the skymode central coords
        beta_sm = params_skymode[skymode]['beta']
        lambda_sm = params_skymode[skymode]['lambda']

        # Get the declination limits of the octant
        beta_lowerlim = np.min([0.,np.sign(params_skymode[(1,0)]['beta'])*skymode[0]*np.pi]) # This should be lisa frame though so you gotta change this to +pi/3 but need to check orientations first
        beta_upperlim = np.max([0.,np.sign(params_skymode[(1,0)]['beta'])*skymode[0]*np.pi])

        # Get the right accension limits of the octant
        lambda_lowerlim = lambda_sm - np.pi/4.
        lambda_upperlim = lambda_sm + np.pi/4.

        # wrap lambda limits into [-pi,pi]
        if lambda_lowerlim<-np.pi:
            lambda_lowerlim = lambda_lowerlim+2.*np.pi
        if lambda_upperlim>np.pi:
            lambda_upperlim = lambda_upperlim-2.*np.pi

        # Collect the samples in the current octant, take care for limits that have wrapped around
        if lambda_lowerlim<lambda_upperlim:
            octant_betas = [beta_posteriors[ii] for ii in range(len(beta_posteriors)) if lambda_posteriors[ii]>=lambda_lowerlim and lambda_posteriors[ii]<=lambda_upperlim]
        else:
            octant_betas = [beta_posteriors[ii] for ii in range(len(beta_posteriors)) if lambda_posteriors[ii]>=lambda_lowerlim or lambda_posteriors[ii]<=lambda_upperlim]
        if beta_lowerlim<beta_upperlim:
            octant_betas = [p for p in octant_betas if p>=beta_lowerlim and p<=beta_upperlim]
        else:
            octant_betas = [p for p in octant_betas if p>=beta_lowerlim or p<=beta_upperlim]

        # Change this into a probability
        OctantPostProbs[skymode] = np.size(octant_betas)/np.size(beta_posteriors)

    return lnLikeRatios,OctantPostProbs

def read_h5py_file(FileName):

    # Check filenames
    JsonFileLocAndName,H5FileLocAndName = CompleteLisabetaDataAndJsonFileNames(FileName)

    # Load data from h5 file
    f = h5py.File(H5FileLocAndName,'r')

    # Get basic source params to differentiate meta-data later
    from SYNEX.SYNEX_PTMC import list_params as full_list_params

    # Take out the dictionary
    infer_params = {}
    inj_param_vals = {}
    static_params = {}
    meta_data = {}
    list_params = f.keys()
    for key in list_params:
        if key == "fishercov":
            meta_data[key] = f[key]
        elif key == "source_params_Lframe" or key == "source_params_SSBframe":
            inj_param_vals[key] = f[key]
        else:
            param_vals = f[key][:]
            if np.array([np.diff(param_vals)==0]).all():
                static_params[key] = param_vals
            elif key in full_list_params:
                infer_params[key] = param_vals
            else:
                meta_data[key] = param_vals
    # f.close()
    return [infer_params, inj_param_vals, static_params, meta_data]

def RunInference(source_or_kwargs, detector, inference_params, PlotInference=False,PlotSkyMap=False,**RunTimekwargs):
    """
    Function to handle lisabeta inference.
    """
    # Using MPI or not
    use_mpi=False
    is_master = True
    mapper = map
    if MPI is not None:
        MPI_size = MPI.COMM_WORLD.Get_size()
        MPI_rank = MPI.COMM_WORLD.Get_rank()
        comm = MPI.COMM_WORLD
        use_mpi = (MPI_size > 1)
        if use_mpi:
            print("MPI rank/size: %d / %d" % (MPI_rank, MPI_size), flush=True)
            pool = ptemcee.mpi_pool.MPIPool(debug=False)
            is_master = pool.is_master()
            mapper = pool.map
        else:
            print("No MPI", flush=True)
            is_master = True
            mapper = map

    # Init source if we need to
    if is_master:
        # See if we have a source class or kwargs -- only master node since this can get memory heavy
        source=SYSs.SMBH_Merger(**source_or_kwargs) if isinstance(source_or_kwargs,dict) else source_or_kwargs
        # Write params to json file
        print("Creating json file...")
        WriteParamsToJson(source,detector,inference_params,is_master,**RunTimekwargs)
        sourceJsonFile=source.JsonFile
    else:
        sourceJsonFile=None # So we only have master node running later if one source and/or detector given

    # Send source to workers
    if use_mpi:
        comm.Barrier()
        sourceJsonFile = comm.bcast(sourceJsonFile, root=0)

    # Start the run. Data will be saved to the 'inference_data' folder by default
    # All processes must execute the run together. mapper (inside ptemcee) will handle coordination between p's.
    # SYP.RunPTEMCEE(source.JsonFile)
    if is_master: print(" --------------- START PTEMCEE --------------- ")
    print("Does Json exist?",MPI_rank,os.path.isfile(sourceJsonFile))
    t1 = time.time()
    command = "python3 " + SYNEX_PATH + "/lisabeta/lisabeta/inference/ptemcee_smbh.py " + sourceJsonFile
    os.system(command)
    t2 = time.time()
    if is_master: print(" ---------------- END PTEMCEE ---------------- ")
    if is_master: print("Time to execute ptemcee: ", round((t2-t1)*10.)/10., "s")
    if is_master: print(" --------------------------------------------- ")

    # if not use_mpi and is_master:
    #     # Create sky_map struct in source object
    #     source.sky_map = source.H5File.split("inference_data")[0] + 'Skymap_files' + source.H5File.split("inference_data")[-1]
    #     source.sky_map = source.sky_map[:-3] + '.fits'
    #     source.CreateSkyMapStruct()
    #     # Update saved source data
    #     source.ExistentialCrisis()

    # Call plotter if asked for
    if is_master and PlotInference:
        # Automatically save the fig if called within inference.
        PlotInferenceData(source.H5File,SaveFig=True)
    if is_master and PlotSkyMap:
        # Automatically save the fig if called within inference.
        PlotSkyMapData(source.sky_map,SaveFig=True)

def WriteParamsToJson(source, detector, inference_params, IsMaster=True, **RunTimekwargs):
    """
    Function to save source and GW detector params to json file for lisabeta inference runs,
    using the saved jsoon file name in source class.
    """
    # See if we have output filenames to asign
    if "out_file" in RunTimekwargs:
        JsonFile,H5File=CompleteLisabetaDataAndJsonFileNames(RunTimekwargs["out_file"])
        source.JsonFile=JsonFile
        source.H5File=H5File
        del RunTimekwargs["out_file"]

    # double check names and paths are ok
    if source.JsonFile==None and source.H5File==None:
        from datetime import date
        today = date.today()
        d = today.strftime("%d_%m_%Y")
        JsonOrOutFile = d + "_InferParams_" + '_'.join(inference_params["infer_params"])
        JsonFile,H5File=CompleteLisabetaDataAndJsonFileNames(JsonOrOutFile)
        source.JsonFile=JsonFile
        source.H5File=H5File
        print("Creating default json name:",JsonFile)
    elif source.JsonFile==None and source.H5File!=None:
        JsonFile,H5File=CompleteLisabetaDataAndJsonFileNames(source.H5File)
        source.JsonFile=JsonFile
        source.H5File=H5File
    elif source.JsonFile!=None and source.H5File==None:
        JsonFile,H5File=CompleteLisabetaDataAndJsonFileNames(source.JsonFile)
        source.JsonFile=JsonFile
        source.H5File=H5File
    elif source.JsonFile!=None and source.H5File!=None:
        JsonFile,H5File=CompleteLisabetaDataAndJsonFileNames(source.JsonFile)
        if source.H5File!=H5File:
            print("Json and H5 filenames are not matched... There is currently no check that the paths are ok here so make sure to pass the entire path if doing this.")

    # import some default parameters defined the ptemcee handler script
    from SYNEX.SYNEX_PTMC import run_params_default # , waveform_params_default

    # Create the default json file parameters
    json_default_dict = {}
    json_default_dict["run_params"] = run_params_default

    # Get additional parameters and add to json param dict
    # these dictionaries are 1. Base source params 2. waveform params 3. extras for fisher function
    # Therefore can rename the dictionaries for clarity, and then delete the key,value unwrap next
    # for the run_params update, since the run_params are handed to the function at call time.
    # But double check that the run_params are not not in the returned waveform params list etc.
    [param_dict, waveform_params, _ ] = ClassesToParams(source,detector,"Inference")
    json_default_dict["source_params"] = param_dict
    json_default_dict["waveform_params"] = waveform_params
    json_default_dict["prior_params"] = inference_params

    # Add some missing kwargs in the LISAnoise subdictionary in waveform params. Not really sure what this does but it is needed.
    # Need to fiure out if this is needed in theFisher stuff too - it will change where we put the defaults etc at detector initialization.
    json_default_dict["waveform_params"]["LISAnoise"]["lowf_add_pm_noise_f0"]=0.0
    json_default_dict["waveform_params"]["LISAnoise"]["lowf_add_pm_noise_alpha"]=2.0

    # Update json dict field defaults where needed
    for key,value in waveform_params.items():
        if key in json_default_dict["run_params"]:
            json_default_dict["run_params"][key] = value # these are not in the
        # elif key in json_default_dict["waveform_params"]:
        #     json_default_dict["waveform_params"][key] = value

    # Change now any keys set in the run time dictionary of kwargs (this could plot flags, run params values, no. of walkers etc)
    # This is NOT meant for binary params or prior params - these need to be specified at the highest script level
    for key,value in RunTimekwargs.items():
        if key in json_default_dict["run_params"]:
            json_default_dict["run_params"][key] = value
        elif key in json_default_dict["waveform_params"]:
            json_default_dict["waveform_params"][key] = value
        elif key in json_default_dict["source_params"]:
            json_default_dict["source_params"][key] = value

    # Set the output file and directory location in the json param list
    H5FilePath="/".join(source.H5File.split("/")[:-1])+"/"
    H5FileName=source.H5File.split("/")[-1]
    H5FileName=".".join(H5FileName.split(".")[:-1])
    json_default_dict["run_params"]["out_dir"] = H5FilePath
    json_default_dict["run_params"]["out_name"] =  H5FileName # this needs to not have the '.h5' added at the end to work

    # Write the json file only if master node or not mpi
    if IsMaster:
        print("Writting params to",source.JsonFile)
        with open(source.JsonFile, 'w') as f:
            json.dump(json_default_dict, f, indent=2)
        f.close()
    else:
        time.sleep(10) # MPI.COMM_WORLD.Barrier()???

def RunFoMOverRange(source,detector,ParamDict,FigureOfMerit='SNR',RunGrid=False,inference_params=None,**InferenceTechkwargs):
    """
    Function to run a Figure of Merit (FoM) over a range of parameter values.

    An example use would be sky localization error, using a full mcmc inference,
    over a grid of values for spins chi1 and chi2.

    Parameters
    ----------

    source : SYNEX 'source' class

    detector : SYNEX 'detector' class

    ParamDict : Dictionary
        Each field has key equal to the parameter name you want to loop over, this can be any field in the viable fields for the source and detector classes, including for example waveform approximant, total binary mass, or TDI type.
        The value attributed to each field is ether:
            - A single value.
                In this case the value is taken as a change to an existing parameter in the source or
                detector class. The inference will run by setting the json values from the source and
                detector classes, but then overwrite them using these set parameters. This is to save
                speed in the iterations over the loop by not having to create many versions of the
                same source and detector classes through runtime.
            - An array of three values.
                The first and second values are the min and max values, and the third value is the
                number of steps. This range specifies over what values of a particular variable the
                FoM should be calculated over. A return value - the value the
                parameter should return to at the end of the loop (if e.g. running several independent
                loops not in a grid) is taken as the midpoint between the min and max value.
            - An array of four values.
                Same as array with three values, but the fourth value if now the return value.
            - An array of values with first element = "custom".
                A custom range e.g. uneaven sampling. THIS IS TO BE IMPLEMENTED.

            NB: the return value in the grid case is not well defined - this will be developed in the future.

        FigureOfMerit : str
            A string indicated which FoM should be taken. Values are either:
             - "SNR"
                 Signal to noise ratio between the noise stored in the detector class, and the waveform
                 generated by the waveform parameters stored in both the detector and source classes.
             - "SkyArea"
                 A Fisher inf matrix estimate of the sky localization , e.g. the error ellipsoid when
                 infering the source location using the detector.
             - "SkyLocInfer"
                 A full MCMC inference of the multimodal posteriors using a set of defined parameter priors
                 that should also be specified (see "inference_params" definition). If singular values are specified
                 in "ParamDict", then the inference routine will take these values preferentially to source and
                 detector values in it's stored dictionaries. The sampler used here is "ptemcee", developed by
                 Sylvain Marsat in tandem with "lisabeta".

        RunGrid : Bool
            A boolean to indicate if a 2D grid should be taken over a set of parameters for e.g. a surface plot in
            the end. If false, the program will create a 1D array of values for the FoM over each param specified in
            ParamDict.

        inference_params : dict
            A dictionary containing the prior information for parameters to be inferred if using the "SkyLocInfer"
            FoM. The dictionary should have the fields:
             - "infer_params" - array of param names to infer
             - "params_range" - array of tuples of dim 2 for [min, max]
             - "prior_type" - array of prior types e.g. "uniform", "sin", "cos", etc.
             - "wrap_params" - array of flags to wrap the parameter or not

             NB: Each array should be organized in order of param names listed in "infer_params".
             NB: The dictionary of infered parameters is assumed to be unchanged through the ensemble
             of loops being asked for, to maintain apples-to-apples in analysis of results.

        InferenceTechkwargs : dict
            A dictionary of extra parameters used for inference. This includes parameters like number of walkers, flag to include thinning, etc.

        TO DO:
        ------
        - mpi.scatter for RunFoMOverRange implementation with mpi. this could handle automatically sending each job to different nodes instead of
        depending on mpi.map in the inference part. This is important since the parallel side of things seems to be more efficient when scattering
        jobs rather than parallelizing the inference...
    """

    if FigureOfMerit=="SkyLocInfer" and not inference_params:
        raise ValueError("You need a dictionary of inference parameter information to run inference for sky localization. Make sure to include fields for prior type, prior ranges, parameter name and a flag to wrap the parameter or not.")
    BaseParamDict = {}
    LoopVars = []
    ReturnValue = {}
    for key,values in ParamDict.items():
        if isinstance(values, (float,int,str)) and len(values)==1:
            BaseParamDict[key] = copy.deepcopy(values)
            ReturnValue[key] = copy.deepcopy(values)
        else:
            LoopVars.append(key)
            if len(values)==4:
                ReturnValue[key] = copy.deepcopy(values[3])
            else:
                ReturnValue[key] = (copy.deepcopy(values[0])+copy.deepcopy(values[1]))/2. # include statements that get the return values from the fields in source or detector if they exist

    # Check how many params for grid case are not more than 2 -- TO BE DEVELOPED
    if RunGrid and len(LoopVars)>2:
        raise ValueError("If you want to run a grid make sure there are only 2 variables being looped over. This will be developed in the future for larger grud dimensions, but for now please only specify 2 param ranges.")

    FoMOverRange = {}

    # Record the variables that will changed, and if it is a grid or linear range
    if isinstance(LoopVars, list):
        FoMOverRange["LoopVariables"] = LoopVars
    else:
        FoMOverRange["LoopVariables"] = LoopVars.tolist()
    FoMOverRange["IsGrid"] = RunGrid

    # Do the loop over each variable
    if RunGrid:
        Var1,Var2 = LoopVars[0],LoopVars[1]
        # Extract the first loop variable
        Var1_min = ParamDict[Var1][0]
        Var1_max = ParamDict[Var1][1]
        Var1_nSteps = ParamDict[Var1][2]
        if Var1 == "maxf" or Var1 == "M" or Var1 == "m1" or Var1 == "m2" or Var1 == "DeltatL_cut":
            if Var1 == "DeltatL_cut":
                Var1Array = np.logspace(np.log10(-Var1_min), np.log10(-Var1_max), Var1_nSteps)
            else:
                Var1Array = np.logspace(np.log10(Var1_min), np.log10(Var1_max), Var1_nSteps)
        else:
            Var1Array= np.linspace(Var1_min, Var1_max, Var1_nSteps)

        # Extract the second loop variable
        Var2_min = ParamDict[Var2][0]
        Var2_max = ParamDict[Var2][1]
        Var2_nSteps = ParamDict[Var2][2]
        if Var2 == "maxf" or Var2 == "M" or Var2 == "m1" or Var2 == "m2" or Var2 == "DeltatL_cut":
            if Var2 == "DeltatL_cut":
                Var2Array = np.logspace(np.log10(-Var2_min), np.log10(-Var2_max), Var2_nSteps)
            else:
                Var2Array = np.logspace(np.log10(Var2_min), np.log10(Var2_max), Var2_nSteps)
        else:
            Var2Array= np.linspace(Var2_min, Var2_max, Var2_nSteps)

        # Create space to output into dictionary
        FoMOverRange[Var1+"_xs"] = Var1Array.tolist()
        FoMOverRange[Var2+"_xs"] = Var2Array.tolist()
        FoMOverRange["grid_ys"] = np.empty((Var2_nSteps,Var1_nSteps))

        if FigureOfMerit == 'SkyLocInfer':
            BaseParamDict.update(InferenceTechkwargs)

        for Loop1ID in range(Var1_nSteps):
            for Loop2ID in range(Var2_nSteps):

                if FigureOfMerit == 'SkyLocInfer':
                    # Make sure json and output file locations exist
                    LisaBetaPath = os.path.dirname(os.path.realpath(__file__))
                    print(LisaBetaPath)
                    if not os.path.isdir(os.path.join(LisaBetaPath+"/../inference_param_files/GridInference_"+Var1+"_"+Var2+"/")):
                        os.mkdir(os.path.join(LisaBetaPath+"/../inference_param_files/GridInference_"+Var1+"_"+Var2+"/"))
                    if not os.path.isdir(os.path.join(OutFileLoc+"/../inference_data/GridInference_"+Var1+"_"+Var2+"/")):
                        os.mkdir(os.path.join(OutFileLoc+"/../inference_data/GridInference_"+Var1+"_"+Var2+"/"))

                # Change the variable that needs to be changed
                if Var1 == "DeltatL_cut":
                    BaseParamDict[Var1] = -Var1Array[Loop1ID]
                else:
                    BaseParamDict[Var1] = Var1Array[Loop1ID]
                if Var2 == "DeltatL_cut":
                    BaseParamDict[Var2] = -Var2Array[Loop2ID]
                else:
                    BaseParamDict[Var2] = Var2Array[Loop2ID]

                # Calculate the figure of merit, remove special case where a tested DeltatL_cut is larger than time in detection band
                [source_param_dict, _ , _ ] = ClassesToParams(source,detector,"Inference", **BaseParamDict)
                fLow, fHigh = lisa.FrequencyBoundsLISATDI_SMBH(source_param_dict)
                if Var1 == "DeltatL_cut" and pytools.funcNewtonianfoft(source_param_dict["m1"], source_param_dict["m2"], Var1Array[Loop1ID])<fLow:
                    print("Time cut is larger than detectable life span - setting FoM to 0")
                    FoM = 1e-30
                elif Var2 == "DeltatL_cut" and pytools.funcNewtonianfoft(source_param_dict["m1"], source_param_dict["m2"], Var2Array[Loop2ID])<fLow:
                    print("Time cut is larger than detectable life span - setting FoM to 0")
                    FoM = 1e-30
                elif FigureOfMerit == "SNR":
                    FoMOverRange["FoM"] = "SNR"
                    FoM = ComputeSNR(source,detector,**BaseParamDict)
                    # lisa.GenerateLISATDISignal_SOBH(params, **waveform_params) -- this is the function that has a field for snr: tdisignal['SNR']
                elif FigureOfMerit == "SkyArea":
                    FoMOverRange["FoM"] = "SkyArea"
                    fishercov = GetFisher_smbh(source, detector, **BaseParamDict)
                    FoM = lisatools.sky_area_cov(fishercov, sq_deg=True, n_sigma=None, prob=0.90)
                elif FigureOfMerit == "SkyLocInfer":
                    FoMOverRange["FoM"] = "SkyLocInfer"
                    # Create the output file name - to change on each iteration so it's clear what injection values are
                    [param_dict, _ , _ ] = ClassesToParams(source,detector,"Inference")
                    param_keys = param_dict.keys()
                    param_keys_mod = param_keys[param_keys!=Var1 and param_keys!=Var2]
                    OutFileName = "GridInference_"+Var1+"_"+Var2+"/InjValues_" + Var1 + "_" + str(Var1Array[Loop1ID]) + "_" + Var2 + "_" + str(Var2Array[Loop2ID])
                    OutFileName = OutFileName + ["_" + param + "_" + value for param,value in zip([param_keys_mod,BaseParamDict[param_keys_mod]])]
                    # Manually write the output file location and get the filepath
                    [InJsonFilePath, OutFileLoc, OutFileName] = WriteParamsToJson(source,detector,inference_params,OutFileName,**BaseParamDict)
                    # Makesure target folder exists
                    if not os.path.isdir(os.path.join(OutFileLoc+"GridInference_"+Var1+"_"+Var2+"/")):
                        os.mkdir(os.path.join(OutFileLoc+"GridInference_"+Var1+"_"+Var2+"/"))
                    # Iteratively running inference so no plots pls
                    DoPlots = False
                    # Run the inference step which will dump data into the output folder
                    RunInference(inference_params=inference_params, Plots=DoPlots, OutFileName=OutFileName,InJsonFilePath=InJsonFilePath)
                    # Now need to convert the data to some useful format
                    FoM = GetPosteriorVal(OutFileLoc+OutFileName) # This is a dictionary of mean, mode, median, quartiles, and standard dev (using mean)
                elif FigureOfMerit == "SkyModeLikelihoodsBETA":
                    FoMOverRange["FoM"] = "SkyModeLikelihoodsBETA" # Include beta so later we can include an option to check lamda, psi and other degeneracies.
                    lnL_skymodes,params_skymode = GetSkyMultiModeProbFromClasses(source, detector, **BaseParamDict)
                    InjMode = lnL_skymodes[(1,0)]
                    InjMode += lnL_skymodes[(1,1)]
                    InjMode += lnL_skymodes[(1,2)]
                    InjMode += lnL_skymodes[(1,3)]
                    BetaRefMode = lnL_skymodes[(-1,0)]
                    BetaRefMode += lnL_skymodes[(-1,1)]
                    BetaRefMode += lnL_skymodes[(-1,2)]
                    BetaRefMode += lnL_skymodes[(-1,3)]
                    FoM = lnL_skymodes[(-1,0)] - lnL_skymodes[(1,0)] # Primary reflection is most important - this is always negative. I think the logic with the the log10 means that taking the difference between all will always capture where reflections happen and not be washed out by others like Sylvain was saying... Need to discusss this with him more- maybe along with discussion of what you changed in PTEMCEE.
                    # FoM = BetaRefMode - InjMode
                elif FigureOfMerit == "SkyModeLikelihoodsLAMBDA":
                    FoMOverRange["FoM"] = "SkyModeLikelihoodsLAMBDA" # Include beta so later we can include an option to check lamda, psi and other degeneracies.
                    lnL_skymodes,params_skymode = GetSkyMultiModeProbFromClasses(source, detector, **BaseParamDict)
                    InjMode = lnL_skymodes[(1,0)]
                    InjMode += lnL_skymodes[(-1,0)]
                    lambdaRefMode1 = lnL_skymodes[(1,1)]
                    lambdaRefMode1 += lnL_skymodes[(-1,1)]
                    lambdaRefMode2 = lnL_skymodes[(1,2)]
                    lambdaRefMode2 += lnL_skymodes[(-1,2)]
                    lambdaRefMode3 = lnL_skymodes[(1,3)]
                    lambdaRefMode3 += lnL_skymodes[(-1,3)]
                    # FoM = lnL_skymodes[(-1,2)] - lnL_skymodes[(1,0)] # antipodal mode is most important
                    FoM = lambdaRefMode1 + lambdaRefMode2 + lambdaRefMode3 - InjMode
                    # FoM = np.max([lnL_skymodes[(1,1)], lnL_skymodes[(-1,1)], lnL_skymodes[(1,2)], lnL_skymodes[(-1,2)], lnL_skymodes[(1,3)], lnL_skymodes[(-1,3)]]) - lnL_skymodes[(1,0)]
                else:
                    raise ValueError("FoM requested: " + FigureOfMerit + " not implemented. Try a different FoM.")
                FoMOverRange["grid_ys"][Loop2ID][Loop1ID] = FoM

        # Change output to list
        FoMOverRange["grid_ys"] = FoMOverRange["grid_ys"].tolist()

        # Return to initial values
        BaseParamDict[Var1] = ReturnValue[Var1]
        BaseParamDict[Var2] = ReturnValue[Var2]

    else:
        if FigureOfMerit == 'SkyLocInfer':
            BaseParamDict.update(InferenceTechkwargs)

        for Var in LoopVars:
            if FigureOfMerit == 'SkyLocInfer':
                # Make sure json and output file locations exist
                LisaBetaPath = os.path.dirname(os.path.realpath(__file__))
                if not os.path.isdir(os.path.join(LisaBetaPath+"/../inference_param_files/RangeInference_"+Var+"/")):
                    os.mkdir(os.path.join(LisaBetaPath+"/../inference_param_files/RangeInference_"+Var+"/"))
                if not os.path.isdir(os.path.join(LisaBetaPath+"/../inference_data/RangeInference_"+Var+"/")):
                    os.mkdir(os.path.join(LisaBetaPath+"/../inference_data/RangeInference_"+Var+"/"))

            # Extract the loop variable
            Var_min = ParamDict[Var][0]
            Var_max = ParamDict[Var][1]
            Var_nSteps = ParamDict[Var][2]
            if Var == "maxf" or Var == "M" or Var == "m1" or Var == "m2": #  or Var == "DeltatL_cut":
                VarArray= np.logspace(np.log10(Var_min), np.log10(Var_max), Var_nSteps)
            else:
                VarArray= np.linspace(Var_min, Var_max, Var_nSteps)

            # Create space to output into dictionary
            FoMOverRange[Var+"_xs"] = VarArray.tolist()
            FoMOverRange[Var+"_ys"] = []

            for LoopID in range(Var_nSteps):
                # Change the variable that needs to be changed
                BaseParamDict[Var] = VarArray[LoopID]

                # Calculate the figure of merit, remove special case where a tested DeltatL_cut is larger than time in detection band
                [source_param_dict, _ , _ ] = ClassesToParams(source,detector,"Inference", **BaseParamDict)
                fLow, fHigh = lisa.FrequencyBoundsLISATDI_SMBH(source_param_dict)
                if Var == "DeltatL_cut" and pytools.funcNewtonianfoft(source_param_dict["m1"], source_param_dict["m2"], -VarArray[LoopID])<fLow:
                    print("Time cut is larger than detectable life span - setting FoM to 0")
                    FoM = 0.
                elif FigureOfMerit == "SNR":
                    FoMOverRange["FoM"] = "SNR"
                    FoM = ComputeSNR(source,detector,**BaseParamDict)
                elif FigureOfMerit == "SkyArea":
                    FoMOverRange["FoM"] = "SkyArea"
                    fishercov = GetFisher_smbh(source, detector, **BaseParamDict)
                    FoM = lisatools.sky_area_cov(fishercov, sq_deg=True, n_sigma=None, prob=0.90)
                elif FigureOfMerit == "SkyLocInfer":
                    FoMOverRange["FoM"] = "SkyLocInfer"
                    # Create the output file name - to change on each iteration so it's clear what injection values are
                    [source_param_dict, _ , _ ] = ClassesToParams(source,detector,"Inference")
                    for key in source_param_dict: # Update values here that were modified at run time for speed - note though BaseParamDict also has parameters like thinning etc so can't just use dict.update.
                        if key in BaseParamDict:
                            source_param_dict[key] = BaseParamDict[key]
                    _=source_param_dict.pop(Var, None) # source_param_dict contains all parameter values in this inference run, but Var is not updated, so need to pop
                    JsonFileName = "RangeInference_"+Var+"/InjValues_" + Var + "_" + str(round(10.*VarArray[LoopID])/10.)
                    for param,value in source_param_dict.items():
                        JsonFileName = JsonFileName + "_" + param + "_" + str(round(10.*value)/10.)

                    # Manually write the output file location and get the filepath
                    [JsonFileAndPath, OutFileLoc, OutFileName] = WriteParamsToJson(source,detector,inference_params,JsonFileName,**BaseParamDict)

                    # Iteratively running inference so no plots pls
                    DoPlots = False
                    # Run the inference step which will dump data into the output folder
                    RunInference(inference_params=inference_params, Plots=DoPlots, OutFileName=OutFileName,JsonFileAndPath=JsonFileAndPath)
                    # Now need to convert the data to some useful format
                    FoM = GetPosteriorVal(OutFileLoc+OutFileName) # This is a dictionary of mean, mode, median, quartiles, and standard dev (using mean)
                elif FigureOfMerit == "SkyModeLikelihoodsBETA":
                    FoMOverRange["FoM"] = "SkyModeLikelihoodsBETA" # Include beta so later we can include an option to check lamda, psi and other degeneracies.
                    lnL_skymodes,params_skymode = GetSkyMultiModeProbFromClasses(source, detector, **BaseParamDict)
                    InjMode = lnL_skymodes[(1,0)]
                    InjMode += lnL_skymodes[(1,1)]
                    InjMode += lnL_skymodes[(1,2)]
                    InjMode += lnL_skymodes[(1,3)]
                    BetaRefMode = lnL_skymodes[(-1,0)]
                    BetaRefMode += lnL_skymodes[(-1,1)]
                    BetaRefMode += lnL_skymodes[(-1,2)]
                    BetaRefMode += lnL_skymodes[(-1,3)]
                    FoM = BetaRefMode - InjMode
                elif FigureOfMerit == "SkyModeLikelihoodsLAMBDA":
                    FoMOverRange["FoM"] = "SkyModeLikelihoodsLAMBDA" # Include beta so later we can include an option to check lamda, psi and other degeneracies.
                    lnL_skymodes,params_skymode = GetSkyMultiModeProbFromClasses(source, detector, **BaseParamDict)
                    InjMode = lnL_skymodes[(1,0)]
                    InjMode += lnL_skymodes[(-1,0)]
                    lambdaRefMode1 = lnL_skymodes[(1,1)]
                    lambdaRefMode1 += lnL_skymodes[(-1,1)]
                    lambdaRefMode2 = lnL_skymodes[(1,2)]
                    lambdaRefMode2 += lnL_skymodes[(-1,2)]
                    lambdaRefMode3 = lnL_skymodes[(1,3)]
                    lambdaRefMode3 += lnL_skymodes[(-1,3)]
                    FoM = lambdaRefMode1 + lambdaRefMode2 + lambdaRefMode3 - InjMode
                    # FoM = lambdaRefMode1 - InjMode
                else:
                    raise ValueError("FoM requested: " + FigureOfMerit + " not implemented. Try a different FoM.")

                FoMOverRange[Var+"_ys"].append(FoM)

            # Return to initial value
            BaseParamDict[Var] = ReturnValue[Var]

    return FoMOverRange

def GetSkyMultiModeProbFromClasses(source, detector, **kwargs):
    # Get the parameters out of the classes and assign if given in kwargs dict
    [param_dict, waveform_params, extra_params] = ClassesToParams(source,detector,"Inference",**kwargs)

    # Define the likelihood class if not already existant
    likelihoodClass = lisa.LikelihoodLISASMBH(param_dict, **waveform_params)
    return lisatools.func_loglikelihood_skymodes(likelihoodClass)

def GetSkyMultiModeProbFromJson(FileName, **kwargs):

    # Check filenames
    JsonFileLocAndName,H5FileLocAndName = CompleteLisabetaDataAndJsonFileNames(FileName)

    # Read contents of file
    with open(JsonFileAndPath, 'r') as input_file:
        input_params = json.load(input_file)

    # Extract params and complete
    source_params = input_params['source_params']
    source_params = pytools.complete_mass_params(source_params)
    prior_params = input_params['prior_params']
    input_params['waveform_params'] = inference.waveform_params_json2py(
                                                input_params['waveform_params'])
    waveform_params = input_params['waveform_params'].copy()

    # Define the likelihood class if not already existant
    likelihoodClass = lisa.LikelihoodLISASMBH(source_params, **waveform_params)

    # Return
    return lisatools.func_loglikelihood_skymodes(likelihoodClass)

def GetFisher_smbh(source, detector, **kwargs):
    """ Function to grab the fisher cov matrix from lisabeta

    Parameters
    ---------
    Source : SYNEX source object.

    Detector : SYNEX detector object.

    freqs=None
    steps=default_steps
    list_params=default_list_fisher_params
    list_fixed_params=[]
    Lframe=False
    prior_invcov=None
    """
    # Get the parameters out of the classes and assign if given in kwargs dict
    [param_dict, waveform_params, extra_params] = ClassesToParams(source,detector,"Fisher",**kwargs)

    # Call the fisher function
    return lisa_fisher.fisher_covariance_smbh(param_dict, **waveform_params)

def GetSMBHGWDetection(source, detector, **kwargs):
    """ Function to grab the measured (TDI) waveform from a SYNEX binary object and SYNEX GW detector object

    Parameters
    ---------
    Source : SYNEX source object.

    Detector : SYNEX detector object.
    """
    # Get the parameters out of the classes and assign if given in kwargs dict
    [param_dict, waveform_params, extra_params] = ClassesToParams(source, detector, "Fisher", **kwargs) # this might need to be changed to 'Base' or something...

    # This is a workaround for flexibility with the input dicts. Not sure how fast it is- to be improved later!
    s = ' '
    for key,value in waveform_params.items():
        if isinstance(value, str):
            s = s + ', ' +  key + '="' + str(value) + '"'
        else:
            s = s + ', ' +  key + '=' + str(value)
    return eval('lisa.GenerateLISATDI_SMBH(param_dict' + s + ',extra_params=extra_params)')

def GetSMBHGWDetection_FreqSeries(source, detector, freqs=None, **kwargs):
    """ Function to grab the measured (TDI) waveform from a SYNEX binary object and SYNEX GW detector object

    Parameters
    ---------
    Source : SYNEX source object.

    Detector : SYNEX detector object.
    """

    # Get the parameters out of the classes and assign if given in kwargs dict
    [param_dict, waveform_params, extra_params] = ClassesToParams(source, detector, "Fisher", **kwargs) # this might need to be changed to 'Base' or something...

    # Generate frequencies if not already specified or if specified as a list of requirements, e.g. ['linear', fmin, fmax] (or something)
    if freqs is None:
        freqs = ['linear', None]
    if not isinstance(freqs, np.ndarray):
        freqs = GenerateFreqs(freqs, param_dict, **waveform_params)

    # This is a workaround for flexibility with the input dicts. Not sure how fast it is- to be improved later!
    s = ' '
    for key,value in waveform_params.items():
        if isinstance(value, str):
            s = s + ', ' +  key + '="' + str(value) + '"'
        else:
            s = s + ', ' +  key + '=' + str(value)
    print(s)
    # return lisa.GenerateLISATDIFreqseries_SMBH(param_dict, freqs,
    # timetomerger_max=1.0, minf=1e-05, t0=0.0, tref=0.0, phiref=0.0, fref_for_phiref=0.0, fref_for_tref=0.0,
    # force_phiref_fref=True, toffset=0.0, acc=0.0001, approximant="IMRPhenomHM", DeltatL_cut=14400.0, TDI="TDIAET",
    # order_fresnel_stencil=0, LISAconst="Proposal", responseapprox="full", frozenLISA=False, TDIrescaled=True,
    # LISAnoise={'InstrumentalNoise': 'SciRDv1', 'WDbackground': True, 'WDduration': 3.0, 'lowf_add_pm_noise_f0': 0.0,
    # 'lowf_add_pm_noise_alpha': 2.0)
    return eval('lisa.GenerateLISATDIFreqseries_SMBH(param_dict, freqs' + s + ')') #  + ', **waveform_params)') # ',extra_params=extra_params)')

def ComputeSNR(source, detector, freqs=None, Lframe=False, ReturnAllVariable=False, **kwargs):
    # Sort into dictionaries params, freqs
    [params, waveform_params, extra_params] = ClassesToParams(source,detector,"Fisher",**kwargs)

    if Lframe:
        raise ValueError('Lframe not implemented yet, sorry.')

    # Parameters
    m1 = params['m1']
    m2 = params['m2']
    M = m1 + m2
    q = m1 / m2
    params['M'] = M
    params['q'] = q

    LISAnoise = waveform_params.pop('LISAnoise', pyLISAnoise.LISAnoiseSciRDv1)
    TDI = waveform_params.get('TDI', 'TDIAET')
    TDIrescaled = waveform_params.get('TDIrescaled', True)
    LISAconst = waveform_params.get('LISAconst', pyresponse.LISAconstProposal)

    # Default frequencies (safe and slow):
    # linear spacing, adjust deltaf=1/(2T) with T approximate duration
    if freqs is None:
        freqs = ['linear', None]
    if not isinstance(freqs, np.ndarray):
        freqs = GenerateFreqs(freqs, params, **waveform_params)

    # Generate tdi freqseries
    tdifreqseries_base = lisa.GenerateLISATDIFreqseries_SMBH(params, freqs, **waveform_params)

    # Compute noises
    noise_evaluator = pyLISAnoise.initialize_noise(LISAnoise,
                                        TDI=TDI, TDIrescaled=TDIrescaled,
                                        LISAconst=LISAconst)
    Sn1_vals, Sn2_vals, Sn3_vals = pyLISAnoise.evaluate_noise(
                          LISAnoise, noise_evaluator, freqs,
                          TDI=TDI, TDIrescaled=TDIrescaled, LISAconst=LISAconst)

    # Modes and channels
    modes = tdifreqseries_base['modes']
    channels = ['chan1', 'chan2', 'chan3']
    Snvals = {}
    Snvals['chan1'] = Sn1_vals
    Snvals['chan2'] = Sn2_vals
    Snvals['chan3'] = Sn3_vals
    h_full = {}
    for chan in channels:
        h_full[chan] = np.zeros_like(freqs, dtype=complex)
        for lm in modes:
            h_full[chan] += tdifreqseries_base[lm][chan]

    # Compute SNR
    # We do not assume that deltaf is constant
    df = np.diff(freqs)
    SNR = 0.
    for chan in channels:
        SNR += 4*np.sum(df * np.real(h_full[chan] * np.conj(h_full[chan]) / Snvals[chan])[:-1])

    if ReturnAllVariable:
        return np.sqrt(SNR), freqs, h_full, Snvals
    else:
        # wftdi = GetSMBHGWDetection(source, detector, **kwargs)
        # return lisa.computeSNR(wftdi) # ,LISAconstellation=pyresponse.LISAconstProposal,LISAnoise=pyLISAnoise.LISAnoiseProposal,LDCnoise=None,TDIrescaled=False,Nf=None)
        return np.sqrt(SNR)

def GenerateFreqs(freqs, params, **waveform_params):
    """
    Function to automatiucally generate the frequencies needed for caluclations like SNR.

    Params
    ------
        freqs : list
        Object to specify what kind of frequencies is required. For example:
            freqs = ['linear', None]
        requests linear spacing, with no pre-required bounds. The generator will find the total
        time of the signal from the params and waveform_params dictionaries, and then generate a
        frequeny range that covers the full waveform.

        params : dict
        Dictionary of parameters describving the source (mass, spin etc).

        waveform_params : dict
        Dictionary of parameters needed for waveform generation (approximant, modes, etc.)
    """

    # Determine (2,2) frequency bounds based on frequency and time limits
    fLow, fHigh = lisa.FrequencyBoundsLISATDI_SMBH(params, **waveform_params)

    # Set of harmonics to be returned, default for approximants
    # TODO: add error raising if incompatibility with mode content of approx.
    modes = waveform_params.get('modes', None)
    if modes is None:
        if (not 'approximant' in waveform_params) or waveform_params['approximant']=='IMRPhenomD':
            modes = [(2,2)]
        elif waveform_params['approximant']=='IMRPhenomHM':
            modes = [(2,2), (2,1), (3,3), (3,2), (4,4), (4,3)]
    mmax = max([lm[1] for lm in modes])

    # fHigh is for (2,2) -- we use one single frequency array for all modes
    fHigh = mmax/2. * fHigh

    # Format: freqs = ['log', npt], npt=None for auto setting with T
    if freqs[0]=='log':
        if freqs[1] is not None:
            nptlog = int(freqs[1])
        else:
            # Estimate length and choose deltaf according to Nyquist
            # NOTE: No safety margin, so set timetomerger_max must carefully
            T = pytools.funcNewtoniantoff(params["m1"], params["m2"], fLow)
            deltaf = 1./(2*T)
            nptlog = np.ceil(np.log(fHigh/fLow) / np.log(1. + deltaf/fLow))
        freqs = pytools.logspace(fLow, fHigh, nptlog)

    # Format: freqs = ['linear', deltaf], deltaf=None for auto setting with T
    if freqs[0]=='linear':
        if freqs[1] is not None:
            deltaf = freqs[1]#int(freqs[1])
        else:
            # Estimate length and choose deltaf according to Nyquist
            # NOTE: No safety margin, so set timetomerger_max must carefully
            T = pytools.funcNewtoniantoff(params["m1"], params["m2"], fLow)
            deltaf = 1./(2*T)
        freqs = np.arange(fLow, fHigh, deltaf)

    return freqs

def ComputeDetectorNoise(source, detector, freqs=None, Lframe=False, ReturnAllVariable=False, **kwargs):
    # Grab the variables from classes
    [params, waveform_params, extra_params] = ClassesToParams(source,detector,"Fisher",**kwargs)

    LISAnoise = waveform_params.get('LISAnoise', pyLISAnoise.LISAnoiseSciRDv1)
    TDI = waveform_params.get('TDI', 'TDIAET')
    TDIrescaled = waveform_params.get('TDIrescaled', True)
    LISAconst = waveform_params.get('LISAconst', pyresponse.LISAconstProposal)

    # Convert input parameters to either Lframe or SSBframe
    if not params.get('Lframe', False) and Lframe:
        params_base = lisatools.convert_SSBframe_to_Lframe(
                            params,
                            t0=waveform_params['t0'],
                            frozenLISA=waveform_params['frozenLISA'])
    elif params.get('Lframe', False) and not Lframe:
        params_base = lisatools.convert_Lframe_to_SSBframe(
                            params,
                            t0=waveform_params['t0'],
                            frozenLISA=waveform_params['frozenLISA'])
    else:
        params_base = params.copy()
    params_base = pytools.complete_mass_params(params_base)
    params_base = pytools.complete_spin_params(params_base)

    # Default frequencies (safe and slow):
    # linear spacing, adjust deltaf=1/(2T) with T approximate duration
    if freqs is None:
        freqs = ['linear', None]
    # Determine the frequencies to use for the overlaps, if not given as input
    if not isinstance(freqs, np.ndarray):
        freqs = GenerateFreqs(freqs, params_base, **waveform_params)

    # Compute noises
    noise_evaluator = pyLISAnoise.initialize_noise(LISAnoise,
                                        TDI=TDI, TDIrescaled=TDIrescaled,
                                        LISAconst=LISAconst)
    Sn1_vals, Sn2_vals, Sn3_vals = pyLISAnoise.evaluate_noise(
                          LISAnoise, noise_evaluator, freqs,
                          TDI=TDI, TDIrescaled=TDIrescaled, LISAconst=LISAconst)

    # Compute derivatives dh
    Snvals = {}
    Snvals['chan1'] = Sn1_vals
    Snvals['chan2'] = Sn2_vals
    Snvals['chan3'] = Sn3_vals

    detector.Snvals = Snvals

def fmaxFromTimeToMerger(source, detector, T_obs_end_to_merger=4.*60.*60., ReturnVal=False, **kwargs):
    # Get the relevant waveform data from classes, needed to ensure any other changes variables are propagated
    [param_dict, waveform_params, _ ] = ClassesToParams(source,detector,"Inference",**kwargs)

    # Calculate the frequency at which the time to merger is equal to what we want
    F = pytools.funcNewtonianfoft(param_dict["m1"], param_dict["m2"], T_obs_end_to_merger)

    if ReturnVal:
        # return the calculated value too if requested (used in RunFoMOverRange function)
        return F
    else:
        # Else modify the stored fmax in source class directly
        source.maxf = F

def GetSourceFromLisabetaData(FileName, **kwargs):

    # Check filenames
    JsonFileLocAndName,H5FileLocAndName=CompleteLisabetaDataAndJsonFileNames(FileName)

    # Extract data
    with open(JsonFileLocAndName, 'r') as f:
        input_params = json.load(f)
    f.close()

    # check if it was run on a different system, e.g. the cluster, so paths are wrong
    check_path = input_params["run_params"]["out_dir"].split("/SYNEX/")[-1]
    check_path_true = H5FileLocAndName.split("/SYNEX/")[-1]
    if check_path!=check_path_true:
        input_params["run_params"]["out_dir"] = "/".join(H5FileLocAndName.split("/")[:-1])+"/" # Needs trailing slash...
    # For now we do not re-write the file because this might break stuff...
    # print("Re-writting params to",JsonFileLocAndName)
    # with open(source.JsonFile, 'w') as f:
    #     json.dump(input_params, f, indent=2)
    # f.close()

    return ParamsToClasses(input_params,CollectionMethod="Inference",**kwargs)






#                 ##########################################
#                 #                                        #
#                 #   SYNEX -- GWEMOPT Utility functions   #
#                 #                                        #
#                 ##########################################


def TileSkyArea(source_or_kwargs,detectors=None,base_telescope_params=None,cloning_params=None):
    """
    Global function to tile skyarea. This will handle all extra tiling we might add,
    including handling series inputs over lists of sources and/or detectors.

    TO DO:
    ------
    1. Case where we inject a list of sources...
    2. Sort 'out_dirs' list if we have an injected long list of detectors where
       more than one variable is changing. Maybe a variable for cloning_params?
    """
    # Using MPI or not
    if MPI is not None:
        MPI_size = MPI.COMM_WORLD.Get_size()
        MPI_rank = MPI.COMM_WORLD.Get_rank()
        comm = MPI.COMM_WORLD
    else:
        MPI_size=1
        MPI_rank=0

    # Populate detectors if we have cloning parameters
    if cloning_params!=None:
        # Create base_telescope_params from default vals or whatever was handed to us
        # We preferentially choose base_telescope_params that was handed to us
        if base_telescope_params==None:
            if detectors==None:
                base_telescope_params=SYDs.Athena().__dict__
            elif isinstance(detectors,list):
                if len(detectors)>1:
                    print("Cloning first detector in list only...")
                base_telescope_params=detectors[0].__dict__
            else:
                base_telescope_params=detectors.__dict__

        # Base detector savefile and name
        BaseExFileName=base_telescope_params["ExistentialFileName"] if "ExistentialFileName" in base_telescope_params else "Athena_Base"
        if 'detector_config_struct' in base_telescope_params:
            BaseTelescopeName=base_telescope_params["detector_config_struct"]["telescope"]
        elif "telescope" in base_telescope_params:
            BaseTelescopeName=base_telescope_params["telescope"]
        else:
            BaseTelescopeName="Athena"

        # Start populating detector list -- list of lists, each embedded list corresponds to each requested cloning parameter
        out_dirs=[]
        detectors=[]
        for key,values in cloning_params.items():
            # In case single value
            if not isinstance(values,(list,np.ndarray)): values=[values]
            # Work out how many per cpu if we are on cluster
            Nvals=len(values)
            NValsPerCore=int(Nvals//MPI_size)
            CoreLenVals=[NValsPerCore+1 if ii<len(values)%MPI_size else NValsPerCore for ii in range(MPI_size)]
            CPU_ENDs=list(np.cumsum(CoreLenVals))
            CPU_STARTs=[0]+CPU_ENDs[:-1]
            # Make list of detector properties depending on MPI_rank
            dict_list = [base_telescope_params.copy() for _ in range(CPU_STARTs[MPI_rank],CPU_ENDs[MPI_rank])]
            for ii in range(CPU_STARTs[MPI_rank],CPU_ENDs[MPI_rank]): dict_list[ii-CPU_STARTs[MPI_rank]].update({"ExistentialFileName":BaseExFileName,
                          "NewExistentialFileName":".".join(BaseExFileName.split(".")[:-1])+"_"+key+"_"+str(ii+1)+"."+BaseExFileName.split(".")[-1],
                          key:values[ii],
                          "cloning key":key,
                          "cloning value":values[ii],
                          "telescope":BaseTelescopeName+"_"+key+"_"+str(ii+1)})
            if MPI_size>1: print("CPU",MPI_rank+1,"/",MPI_size,"has",len(dict_list),"detector objects for",key,flush=True)
            detectors+=[SYDs.Athena(**dict_ii) for dict_ii in dict_list]
            out_dirs+=[key]*len(dict_list)
        if MPI_size>1: print("CPU",MPI_rank+1,"/",MPI_size,"has",len(dict_list),"total detector objects.",flush=True)
    else:
        # No cloning params so see if we have input detectors
        if detectors==None and base_telescope_params==None:
            detectors=[SYDs.Athena()] if MPI_rank==0 else None
        elif detectors==None and base_telescope_params!=None:
            detectors=[SYDs.Athena(**base_telescope_params)] if MPI_rank==0 else None
        elif detectors!=None and not isinstance(detectors,list):
            detectors=[detectors] if MPI_rank==0 else None # if only one detector given
        elif detectors!=None and isinstance(detectors,list):
             # List of detectors given: split over all cpu.
            Nvals=len(detectors)
            NValsPerCore=int(Nvals//MPI_size)
            CoreLenVals=[NValsPerCore+1 if ii<len(values)%MPI_size else NValsPerCore for ii in range(MPI_size)]
            CPU_ENDs=list(np.cumsum(CoreLenVals))
            CPU_STARTs=[0]+CPU_ENDs[:-1]
            detectors=[detectors[ii] for ii in range(CPU_STARTs[MPI_rank],CPU_ENDs[MPI_rank])]

        # Default of no output directory extension
        out_dirs=[None]*len(detectors) if isinstance(detectors,list) else None

    # Init source if we need to
    if MPI_rank==0:
        # See if we have a source class or kwargs -- only master node since this can get memory heavy
        source=SYSs.SMBH_Merger(**source_or_kwargs) if isinstance(source_or_kwargs,dict) else source_or_kwargs
    else:
        source=None # So we only have master node running later if one source and/or detector given

    # Send source to workers
    if MPI_size>1: source = comm.bcast(source, root=0)

    # Calculate source flux data -- NB CTR is telescope dependent so included inside loop for coverage info later if ARF file changes
    if source!=None and not hasattr(source,"EM_Flux_Data"): source.GenerateEMFlux(fstart22=1e-4,TYPE="const",**{})

    # Make sure only cases where detectors and source are defined are run (i.e. if we handed cluster one object with many cluster).
    # We will optimise this later to ask gwemopt to run in parallel once we understand if this will speed up calculations for a single
    # telescope or if they employ parallel code only when you ask it to consider several (unique!) telescopes in the coverage scheme.
    if detectors!=None and source!=None:
        # Finish source EM calculations based on broadcast detectors only if ARF file does not change in detectors
        if len(set([detector.ARF_file_loc_name for detector in detectors]))==1: source.GenerateCTR(detectors[0].ARF_file_loc_name,gamma=1.7)

        # Loop over list of detectors
        for i in range(len(detectors)): go_params, map_struct, tile_structs, coverage_struct, detectors[i] = TileWithGwemopt(source,detectors[i],out_dirs[i])

    return detectors

def TileWithGwemopt(source,detector,outDirExtension=None):
    """
    Function to handle tiling through gwemopt.
    """

    # Get the right dicts to use
    t_tile_0=time.time()
    go_params,map_struct=PrepareGwemoptDicts(source,detector,outDirExtension)
    # go_params["doParallel"]=True
    # go_params["Ncores"]=2

    # Get segments -- need to understand what this is. Line 469 of binary file 'gwemopt_run'
    import SYNEX.segments_athena as segs_a
    go_params = segs_a.get_telescope_segments(go_params)

    # Get tile_structs
    if go_params["tilesType"]=="MaxProb":
        # Get tiling structs
        moc_structs = gwemopt.moc.create_moc(go_params, map_struct=map_struct)
        tile_structs = gwemopt.tiles.moc(go_params,map_struct,moc_structs)
    elif go_params["tilesType"]=="moc":
        moc_structs = gwemopt.moc.create_moc(go_params, map_struct=map_struct)
        tile_structs = gwemopt.tiles.moc(go_params,map_struct,moc_structs,doSegments=False) # doSegments=False ?? Only for 'moc'... Otherwise it calls gwemopt.segments.get_segments_tiles
        for telescope in tile_structs.keys():
            tile_structs[telescope] = segs_a.get_segments_tiles(go_params, go_params["config"][telescope], tile_structs[telescope])
    elif go_params["tilesType"]=="greedy":
        tile_structs = gwemopt.tiles.greedy(go_params,map_struct)
        go_params["Ntiles"] = []
        for telescope in go_params["telescopes"]:
            tile_structs[telescope] = segs_a.get_segments_tiles(go_params, go_params["config"][telescope], tile_structs[telescope]) # replace segs with our own
            go_params["config"][telescope]["tesselation"] = np.empty((0,3))
            tiles_struct = tile_structs[telescope]
            for index in tiles_struct.keys():
                ra, dec = tiles_struct[index]["ra"], tiles_struct[index]["dec"]
                go_params["config"][telescope]["tesselation"] = np.append(go_params["config"][telescope]["tesselation"],[[index,ra,dec]],axis=0)
            go_params["Ntiles"].append(len(tiles_struct.keys()))
    elif go_params["tilesType"]=="hierarchical":
        tile_structs = gwemopt.tiles.hierarchical(go_params,map_struct) # ,doSegments=False)
        go_params["Ntiles"] = []
        for telescope in go_params["telescopes"]:
            tile_structs[telescope] = segs_a.get_segments_tiles(go_params, go_params["config"][telescope], tile_structs[telescope]) # replace segs with our own
            go_params["config"][telescope]["tesselation"] = np.empty((0,3))
            tiles_struct = tile_structs[telescope]
            for index in tiles_struct.keys():
                ra, dec = tiles_struct[index]["ra"], tiles_struct[index]["dec"]
                go_params["config"][telescope]["tesselation"] = np.append(go_params["config"][telescope]["tesselation"],[[index,ra,dec]],axis=0)
            go_params["Ntiles"].append(len(tiles_struct.keys()))
    elif go_params["tilesType"]=="ranked":
        # Get tiling structs
        moc_structs = gwemopt.rankedTilesGenerator.create_ranked(go_params,map_struct)
        tile_structs = gwemopt.tiles.moc(go_params,map_struct,moc_structs,doSegments=False)
        for telescope in tile_structs.keys():
            tile_structs[telescope] = segs_a.get_segments_tiles(go_params, go_params["config"][telescope], tile_structs[telescope])
    elif go_params["tilesType"]=="galaxy":
        # Really not sure how this works and where segments are calculated... Use this method with care.
        map_struct, catalog_struct = gwemopt.catalog.get_catalog(go_params, map_struct)
        tile_structs = gwemopt.tiles.galaxy(go_params,map_struct,catalog_struct)
        for telescope in go_params["telescopes"]:
            # tile_structs[telescope] = segs_a.get_segments_tiles(go_params, go_params["config"][telescope], tile_structs[telescope])
            go_params["config"][telescope]["tesselation"] = np.empty((0,3))
            tiles_struct = tile_structs[telescope]
            for index in tiles_struct.keys():
                ra, dec = tiles_struct[index]["ra"], tiles_struct[index]["dec"]
                go_params["config"][telescope]["tesselation"] = np.append(go_params["config"][telescope]["tesselation"],[[index,ra,dec]],axis=0)

    # Allocate times to tiles
    tile_structs, coverage_struct = gwemopt.coverage.timeallocation(go_params, map_struct, tile_structs)

    # Get info for run
    detector = GetCoverageInfo(go_params, map_struct, tile_structs, coverage_struct, detector, source=source, verbose=False)

    # Add info about source to coverage summary
    SourceInfo={
                "source save file":source.ExistentialFileName,
                "source gpstime":source.gpstime,
                "source sky_map":source.sky_map,
                "source JsonFile":source.JsonFile,
                "source H5File":source.H5File,
                "source Is3D":source.do3D,
                "source true loc":[source.true_ra,source.true_dec]
                }
    if source.do3D:
        SourceInfo["source true loc"]+=[source.true_distance]
    detector.detector_source_coverage.update(SourceInfo)

    # Save detector
    detector.ExistentialCrisis()

    return go_params, map_struct, tile_structs, coverage_struct, detector

def PrepareGwemoptDicts(source,detector,outDirExtension=None):
    """
    Function to prepare GWELOPT dictionary from a source object and a detector objct.

    Need to adjust for case of updating skymaps.

    Include variable Tobs- gwemopt uses the option of several available windows... Maybe for scheduling time wtih telescope(s) or map updates.
    """
    # Prepare go_params from detector first
    go_params = copy.deepcopy(detector.detector_go_params)
    telescopes = [detector.detector_config_struct["telescope"]] # They want this to be a list
    go_params["telescopes"]=telescopes
    go_params["exposuretimes"]=np.array([detector.detector_config_struct["exposuretime"]]) # This line is necessary and not really sure why...
    go_params["config"]={detector.detector_config_struct["telescope"] : detector.detector_config_struct}

    # Ideal pixelation for detector FoV has pixel area <= FOV/4 (idk this for sure, I just made it up...)
    nside_arr=np.array([2**i for i in range(11,3,-1)])
    area_arr=hp.nside2pixarea(nside_arr,degrees=True)
    nside=nside_arr[np.searchsorted(area_arr, detector.detector_config_struct["FOV"]/4., side='left')-1]

    # Update missing params in go_params contained in source
    go_params["nside"]=nside # hp.pixelfunc.get_nside(source.map_struct["prob"])
    if source.gpstime!=None:
        go_params["gpstime"]=source.gpstime
    else:
        go_params["gpstime"]=detector.detector_config_struct["gps_science_start"]-source.DeltatL_cut/86400.
        source.gpstime=go_params["gpstime"]
    go_params["do3D"]=source.do3D
    go_params["true_ra"]=source.true_ra
    go_params["true_dec"]=source.true_dec
    go_params["true_distance"]=source.true_distance
    go_params["Tobs"]=np.array([0.,-source.DeltatL_cut/86400.]) # np.array([0.,100.]) # in pairs of [Tstart,Tend] for times in DAYS. We CAN pass more than one pair here.

    # Update the gwemopt output dir to comply with lisabeta h5 and json filenames
    # This will overwrite existing files with the same event name... Maybe later include checks?
    if not go_params["tilesType"] in ["hierarchical","greedy"]:
        ### DONT add subfolder structure to 'outputDir' if sampler based tilesType is requested:
        ### The pathlength becomes too long for pymultinest so filenames get truncated, leading to 'file not found' errors...
        ### "greedy" acutally uses emcee -- does this have the same problem?
        go_params["outputDir"]+="/"+".".join(source.JsonFile.split("/")[-1].split(".")[:-1])
        if not outDirExtension==None:
            go_params["outputDir"]+="_"+outDirExtension
    pathlib.Path(go_params["outputDir"]).mkdir(parents=True, exist_ok=True)
    if os.path.isfile(go_params["outputDir"]):
        print("WARNING: outputDir",go_params["outputDir"].split("/SYNEX/")[-1],"already exists... Files within will be overwriten.")

    # Now get map struct
    map_struct=source.map_struct

    # Run through the GWEMOPT checker for uniformity between go_params and sky_map (rescale nside etc)
    # We can also apply rotations etc here
    map_struct=gou.read_skymap(go_params,source.do3D,map_struct=map_struct)

    return go_params,map_struct

def GetCoverageInfo(go_params, map_struct, tile_structs, coverage_struct, detector, source=None, verbose=True):
    """
    Extract coverage information, including a count for tiles that cover source and
    tiles that cover source with exposure time at least the minimum for a threshold
    photon count.
    """
    # gwemopt.scheduler.summary(go_params, map_struct, coverage_struct)
    # gwemopt.plotting.coverage(go_params, map_struct, coverage_struct)

    # Extract some data structures
    source_pix = hp.ang2pix(go_params["nside"],np.rad2deg(source.lamda),np.rad2deg(source.beta),lonlat=True)
    cov_data=coverage_struct["data"] # Data is N*9 with columns ra,dec,mjd (exposure start),mag,exposureTime,field,tile_probs,airmass,program_id
    cov_ipix=coverage_struct["ipix"] # list because each sub-list is of variable length
    map_probs=map_struct["prob"]

    # Tiling info
    telescope=detector.detector_config_struct["telescope"]
    tile_struct=tile_structs[telescope]
    tile_keys_list=list(tile_struct.keys())
    pix_record=np.unique(np.concatenate([tile_struct[tile_id]["ipix"] for tile_id in tile_keys_list if len(tile_struct[tile_id]["ipix"])>0],axis=0))
    probs1=np.array([tile_struct[tile_id]["prob"] for tile_id in tile_keys_list])
    probs2=np.array([np.sum(map_probs[tile_struct[tile_id]["ipix"]]) for tile_id in tile_keys_list])
    TotProb1=np.sum(probs1) # Tile struct probabilities
    TotProb2=np.sum(probs2) # Mp struct probs summed using ipix in tile struct
    TotProb3=np.sum(map_probs[pix_record]) # Unique probabilities
    detector.detector_tile_struct=tile_struct
    if verbose:
        print("\n\n")
        print("#"*20+"  "+telescope+"  "+"#"*20+"\n")
        print("Total prob (tile_struct, map_struct, map_struct drop dups):",TotProb1,TotProb2,TotProb3)
        print("Params/Config checks:", telescope, go_params["config"][telescope]["tot_obs_time"], go_params["config"][telescope]["exposuretime"], go_params["config"][telescope]["tot_obs_time"]/go_params["config"][telescope]["exposuretime"], "Tot tiles avail:", len(tile_struct.keys()))
        print("\n")

    # Coverage info
    cov_ipix_array_unique=np.unique([pix for cov_ipix_el in cov_ipix for pix in cov_ipix_el])
    prob1=np.sum(cov_data[:,6])
    prob2=np.sum([map_probs[pix] for tile_pix in cov_ipix for pix in tile_pix]) # Probs according to coverage struct ipix
    prob3=np.sum([map_probs[pix] for pix in cov_ipix_array_unique]) # Same as prob2 but only unique pixels
    ex_Times = cov_data[:,4]
    ex_Times_unique = np.unique(ex_Times)
    detector.detector_coverage_struct=coverage_struct
    detector.detector_coverage_struct.update({"tile exposure times":ex_Times,"unique tile exposure times":ex_Times_unique,
                                              "normalised tile exposure times":[T_e/detector.detector_config_struct["exposuretime"] for T_e in ex_Times]})

    # Print summary if asked for
    if verbose:
        print(len(ex_Times),"/",len(cov_data[:,4]), "coverage tiles")
        print("Detector coverage prob (cov_struct data, map_struct, map_struct dup drop):",prob1,prob2,prob3)
        print(len(ex_Times_unique),"/",len(cov_data[:,4]), "unique coverage tiles")
        print(len(ex_Times_unique),"unique exposure times with (min, mean, max):",np.min(ex_Times_unique),np.mean(ex_Times_unique),np.max(ex_Times_unique))
        print("\n")

    ## Source tile info -- was it covered?
    UniqueSourceCoverageCount=0
    cov_source_tile=[]
    cov_source_pix=[]
    TotExpTime=0.
    for tile_ii,tile_pix in enumerate(cov_ipix):
        TotExpTime+=cov_data[tile_ii,4]
        # print(tile_ii, source_pix, tile_pix, cov_data[tile_ii,4], detector.detector_config_struct["exposuretime"], cov_data[tile_ii,4]/detector.detector_config_struct["exposuretime"])
        if source_pix in tile_pix:
            cov_source_tile += [tile_ii] # [tile_ii for tile_ii,tile_pix in enumerate(cov_ipix) if source_pix in tile_pix]
            cov_source_pix += list(tile_pix)
            if all([pix not in cov_source_pix for pix in tile_pix]):
                UniqueSourceCoverageCount+=1
    print("Total exposure time (d):",TotExpTime/(24*60*60))
    SourceTile_prob1=[cov_data[i,6] for i in cov_source_tile]
    SourceTile_prob2=[map_probs[pix] for pix in cov_source_pix]
    SourceTile_accum_prob1=np.sum(SourceTile_prob1)
    SourceTile_accum_prob2=np.sum(SourceTile_prob2)
    detector.detector_source_coverage={"telescope":telescope,"Source Tile Prob (by cov tiles)":SourceTile_prob1,
                        "Source Tile Prob (by cov pixels)":SourceTile_prob2,"tile ID containing source pixel":cov_source_tile,
                        "tile pixels containing source pixel":cov_source_pix, "No. source coverages":len(SourceTile_prob1),
                        "No. unique source coverages":UniqueSourceCoverageCount}

    # Print summary if asked for
    if verbose:
        print("--- Source Coverage ---")
        print("\n")
        print("Source covered by",telescope, len(cov_source_tile),"times.\n")
        print(" "*10+"cov_struct p data:")
        print(SourceTile_prob1,"\n")
        print(" "*10+"map_struct pixel p's:")
        print(SourceTile_prob2,"\n")
        print("Source tile accumulated p:",SourceTile_accum_prob1,SourceTile_accum_prob2)

    print("Coverage stuct exposures within range:",Time(cov_data[0,2],format='mjd', scale='utc').isot,Time(cov_data[-1,2],format='mjd', scale='utc').isot)
    print("start times checks:",Time(source.gpstime, format='gps', scale='utc').isot,Time(go_params["gpstime"], format='gps', scale='utc').isot)

    # Now add cuts to coverage info depending on photon flux
    if source!=None and len(cov_source_tile)>0:
        # Make sure source has relevant data and create it if not
        # if not hasattr(source,"EM_Flux_Data"): source.GenerateEMFlux(fstart22=1e-4,TYPE="const",**{})
        # if not hasattr(source,"CTR_Data"): source.GenerateCTR(detector.ARF_file_loc_name,gamma=1.7) # Should be include gamma as a source param? Would we ever want to change this at run time?
        source.GenerateEMFlux(fstart22=1e-4,TYPE="const",**{})
        source.GenerateCTR(detector.ARF_file_loc_name,gamma=1.7)

        # Calculate the exposuretimes for each tile that covers source
        SourceTileExpTimes=[cov_data[tile_ii,4] for tile_ii in cov_source_tile]
        SourceTileStartTimes=[86400.*(cov_data[tile_ii,2]-Time(go_params["gpstime"], format='gps', scale='utc').mjd) for tile_ii in cov_source_tile]

        # Get CTR data out from source and cut to Tobs -- times here are seconds to merger (<0)
        CTRs=source.CTR_Data["CTR"]
        CTR_times=list(86400.*(-np.array(source.EM_Flux_Data["xray_time"])/86400. - Time(source.gpstime, format='gps', scale='utc').mjd + Time(go_params["gpstime"], format='gps', scale='utc').mjd)) # Just in case these are different -- i.e. if we have a run with multiple sources within Tobs
        CTRs=[CTR for CTR,t in zip(CTRs,CTR_times) if t>=(go_params["Tobs"][0]*86400.) and t<=(go_params["Tobs"][-1]*86400.)]
        CTR_times=[t for t in CTR_times if t>=(go_params["Tobs"][0]*86400.) and t<=(go_params["Tobs"][-1]*86400.)] # because times are seconds TO MERGER (-ve)

        # Integrate CTR for each tile covering source to check no. of source photons received
        TileListTimes=[[ts+i*dur/49 for i in range(50)] for ts,dur in zip(SourceTileStartTimes,SourceTileExpTimes)]
        SourceTilePhotonCounts=[np.trapz(np.interp(TileTs,CTR_times,CTRs),TileTs) for TileTs in TileListTimes]

        detector.detector_source_coverage.update({"Source tile exposuretimes (s)":SourceTileExpTimes,
                                                  "Source photon counts":SourceTilePhotonCounts,
                                                  "Source tile start times (s)":SourceTileStartTimes,
                                                  "Start time (mjd)":Time(go_params["gpstime"], format='gps', scale='utc').mjd})
    else:
        detector.detector_source_coverage.update({"Source tile exposuretimes (s)":[],
                                                  "Source photon counts":[0],
                                                  "Source tile start times (s)":[],
                                                  "Start time (mjd)":Time(go_params["gpstime"], format='gps', scale='utc').mjd})

    # Return values
    return detector

def WriteSkymapToFile(map_struct,SkyMapFileName,go_params=None,PermissionToWrite=True):
    """
    So far have not included case where this is 3D. This will be updated soon - we have
    distance posteriors... Just need to understand how GWEMOPT expects this to be tabulated.

    Note: filename is JUST the name without the path- the path is calculated in situ.
    """
    from astropy.table import Table, Column
    import os.path

    # get path and check filename...
    # TO DO: add this to a helper function to check paths of data directories are congruent.... Can also do this for tile functions in the same function but different to lisabeta function, but this function reference lisabeta one to harmonize the architectures..
    if len(SkyMapFileName.split("Skymap_files"))==1 and len(SkyMapFileName.split("SYNEX"))>1:
        raise ValueError("Sorry but for continuity across functions please direct skymap directories to be within './SYNEX/Skymap_files/...'")
    elif len(SkyMapFileName.split("Skymap_files"))>1 and len(SkyMapFileName.split("SYNEX"))==1:
        print("Directing given filename to './SYNEX/Skymap_files/...'")
        SkyMapFileName = SYNEX_PATH + "/Skymap_files" + SkyMapFileName.split("Skymap_files")[-1]
    elif len(SkyMapFileName.split("Skymap_files"))==1 and len(SkyMapFileName.split("SYNEX"))==1:
        print("WARNING: Adding ./SYNEX/Skymap_files/ to filename provided...")
        SkyMapFileName = SYNEX_PATH + "/Skymap_files/" + SkyMapFileName
    if SkyMapFileName[-5:]!='.fits':
        SkyMapFileName = SkyMapFileName+'.fits'

    # Check if directory exists and create if it doesnt
    SkyMapFilePath = "/".join(SkyMapFileName.split("/")[:-1])
    pathlib.Path(SkyMapFilePath).mkdir(parents=True, exist_ok=True)

    # Write to fits file
    if "distmu" in map_struct:
        data_to_save = np.vstack((map_struct["prob"],map_struct["distmu"],map_struct["distsigma"],map_struct["distnorm"]))
        if PermissionToWrite:
            hp.write_map(SkyMapFileName, data_to_save, overwrite=True)
        if go_params!=None:
            go_params["do3D"]=True
    elif PermissionToWrite:
        hp.write_map(SkyMapFileName, map_struct["prob"], overwrite=True)

    if go_params!=None:
        # update the filename stored in go_params
        go_params["skymap"] = SkyMapFileName
        return go_params,map_struct
    else:
        return SkyMapFileName,map_struct






#                      ##########################
#                      #                        #
#                      #   Plotting functions   #
#                      #                        #
#                      ##########################

##### Need to organize these better and rename #####

def PlotSkyMapData(source,SaveFig=False,plotName=None,DO_CONTOURS=False):
    """
    Plotting tool to either take a filename or a source object and plot the skyap associated.
    Adapted from 'gwemopt.plotting.skymap.py'
    """
    unit='Gravitational-wave probability'
    cbar=False

    lons = np.arange(-150.0,180,30.0) #
    lats = np.zeros(lons.shape)

    if np.percentile(source.map_struct["prob"],99) > 0:
        hp.mollview(source.map_struct["prob"],title='',unit=unit,cbar=cbar,min=np.percentile(source.map_struct["prob"],1),max=np.percentile(source.map_struct["prob"],99),cmap=cmap)
    else:
        hp.mollview(source.map_struct["prob"],title='',unit=unit,cbar=cbar,min=np.percentile(source.map_struct["prob"],1),cmap=cmap)

    hp.projplot(source.true_ra, source.true_dec, lonlat=True, coord='G', marker='D', markersize=1.7, c='blue', linestyle='None', label='true location')
    plt.legend()

    add_edges()
    if SaveFig:
        # Default save name
        if plotName==None:
            plotName = SYNEX_PATH+"/Plots"+source.sky_map.split("Skymap_files")[-1]
            plotName = plotName[:-5]+"_prob.pdf"
        # Check extension is there - move this to a general function please
        if plotName[-4:]!=".pdf":
            plotName = plotName + ".pdf"
        # Check if directory tree exists
        PlotPath="/".join(plotName.split("/")[:-1])
        pathlib.Path(PlotPath).mkdir(parents=True, exist_ok=True)
        # Save
        plt.savefig(plotName,dpi=200)
    else:
        plt.show()
    plt.close('all')

    # Extra plotting for 3D skymaps
    if "distmu" in source.map_struct:
        fin = np.copy(source.map_struct["distmu"])
        fin[~np.isfinite(fin)] = np.nan
        hp.mollview(source.map_struct["distmu"],unit='Distance [Mpc]',min=np.nanpercentile(fin,10),max=np.nanpercentile(fin,90))
        add_edges()
        if SaveFig:
            plotName = plotName.split("_prob")[0]+'_dist.pdf'
            plt.savefig(plotName,dpi=200)
        else:
            plt.show()
        plt.close('all')
    if "distmed" in source.map_struct:
        fin = np.copy(source.map_struct["distmed"])
        fin[~np.isfinite(fin)] = np.nan
        hp.mollview(source.map_struct["distmed"],unit='Distance [Mpc]',min=np.nanpercentile(fin,10),max=np.nanpercentile(fin,90))
        add_edges()
        if SaveFig:
            plotName = plotName.split("_dist")[0]+'_dist_median.pdf'
            plt.savefig(plotName,dpi=200)
        else:
            plt.show()
        plt.close('all')

def PlotTilesArea_old(TileFileName,n_tiles=10):
    """
    Function to plot a sample of tiles from a saved dictionary of tiles
    """
    # Load Tile dictionary
    with open(TileFileName, 'rb') as f:
        TileDict = pickle.load(f)
    from matplotlib.patches import Rectangle

    # Set the default color cycle
    TileColor=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
    n_tile_colours = 10

    # Get background data for liklihoods
    fig, InjParam_InjVals, SampleModes, X, lambda_bins, Y, beta_bins, Z = PlotInferenceLambdaBeta(TileDict["LISA Data File"], bins=50, SkyProjection=False, SaveFig=False, return_data=True)
    # fig = PlotInferenceLambdaBeta(TileDict["LISA Data File"], bins=50, SkyProjection=False, SaveFig=False, return_data=True)
    ax = plt.gca()

    # Add n_tiles tiles
    for ii in range(1,n_tiles+1):
        x_lower = TileDict[str(ii)]["lambda_range"][0]
        y_lower = TileDict[str(ii)]["beta_range"][0]
        x_upper = TileDict[str(ii)]["lambda_range"][1]
        y_upper = TileDict[str(ii)]["beta_range"][1]
        width = x_upper-x_lower
        height = y_upper-y_lower
        if ii<n_tile_colours:
            ax.add_patch(Rectangle((x_lower, y_lower), width, height, linewidth=1,
                             facecolor="none", edgecolor=TileColor[ii]))
        else:
            ax.add_patch(Rectangle((x_lower, y_lower), width, height, linewidth=1,
                             facecolor="none", edgecolor=TileColor[ii%n_tile_colours]))

    # Formatting stuff
    plt.xlabel("RA [rad]")
    plt.ylabel("Dec [rad]")
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi/2.,np.pi/2.])
    plt.show()

def PlotLikeRatioFoMFromJsonWithInferenceInlays(FoMJsonFileAndPath, BF_lim=20., SaveFig=False, InlayType="histogram"):
    """
    NB: This function has hard coded paths in it and will likely be removed wince we have moved on from needing it
    """

    with open(FoMJsonFileAndPath) as f:
        FoMOverRange = json.load(f)
    f.close()

    # Set the Bayes limit
    # For inc vs maxf: beta BF_lim=31, lambda BF_lim=15.5
    # For inc vs beta: beta BF_lim=15.5, lambda BF_lim= anything less than 11,000
    # For M vs z: beta BF_lim= between 10.5 and 31, lambda BF_lim= between 50 and several 1000

    # Get the variables
    Vars = FoMOverRange["LoopVariables"]

    # Check that you loaded a grided FoMOverRange file
    if not FoMOverRange["IsGrid"]:
        raise ValueError("Must be a full grid FoMOverRange file, otherwise this code does not work.")
    else:
        ################### PLOT THE MAIN DATA ###################

        X,Y = np.meshgrid(FoMOverRange[Vars[0]+"_xs"], FoMOverRange[Vars[1]+"_xs"])
        if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA" or FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
            Z = np.log10(abs(np.array(FoMOverRange["grid_ys"])+BF_lim))
        else:
            Z = np.log10(FoMOverRange["grid_ys"])

        # Master figure and axes
        fig, ax = plt.subplots(constrained_layout=True)
        im = ax.pcolormesh(X, Y, Z, shading='gouraud', vmin=Z.min(), vmax=Z.max())
        cbar = fig.colorbar(im, ax=ax)


        ################### DEFINE X-VARIABLE INLAYS ###################

        if Vars[0] == "M":
            ax.set_xscale('log')
            plt.xlabel(r'M$_{\mathrm{tot}} \; (M_{\odot})$')
            VarStringForSaveFile = "Mass_"

            # Positioning arguments for inlay inference plots
            pos_lefts = np.array([0.13, 0.36, 0.61])
            pos_bottoms = np.array([0.42, 0.45, 0.47])
            pos_widths = np.array([0.2, 0.2, 0.2])
            pos_heights = np.array([0.2, 0.2, 0.2])

            # Inference plot files and strings
            VarsTestedStrs = ["3e6", "1e7", "8e7"]
            DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/MassGridTest_NoMpi_"
            jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/MassGridTest_NoMpi_"

            # Arrow coords and widths
            Arrow_xs = [3.e6, 1.e7, 8e7]
            Arrow_ys = [3., 3., 3.]
            Arrow_dxs = [2.e6-Arrow_xs[0], 2.e7-Arrow_xs[1], 9.e7-Arrow_xs[2]]
            Arrow_dys = [5.-Arrow_ys[0], 5.6-Arrow_ys[1], 4.2-Arrow_ys[2]]

            # Record the coord of the inferrence point
            ax.plot(Arrow_xs, Arrow_ys, 'ro', linestyle='None')
        elif Vars[0] == "inc":
            ax.set_xlim([0.,np.pi])
            plt.xlabel(r'$\iota \; (\mathrm{rad.})$')
            VarStringForSaveFile = "Inc_"

            # Inference plot files and strings
            VarsTestedStrs = ["0", "3PiBy6", "Pi"]
            DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/IncGridTest_NoMpi_"
            jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/IncGridTest_NoMpi_"

            # Arrow coords and widths
            if Vars[1] == "beta":
                # Positioning arguments for inlay inference plots
                pos_lefts = np.array([0.13, 0.36, 0.61])
                pos_bottoms = np.array([0.42, 0.45, 0.47])
                pos_widths = np.array([0.2, 0.2, 0.2])
                pos_heights = np.array([0.2, 0.2, 0.2])

                Arrow_xs = [0.000000002, 3.*np.pi/6., np.pi-0.000000002]
                Arrow_ys = [-0.6785264548614762, -0.6785264548614762, -0.6785264548614762] # [-3.*np.pi/8., -3.*np.pi/8., -3.*np.pi/8.]
                Arrow_dxs = [0.2-Arrow_xs[0], 1.8-Arrow_xs[1], 3.-Arrow_xs[2]]
                Arrow_dys = [-0.45-Arrow_ys[0], -0.41-Arrow_ys[1], -0.32-Arrow_ys[2]]
            elif Vars[1] == "maxf":
                # Positioning arguments for inlay inference plots
                pos_lefts = np.array([0.13, 0.36, 0.61])
                pos_bottoms = np.array([0.71, 0.71, 0.71])
                pos_widths = np.array([0.2, 0.2, 0.2])
                pos_heights = np.array([0.2, 0.2, 0.2])

                Arrow_xs = [0.000000002, 3.*np.pi/6., np.pi-0.000000002]
                Arrow_ys = [0.5, 0.5, 0.5]
                Arrow_dxs = [0.2-Arrow_xs[0], 0., 3.-Arrow_xs[2]] # 1.8-Arrow_xs[1]
                Arrow_dys = [0.2-Arrow_ys[0], 0.2-Arrow_ys[1], 0.2-Arrow_ys[2]]

            # Record the coord of the inferrence point
            ax.plot(Arrow_xs, Arrow_ys, 'ro', linestyle='None')
        elif Vars[0] == "maxf":
            ax.set_xscale('log')
            plt.xlabel(r'$\mathrm{f}_{\mathrm{max}} \; (\mathrm{Hz})$')
            VarStringForSaveFile = "Fmax_"
        elif Vars[0] == "beta":
            plt.xlabel(r'$\beta_{SSB} \; (\mathrm{rad.})$')
            VarStringForSaveFile = "Beta_"
        elif Vars[0] == "z":
            plt.xlabel(r'z')
            VarStringForSaveFile = "Red_"

        ################### ADD X-VARIABLE INLAYS ###################
        # Add arrows
        for iArrow in range(len(Arrow_xs)):
            mpl.pyplot.arrow(Arrow_xs[iArrow], Arrow_ys[iArrow], Arrow_dxs[iArrow], Arrow_dys[iArrow], width=0.000001)

        PostVals = {}
        for i in range(len(VarsTestedStrs)):
            DataFile = DataFilePath + VarsTestedStrs[i] + ".h5"
            jsonFile = jsonFilePath + VarsTestedStrs[i] + ".json"
            DataFileRaw = DataFile[:-3]+"_raw.h5"

            # Make the inference plots from scratch since we don't pass back axes handles in the utils functions
            DataFileLocAndName = DataFile
            hist_n_bins=1000
            [infer_params, inj_param_vals, static_params, meta_data] = read_h5py_file(DataFileLocAndName)
            # print(inj_param_vals["source_params_Lframe"]["beta"][0], inj_param_vals["source_params_SSBframe"]["beta"][0])
            if not inj_param_vals: # Raw files at first didn't record this, so make sure it's there...
                # Get the inj values from the processed data file instead
                DataFileLocAndName_NotRaw = DataFileLocAndName[:-7] + ".h5"
                [_,inj_param_vals,_,_] = read_h5py_file(DataFileLocAndName_NotRaw)
            ndim = len(infer_params.keys())
            labels = list(infer_params.keys())
            if np.size(infer_params[labels[0]][0])>1:
                nsamples = len(infer_params[labels[0]][0])
            else:
                nsamples = len(infer_params[labels[0]])
            print("Posterior sample length: " + str(nsamples) + ", number of infered parameters: " + str(ndim))
            # Grab data for infered parameters
            data = np.empty([ndim,nsamples])
            SampleModes = []
            for ii in range(ndim):
                if np.size(infer_params[labels[0]][0])>1:
                    data[ii][:] = infer_params[labels[ii]][0]
                else:
                    data[ii][:] = infer_params[labels[ii]]
                histn,histbins = np.histogram(data[ii,:], bins=hist_n_bins)
                histbins = histbins[:-1] + 0.5*(histbins[2]-histbins[1]) # change from left bins edges to middle of bins
                mode = histbins[histn==max(histn)]
                if len(mode)>1:
                    mode = mode[1] # Doe now take the first on the list, but if there are several we need to work out what to do there...
                SampleModes.append(mode)
            data = np.transpose(np.array(data)) # should have shape [nsamples, ndim]

            # Get injected values
            InjParam_InjVals = []
            for key in infer_params.keys():
                InjParam_InjVals.append(inj_param_vals["source_params_Lframe"][key][0]) # Lframe is right.

            if InlayType=="scatter":
                # Scatter plot of labmda and beta posterior chains, marginalized over all other params
                ax_inlay = plt.axes([pos_lefts[i], pos_bottoms[i], pos_widths[i], pos_heights[i]], facecolor='w')
                ax_inlay.set_facecolor('yellow')
                ax_inlay.scatter(data[:,0], data[:,5], marker=".", linestyle="None")

                # Add injected and mode vertical and horizontal lines
                ax_inlay.axvline(InjParam_InjVals[0], color="green", linestyle=":")
                ax_inlay.axvline(SampleModes[0], color="blue", linestyle=":")
                ax_inlay.axhline(InjParam_InjVals[5], color="green", linestyle=":")
                ax_inlay.axhline(SampleModes[5], color="blue", linestyle=":")

                # Add points at injected and mode values
                ax_inlay.plot(InjParam_InjVals[0], InjParam_InjVals[5], "sg")
                ax_inlay.plot(SampleModes[0], SampleModes[5], "sb")

                # Labels and grid
                plt.grid()
                ax_inlay.set_xlabel(r'$\beta_{L}$') # (labels[0])
                ax_inlay.set_ylabel(r'$\lambda_{L}$') # (labels[5])
            elif InlayType=="histogram":
                # Histogram of labmda or beta posterior chains, marginalized over all other params
                ax_inlay = plt.axes([pos_lefts[i], pos_bottoms[i], pos_widths[i], pos_heights[i]], facecolor='w')
                ax_inlay.set_facecolor('yellow')
                if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA":
                    iiData = 0
                elif FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
                    iiData = 5
                ax_inlay.hist(data[:,iiData], 100)

                # Add injected and mode vertical lines
                ax_inlay.axvline(InjParam_InjVals[iiData], color="green", linestyle=":")
                ax_inlay.axvline(SampleModes[iiData], color="blue", linestyle=":")

                # Labels
                if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA":
                    ax_inlay.set_xlabel(r'$\beta_{L}$') # (labels[0])
                elif FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
                    ax_inlay.set_xlabel(r'$\lambda_{L}$') # (labels[5])

            # Move to position
            plt.xticks([])
            plt.yticks([])

        ################### DEFINE Y-VARIABLE INLAYS ###################

        if Vars[1] == "M":
            ax.set_yscale('log')
            ax.set_ylabel(r'M$_{\mathrm{tot}} \; (M_{\odot})$')
            VarStringForSaveFile += "Mass"
        elif Vars[1] == "inc":
            ax.set_ylabel(r'$\iota \; (\mathrm{rad.})$')
            VarStringForSaveFile += "Iota"
        elif Vars[1] == "maxf":
            ax.set_yscale('log')
            ax.set_ylabel(r'$\mathrm{f}_{\mathrm{max}} \; (\mathrm{Hz})$')
            VarStringForSaveFile += "Fmax"

            # Positioning arguments for inlay inference plots
            pos_lefts = np.array([0.61, 0.39, 0.17])
            pos_bottoms = np.array([0.18, 0.29, 0.41])
            pos_widths = np.array([0.2, 0.2, 0.2])
            pos_heights = np.array([0.2, 0.2, 0.2])

            # Inference plot files and strings
            VarsTestedStrs = ["4eM4", "2eM3", "3eM2"]
            DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/maxfGridTest_NoMpi_"
            jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/maxFGridTest_NoMpi_"

            # Arrow coords and widths
            Arrow_xs = [np.pi/10., np.pi/10., np.pi/10.]
            Arrow_ys = [4.e-4, 2.e-3, 3.e-2] # [-3.*np.pi/8., 3.*np.pi/8.]
            Arrow_dxs = [2.2-Arrow_xs[0], 1.2-Arrow_xs[1], 0.3-Arrow_xs[2]] # [2.*np.pi/10.-Arrow_xs[0], 2.5*np.pi/10.-Arrow_xs[1], 2.5*np.pi/10.-Arrow_xs[1]]
            Arrow_dys = [0.0003, 0.0005, -0.02] # [0.0005-Arrow_ys[0], 0.2-Arrow_ys[1]]

            # Record the coord of the inferrence point
            ax.plot(Arrow_xs, Arrow_ys, 'ro', linestyle='None')
        elif Vars[1] == "beta":
            ax.set_ylabel(r'$\beta_{SSB} \; (\mathrm{rad.})$')
            VarStringForSaveFile += "Beta"

            # Positioning arguments for inlay inference plots
            pos_lefts = np.array([0.25, 0.3])
            pos_bottoms = np.array([0.18, 0.75])
            pos_widths = np.array([0.2, 0.2])
            pos_heights = np.array([0.2, 0.2])

            # Inference plot files and strings
            VarsTestedStrs = ["Neg3PiBy8", "3PiBy8"]
            DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/BetaGridTest_NoMpi_"
            jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/BetaGridTest_NoMpi_"

            # Arrow coords and widths
            Arrow_xs = [np.pi/10., np.pi/10.]
            Arrow_ys = [-0.6785264548614762, 0.30074618202252795] # [-3.*np.pi/8., 3.*np.pi/8.]
            Arrow_dxs = [2.*np.pi/10.-Arrow_xs[0], 2.5*np.pi/10.-Arrow_xs[1]]
            Arrow_dys = [-2.*np.pi/8.-Arrow_ys[0], 2.*np.pi/8.-Arrow_ys[1]]

            # Record the coord of the inferrence point
            ax.plot(Arrow_xs, Arrow_ys, 'ro', linestyle='None')
        elif Vars[1] == "z":
            ax.set_ylabel(r'z')
            VarStringForSaveFile += "Red"

            # Positioning arguments for inlay inference plots
            pos_lefts = np.array([0.38, 0.38])
            pos_bottoms = np.array([0.18, 0.75])
            pos_widths = np.array([0.2, 0.2])
            pos_heights = np.array([0.2, 0.2])

            # Inference plot files and strings
            VarsTestedStrs = ["2", "8"]
            DataFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_data/RedGridTest_NoMpi_"
            jsonFilePath = "/Users/jonathonbaird/Documents/LabEx_PostDoc/SYNEX/inference_param_files/RedGridTest_NoMpi_"

            # Arrow coords and widths
            Arrow_xs = [6.e6, 6.e6]
            Arrow_ys = [2., 8.]
            Arrow_dxs = [1.3e7-Arrow_xs[0], 1.5e7-Arrow_xs[1]]
            Arrow_dys = [1.48-Arrow_ys[0], 8.75-Arrow_ys[1]]

            # Record the coord of the inferrence point
            ax.plot(Arrow_xs, Arrow_ys, 'ro', linestyle='None')
        if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA":
            cbar.set_label(r'$\log_{10}(|\mathcal{B}_{\beta}|)$') # (r'$\log_{10}(|\Sigma_{b=0}^{3}\log_e(L_{(-1,b)}) - \log_e(L_{(1,b)})+20|)$')
        elif FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
            cbar.set_label(r'$\log_{10}(|\mathcal{B}_{\lambda}|)$') # r'$\log_{10}(|\Sigma_{a=1,-1}\log_e(L_{(a,1)}) + \log_e(L_{(a,2)}) + \log_e(L_{(a,3)}) - \log_e(L_{(a,0)})+20|)$')
        else:
            cbar.set_label(r'$\log_{10}(\Delta \Omega \; (\mathrm{sq. deg.}))$')

        ################### ADD Y-VARIABLE INLAYS ###################
        # Set active axes back to master axes
        plt.sca(ax)

        # Add arrows
        for iArrow in range(len(Arrow_xs)):
            mpl.pyplot.arrow(Arrow_xs[iArrow], Arrow_ys[iArrow], Arrow_dxs[iArrow], Arrow_dys[iArrow], Figure=fig, width=0.000001)

        PostVals = {}
        for i in range(len(VarsTestedStrs)):
            DataFile = DataFilePath + VarsTestedStrs[i] + ".h5"
            jsonFile = jsonFilePath + VarsTestedStrs[i] + ".json"
            DataFileRaw = DataFile[:-3]+"_raw.h5"

            # Make the inference plots from scratch since we don't pass back axes handles in the utils functions
            DataFileLocAndName = DataFile
            hist_n_bins=1000
            [infer_params, inj_param_vals, static_params, meta_data] = read_h5py_file(DataFileLocAndName)
            # print(inj_param_vals["source_params_Lframe"]["beta"][0], inj_param_vals["source_params_SSBframe"]["beta"][0])
            # print(z_at_value(Planck13.distmod, inj_param_vals["source_params_Lframe"]["dist"][0]))
            # from astropy.cosmology import WMAP9 as cosmo
            # print(cosmo.luminosity_distance(3.).to("Mpc").value, inj_param_vals["source_params_SSBframe"]["dist"][0])
            if not inj_param_vals: # Raw files at first didn't record this, so make sure it's there...
                # Get the inj values from the processed data file instead
                DataFileLocAndName_NotRaw = DataFileLocAndName[:-7] + ".h5"
                [_,inj_param_vals,_,_] = read_h5py_file(DataFileLocAndName_NotRaw)
            ndim = len(infer_params.keys())
            labels = list(infer_params.keys())
            if np.size(infer_params[labels[0]][0])>1:
                nsamples = len(infer_params[labels[0]][0])
            else:
                nsamples = len(infer_params[labels[0]])
            print("Posterior sample length: " + str(nsamples) + ", number of infered parameters: " + str(ndim))
            # Grab data for infered parameters
            data = np.empty([ndim,nsamples])
            SampleModes = []
            for ii in range(ndim):
                if np.size(infer_params[labels[0]][0])>1:
                    data[ii][:] = infer_params[labels[ii]][0]
                else:
                    data[ii][:] = infer_params[labels[ii]]
                histn,histbins = np.histogram(data[ii,:], bins=hist_n_bins)
                histbins = histbins[:-1] + 0.5*(histbins[2]-histbins[1]) # change from left bins edges to middle of bins
                mode = histbins[histn==max(histn)]
                if len(mode)>1:
                    mode = mode[1] # Doe now take the first on the list, but if there are several we need to work out what to do there...
                SampleModes.append(mode)
            data = np.transpose(np.array(data)) # should have shape [nsamples, ndim]

            # Get injected values
            InjParam_InjVals = []
            for key in infer_params.keys():
                InjParam_InjVals.append(inj_param_vals["source_params_Lframe"][key][0]) # Lframe is right.

            if InlayType=="scatter":
                # Sctter plot of labmda and beta posterior chains, marginalized over all other params
                ax_inlay = plt.axes([pos_lefts[i], pos_bottoms[i], pos_widths[i], pos_heights[i]], facecolor='w')
                ax_inlay.set_facecolor('yellow')
                ax_inlay.scatter(data[:,0], data[:,5], marker=".", linestyle="None")

                # Add injected and mode vertical and horizontal lines
                ax_inlay.axvline(InjParam_InjVals[0], color="green", linestyle=":")
                ax_inlay.axvline(SampleModes[0], color="blue", linestyle=":")
                ax_inlay.axhline(InjParam_InjVals[5], color="green", linestyle=":")
                ax_inlay.axhline(SampleModes[5], color="blue", linestyle=":")

                # Add points at injected and mode values
                ax_inlay.plot(InjParam_InjVals[0], InjParam_InjVals[5], "sg")
                ax_inlay.plot(SampleModes[0], SampleModes[5], "sb")

                # Labels and grid
                plt.grid()
                ax_inlay.set_xlabel(r'$\beta_{L}$') # (labels[0])
                ax_inlay.set_ylabel(r'$\lambda_{L}$') # (labels[5])

            elif InlayType=="histogram":
                # histogram of labmda or beta posterior chains, marginalized over all other params
                ax_inlay = plt.axes([pos_lefts[i], pos_bottoms[i], pos_widths[i], pos_heights[i]], facecolor='w')
                ax_inlay.set_facecolor('yellow')
                if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA":
                    iiData = 0
                elif FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
                    iiData = 5
                ax_inlay.hist(data[:,iiData], 100)

                # Add injected and mode vertical lines
                ax_inlay.axvline(InjParam_InjVals[iiData], color="green", linestyle=":")
                ax_inlay.axvline(SampleModes[iiData], color="blue", linestyle=":")

                # Labels
                if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA":
                    ax_inlay.set_xlabel(r'$\beta_{L}$') # (labels[0])
                elif FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
                    ax_inlay.set_xlabel(r'$\lambda_{L}$') # (labels[5])

            # Move to position
            plt.xticks([])
            plt.yticks([])

    # Save?
    if SaveFig:
        if FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
            plt.savefig("/Users/jonathonbaird/Documents/LabEx_PostDoc/Figs/lnL_skymodes_LambdaReflections_" + VarStringForSaveFile + "_inlays_BF_" + str(BF_lim) + ".png", facecolor='w', transparent=False)
        elif FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA":
            plt.savefig("/Users/jonathonbaird/Documents/LabEx_PostDoc/Figs/lnL_skymodes_BetaReflections_" + VarStringForSaveFile + "_inlays_BF_" + str(BF_lim) + ".png", facecolor='w', transparent=False)
    plt.show()

def PlotLikeRatioFoMFromJsonWithAllInferencePoints(FoMJsonFileAndPath, BF_lim=20., ModeJumpLimit=0.1, SaveFig=False):
    """
    NB: This function has hard coded paths in it and will likely be removed wince we have moved on from needing it
    """
    # Load data
    with open(FoMJsonFileAndPath) as f:
        FoMOverRange = json.load(f)
    f.close()

    # Set the Bayes limit
    # For inc vs maxf: beta BF_lim=31, lambda BF_lim=15.5
    # For inc vs beta: beta BF_lim=15.5, lambda BF_lim= anything less than 11,000
    # For M vs z: beta BF_lim= between 10.5 and 31, lambda BF_lim= between 50 and several 1000

    # Get the variables
    Vars = FoMOverRange["LoopVariables"]

    # Check that you loaded a grided FoMOverRange file
    if not FoMOverRange["IsGrid"]:
        raise ValueError("Must be a full grid FoMOverRange file, otherwise this code does not work.")
    else:
        ################### PLOT THE MAIN DATA ###################

        X,Y = np.meshgrid(FoMOverRange[Vars[0]+"_xs"], FoMOverRange[Vars[1]+"_xs"])
        if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA" or FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
            Z = np.log10(abs(np.array(FoMOverRange["grid_ys"])+BF_lim))
        else:
            Z = np.log10(FoMOverRange["grid_ys"])

        # Master figure and axes
        fig, ax = plt.subplots(constrained_layout=True)
        im = ax.pcolormesh(X, Y, Z, shading='gouraud', vmin=Z.min(), vmax=Z.max())
        cbar = fig.colorbar(im, ax=ax)

        # Get the path for lisabeta to coordinate search for inference data
        import os
        LisaBetaPath = os.path.dirname(os.path.realpath(__file__))
        InfDataPath = LisaBetaPath + "/../inference_data/"
        jsonFilesPath = LisaBetaPath+"/../inference_param_files/"

        ################### DEFINE X-VARIABLE INLAYS ###################
        import glob
        if Vars[0] == "M":
            ax.set_xscale('log')
            plt.xlabel(r'M$_{\mathrm{tot}} \; (M_{\odot})$')
            VarStringForSaveFile = "Mass_"

            # Get the values tested
            GridJsons = glob.glob(jsonFilesPath+"MassGridTest_NoMpi_*.json")
        elif Vars[0] == "inc":
            ax.set_xlim([0.,np.pi])
            plt.xlabel(r'$\iota \; (\mathrm{rad.})$')
            VarStringForSaveFile = "Inc_"

            # Get the values tested
            GridJsons = glob.glob(jsonFilesPath+"IncGridTest_NoMpi_*.json")
        elif Vars[0] == "maxf":
            ax.set_xscale('log')
            plt.xlabel(r'$\mathrm{f}_{\mathrm{max}} \; (\mathrm{Hz})$')
            VarStringForSaveFile = "Fmax_"

            # Get the values tested
            GridJsons = glob.glob(jsonFilesPath+"maxfGridTest_NoMpi_*.json")
        elif Vars[0] == "beta":
            plt.xlabel(r'$\beta_{SSB} \; (\mathrm{rad.})$')
            VarStringForSaveFile = "Beta_"

            # Get the values tested
            GridJsons = glob.glob(jsonFilesPath+"BetaGridTest_NoMpi_*.json")
        elif Vars[0] == "z":
            plt.xlabel(r'z')
            VarStringForSaveFile = "Red_"

            # Get the values tested
            GridJsons = glob.glob(jsonFilesPath+"RedGridTest_NoMpi_*.json")

        ################### DEFINE Y-VARIABLE Points ###################

        if Vars[1] == "M":
            ax.set_yscale('log')
            ax.set_ylabel(r'M$_{\mathrm{tot}} \; (M_{\odot})$')
            VarStringForSaveFile += "Mass"

            # Get the values tested
            [GridJsons.append(FileString) for FileString in glob.glob(jsonFilesPath+"MassGridTest_NoMpi_*.json")]
        elif Vars[1] == "inc":
            ax.set_ylabel(r'$\iota \; (\mathrm{rad.})$')
            VarStringForSaveFile += "Iota"

            # Get the values tested
            [GridJsons.append(FileString) for FileString in glob.glob(jsonFilesPath+"IncGridTest_NoMpi_*.json")]
        elif Vars[1] == "maxf":
            ax.set_yscale('log')
            ax.set_ylabel(r'$\mathrm{f}_{\mathrm{max}} \; (\mathrm{Hz})$')
            VarStringForSaveFile += "Fmax"

            # Get the values tested
            [GridJsons.append(FileString) for FileString in glob.glob(jsonFilesPath+"maxfGridTest_NoMpi_*.json")]
        elif Vars[1] == "beta":
            ax.set_ylabel(r'$\beta_{SSB} \; (\mathrm{rad.})$')
            VarStringForSaveFile += "Beta"

            # Get the values tested
            [GridJsons.append(FileString) for FileString in glob.glob(jsonFilesPath+"BetaGridTest_NoMpi_*.json")]
        elif Vars[1] == "z":
            ax.set_ylabel(r'z')
            VarStringForSaveFile += "Red"

            # Get the values tested
            [GridJsons.append(FileString) for FileString in glob.glob(jsonFilesPath+"RedGridTest_NoMpi_*.json")]
        if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA":
            cbar.set_label(r'$\log_{10}(|\mathcal{B}_{\beta}|)$') # (r'$\log_{10}(|\Sigma_{b=0}^{3}\log_e(L_{(-1,b)}) - \log_e(L_{(1,b)})+20|)$')
        elif FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
            cbar.set_label(r'$\log_{10}(|\mathcal{B}_{\lambda}|)$') # r'$\log_{10}(|\Sigma_{a=1,-1}\log_e(L_{(a,1)}) + \log_e(L_{(a,2)}) + \log_e(L_{(a,3)}) - \log_e(L_{(a,0)})+20|)$')
        else:
            cbar.set_label(r'$\log_{10}(\Delta \Omega \; (\mathrm{sq. deg.}))$')


        ################### Get the values of each point requested ###################

        VarsTestedStrs = [(GridJson.split('_')[-1]).split('.')[0] for GridJson in GridJsons]
        Gridh5s = [str(InfDataPath + '_'.join((GridJson.split('/')[-1]).split('_')[:-1]) + "_" + (GridJson.split('_')[-1]).split('.')[0] + ".h5") for GridJson in GridJsons]
        TestedVals_xs = []
        TestedVals_ys = []
        PointColours = []
        for iifile in range(len(GridJsons)):
            # Load rach json file that stored param values of each tested point
            with open(GridJsons[iifile]) as f:
                GridJsonData = json.load(f)
            f.close()

            #Grab the right data depending on which parameters are being plotted
            if Vars[0]=="M":
                TestedVals_xs.append(GridJsonData["source_params"]["m1"]+GridJsonData["source_params"]["m2"])
            elif Vars[0]=="z":
                # TestedVals_xs.append(GridJsonData["source_params"]["dist"])
                TestedVals_xs.append(z_at_value(cosmo.luminosity_distance, GridJsonData["source_params"]["dist"]*u.Mpc))
            elif Vars[0]=="maxf":
                TestedVals_xs.append(GridJsonData["waveform_params"]["maxf"])
            else:
                TestedVals_xs.append(GridJsonData["source_params"][Vars[0]])

            if Vars[1]=="M":
                TestedVals_ys.append(GridJsonData["source_params"]["m1"]+GridJsonData["source_params"]["m2"])
            elif Vars[1]=="z":
                # TestedVals_ys.append(GridJsonData["source_params"]["dist"])
                TestedVals_ys.append(z_at_value(cosmo.luminosity_distance, GridJsonData["source_params"]["dist"]*u.Mpc))
            elif Vars[1]=="maxf":
                TestedVals_ys.append(GridJsonData["waveform_params"]["maxf"])
            else:
                TestedVals_ys.append(GridJsonData["source_params"][Vars[1]])

            # Grab an estimate of the posterior modes to check if they are degenerate in either lambda or beta
            PosteriorStats = GetPosteriorStats(Gridh5s[iifile], ModeJumpLimit=ModeJumpLimit) # LogLikeJumpLimit=BF_lim) #

            if len(PosteriorStats["beta"]["mode"]) == 1 and len(PosteriorStats["lambda"]["mode"]) == 1:
                PointColours.append("black")
            elif len(PosteriorStats["beta"]["mode"]) > 1 and len(PosteriorStats["lambda"]["mode"]) == 1:
                PointColours.append("blue")
            elif len(PosteriorStats["beta"]["mode"]) == 1 and len(PosteriorStats["lambda"]["mode"]) > 1:
                PointColours.append("green")
            elif len(PosteriorStats["beta"]["mode"]) > 1 and len(PosteriorStats["lambda"]["mode"]) > 1:
                PointColours.append("red")


        ################### ADD points ###################
        from itertools import compress
        iiToScatter = [PointColour=="black" for PointColour in PointColours]
        ax.scatter(list(compress(TestedVals_xs, iiToScatter)), list(compress(TestedVals_ys, iiToScatter)),
                    c="black", marker='o', label=r'1-mode: No ref.')
        iiToScatter = [PointColour=="blue" for PointColour in PointColours]
        ax.scatter(list(compress(TestedVals_xs, iiToScatter)), list(compress(TestedVals_ys, iiToScatter)),
                    c="blue", marker='o', label=r'2-mode: $\beta$ ref.')
        iiToScatter = [PointColour=="green" for PointColour in PointColours]
        ax.scatter(list(compress(TestedVals_xs, iiToScatter)), list(compress(TestedVals_ys, iiToScatter)),
                    c="green", marker='o', label=r'4-mode: $\lambda$ ref.')
        iiToScatter = [PointColour=="red" for PointColour in PointColours]
        ax.scatter(list(compress(TestedVals_xs, iiToScatter)), list(compress(TestedVals_ys, iiToScatter)),
                    c="red", marker='o', label=r'8-mode: $\beta$ and $\lambda$ ref.')
        ax.legend()

    # Save?
    if SaveFig:
        if FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
            plt.savefig("/Users/jonathonbaird/Documents/LabEx_PostDoc/Figs/lnL_skymodes_LambdaReflections_" + VarStringForSaveFile + "_Points_BF_" + str(BF_lim) + ".png", facecolor='w', transparent=False)
        elif FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA":
            plt.savefig("/Users/jonathonbaird/Documents/LabEx_PostDoc/Figs/lnL_skymodes_BetaReflections_" + VarStringForSaveFile + "_Points_BF_" + str(BF_lim) + ".png", facecolor='w', transparent=False)
    plt.show()

def PlotHistsLambdaBeta(FileName, SaveFig=False, ParamToPlot="beta"): # "lambda" # ["beta", "lambda"]
    """
    Plotting function to output corner plots of inference data.
    Need to add the post-processing stuff from Sylvain's example online?
    """
    # Convert the requested params to a list if it's not already one
    if not type(ParamToPlot) == list:
        ParamToPlot = [ParamToPlot]

    # Unpack the data
    [infer_params, inj_param_vals, static_params, meta_data] = read_h5py_file(FileName)
    if not inj_param_vals: # Raw files at first didn't record this, so make sure it's there...
        # Get the inj values from the processed data file instead
        DataFileLocAndName_NotRaw = DataFileLocAndName[:-7] + ".h5"
        [_,inj_param_vals,_,_] = read_h5py_file(DataFileLocAndName_NotRaw)
    ndim = len(infer_params.keys())
    labels = list(infer_params.keys())
    if np.size(infer_params[labels[0]][0])>1:
        nsamples = len(infer_params[labels[0]][0])
    else:
        nsamples = len(infer_params[labels[0]])
    print("Posterior sample length: " + str(nsamples) + ", number of infered parameters: " + str(ndim))
    # Grab data for infered parameters
    data = np.empty([ndim,nsamples])
    SampleModes = []
    for ii in range(ndim):
        if np.size(infer_params[labels[0]][0])>1:
            data[ii][:] = infer_params[labels[ii]][0]
        else:
            data[ii][:] = infer_params[labels[ii]]
        histn,histbins = np.histogram(data[ii,:], bins=hist_n_bins)
        histbins = histbins[:-1] + 0.5*(histbins[2]-histbins[1]) # change from left bins edges to middle of bins
        mode = histbins[histn==max(histn)]
        if len(mode)>1:
            mode = mode[1] # Doe now take the first on the list, but if there are several we need to work out what to do there...
        SampleModes.append(mode)
    data = np.transpose(np.array(data)) # should have shape [nsamples, ndim]

    # Get injected values
    InjParam_InjVals = []
    for key in infer_params.keys():
        InjParam_InjVals.append(inj_param_vals["source_params_Lframe"][key][0]) # Lframe is right.

    # Sctter plot of lambda and beta posterior chains, marginalized over all other params
    for PlotParam in ParamToPlot:
        print(PlotParam)
        if PlotParam == "beta":
            iiParam = 0
        elif PlotParam == "lambda":
            iiParam = 5
        fig = plt.figure()
        ax = plt.gca()
        plt.hist(data[:,iiParam], 50)
        ax = plt.gca()

        # Add injected and mode vertical lines
        ax.axvline(InjParam_InjVals[iiParam], color="g", linestyle=":")
        ax.axvline(SampleModes[iiParam], color="r", linestyle=":")

        # Labels
        if PlotParam == "beta":
            plt.xlabel(r"$\beta$")
        elif PlotParam == "lambda":
            plt.xlabel(r"$\lambda$")

        # save the figure if asked
        if SaveFig:
            plt.savefig(FileNameAndPath[:-3]+'.png')

        plt.show()

def PlotInferenceData(FileName, SaveFig=False):
    """
    Plotting function to output corner plots of inference data.
    Need to add the post-processing stuff from Sylvain's example online?
    """
    import corner

    # Check filenames
    JsonFileLocAndName,H5FileLocAndName = CompleteLisabetaDataAndJsonFileNames(FileName)

    # Unpack the data
    [infer_params, inj_param_vals, static_params, meta_data] = read_h5py_file(DataFileLocAndName)
    if not inj_param_vals: # Raw files at first didn't record this, so make sure it's there...
        # Get the inj values from the processed data file instead
        DataFileLocAndName_NotRaw = DataFileLocAndName[:-7] + ".h5"
        [_,inj_param_vals,_,_] = read_h5py_file(DataFileLocAndName_NotRaw)
    ndim = len(infer_params.keys())
    labels = list(infer_params.keys())
    if np.size(infer_params[labels[0]][0])>1:
        nsamples = len(infer_params[labels[0]][0])
    else:
        nsamples = len(infer_params[labels[0]])
    print("Posterior sample length: " + str(nsamples) + ", number of infered parameters: " + str(ndim))
    # Grab data for infered parameters
    data = np.empty([ndim,nsamples])
    SampleModes = []
    for ii in range(ndim):
        if np.size(infer_params[labels[0]][0])>1:
            data[ii][:] = infer_params[labels[ii]][0]
        else:
            data[ii][:] = infer_params[labels[ii]]
        histn,histbins = np.histogram(data[ii,:], bins=hist_n_bins)
        histbins = histbins[:-1] + 0.5*(histbins[2]-histbins[1]) # change from left bins edges to middle of bins
        mode = histbins[histn==max(histn)]
        if len(mode)>1:
            mode = mode[1] # Doe now take the first on the list, but if there are several we need to work out what to do there...
        SampleModes.append(mode)
    data = np.transpose(np.array(data)) # should have shape [nsamples, ndim]

    # Get injected values
    InjParam_InjVals = []
    for key in infer_params.keys():
        InjParam_InjVals.append(inj_param_vals["source_params_Lframe"][key][0]) # Lframe is right.

    # Corner plot of posteriors
    figure = corner.corner(data, labels=labels,
                           quantiles=[0.16, 0.5, 0.84],
                           show_titles=True)

    # Extract the axes
    axes = np.array(figure.axes).reshape((ndim, ndim))

    # Loop over the diagonal
    for i in range(ndim):
        ax = axes[i, i]
        ax.axvline(InjParam_InjVals[i], color="g")
        ax.axvline(SampleModes[i], color="r")

    # Loop over the histograms
    for yi in range(ndim):
        for xi in range(yi):
            ax = axes[yi, xi]
            ax.axvline(InjParam_InjVals[xi], color="g")
            ax.axvline(SampleModes[xi], color="r")
            ax.axhline(InjParam_InjVals[yi], color="g")
            ax.axhline(SampleModes[yi], color="r")
            ax.plot(InjParam_InjVals[xi], InjParam_InjVals[yi], "sg")
            ax.plot(SampleModes[xi], SampleModes[yi], "sr")

    # save the figure if asked
    if SaveFig:
        # Put in folder for all lisabeta-related plots
        SaveFile = H5FileLocAndName[:-3]+'.png'
        pathlib.Path(SYNEX_PATH + "/Plots/lisabeta/").mkdir(parents=True, exist_ok=True)
        SaveFile = SYNEX_PATH + "/Plots/lisabeta/" + SaveFile.split("/")[-1]
        plt.savefig(SaveFile)

    plt.show()

def PlotInferenceLambdaBeta(FileName, bins=50, SkyProjection=False, SaveFig=False, return_data=False):
    """
    Plotting function to output corner plots of inference data.
    Need to add the post-processing stuff from Sylvain's example online?
    """
    # Update somee general properties
    # GeneralPlotFormatting()

    # Unpack the data
    [infer_params, inj_param_vals, static_params, meta_data] = read_h5py_file(FileName)
    ndim = len(infer_params.keys())
    labels = list(infer_params.keys())
    if np.size(infer_params[labels[0]][0])>1:
        nsamples = len(infer_params[labels[0]][0])
    else:
        nsamples = len(infer_params[labels[0]])
    print("Posterior sample length: " + str(nsamples) + ", number of infered parameters: " + str(ndim))
    # Grab data for infered parameters
    data = np.empty([ndim,nsamples])
    SampleModes = []
    for ii in range(ndim):
        if np.size(infer_params[labels[0]][0])>1:
            data[ii][:] = infer_params[labels[ii]][0]
        else:
            data[ii][:] = infer_params[labels[ii]]
        histn,histbins = np.histogram(data[ii,:], bins=hist_n_bins)
        histbins = histbins[:-1] + 0.5*(histbins[2]-histbins[1]) # change from left bins edges to middle of bins
        mode = histbins[histn==max(histn)]
        if len(mode)>1:
            mode = mode[1] # Doe now take the first on the list, but if there are several we need to work out what to do there...
        SampleModes.append(mode)
    data = np.transpose(np.array(data)) # should have shape [nsamples, ndim]

    # Get injected values
    InjParam_InjVals = []
    for key in infer_params.keys():
        InjParam_InjVals.append(inj_param_vals["source_params_Lframe"][key][0]) # Lframe is right.

    # Scatter plot of labmda and beta posterior chains, marginalized over all other params
    fig = plt.figure()
    ax = plt.gca()
    if SkyProjection:
        plt.subplot(111, projection="mollweide") # To match gwemopt projections

    # Get 2D hitogram data
    levels = [1.-0.997, 1.-0.95, 1.-0.9, 1.-0.68] # Plot contours for 1 sigma, 90% confidence, 2 sigma, and 3 sigma
    levels_labels=[str(int((1.-l)*1000)/10) for l in levels]

    # Manually do 2D histogram because I don't trust the ones I found online
    # data to hist = [data[:,0], data[:,5]] = beta, lambda. Beta is y data since it is declination.
    hist2D_pops = np.empty([bins,bins])
    areas = np.empty([bins,bins])
    bin_max = max(max(data[:,5]),SampleModes[5]+np.pi/10000.)
    bin_min = min(min(data[:,5]),SampleModes[5]-np.pi/10000.)
    if bin_max>np.pi:
        bin_max=np.pi
    if bin_min<-np.pi:
        bin_max=-np.pi
    if isinstance(bin_min,np.ndarray):
        bin_min = bin_min[0]
    if isinstance(bin_max,np.ndarray):
        bin_max = bin_max[0]
    lambda_bins = np.linspace(bin_min, bin_max, bins+1) # np.linspace(np.min(data[:,5]), np.max(data[:,5]), bins+1)
    bin_max = max(max(data[:,0]),SampleModes[0]+np.pi/20000.)
    bin_min = min(min(data[:,0]),SampleModes[0]-np.pi/20000.)
    if bin_max>np.pi/2.:
        bin_max=np.pi/2.
    if bin_min<-np.pi/2.:
        bin_max=-np.pi/2.
    if isinstance(bin_min,np.ndarray):
        bin_min = bin_min[0]
    if isinstance(bin_max,np.ndarray):
        bin_max = bin_max[0]
    beta_bins = np.linspace(bin_min, bin_max, bins+1) # np.linspace(np.min(data[:,0]), np.max(data[:,0]), bins+1)
    lambda_BinW = np.diff(lambda_bins)
    beta_BinW = np.diff(beta_bins)
    for xii in range(bins):
        list_len = list(range(len(data[:,5])))
        if xii == 0:
            values_ii = [valii for valii in list_len if (lambda_bins[xii]<=data[valii,5]<=lambda_bins[xii+1])]
        else:
            values_ii = [valii for valii in list_len if (lambda_bins[xii]<=data[valii,5]<lambda_bins[xii+1])]
        beta_pops, beta_bins = np.histogram(data[values_ii,0], bins=beta_bins) # beta_bins
        for yii in range(bins):
            hist2D_pops[yii,xii] = beta_pops[yii] # hist2D_pops[xii,yii] = beta_pops[yii]
            areas[yii,xii] = beta_BinW[yii]*lambda_BinW[xii]

    # Define bins information
    lambda_bins = lambda_bins[0:-1] + (lambda_bins[2]-lambda_bins[1] )*0.5
    beta_bins = beta_bins[0:-1] + (beta_bins[2]-beta_bins[1] )*0.5
    X,Y = np.meshgrid(lambda_bins, beta_bins)

    # Rescale bin populations to probabilities
    Z = hist2D_pops/areas
    Z = Z/bins # Z/np.max(Z)
    levels = [l*np.max(Z) for l in levels]

    # function for contour labeling
    labels_dict={}
    for l,s in zip(levels,levels_labels):
        labels_dict[l]=s

    # Plot contour
    contour = plt.contour(X, Y, Z, levels)
    plt.clabel(contour, colors ='k', fmt=labels_dict, fontsize=2)

    # Add injected and mode vertical and horizontal lines
    if not SkyProjection:
        ax.axhline(InjParam_InjVals[0], color="r", linestyle=":")
        ax.axhline(SampleModes[0], color="b", linestyle=":")
        ax.axvline(InjParam_InjVals[5], color="r", linestyle=":")
        ax.axvline(SampleModes[5], color="b", linestyle=":")

    # Add points at injected and mode values
    ax.plot(InjParam_InjVals[5], InjParam_InjVals[0], "sr")
    ax.plot(SampleModes[5], SampleModes[0], "sb")

    # Labels
    if not SkyProjection:
        plt.xlabel(labels[5]) # Lambda
        plt.ylabel(labels[0]) # beta

    # show now or return?
    if not return_data:
        plt.show()
        plt.grid()

        # save the figure if asked
        if SaveFig:
            plt.savefig(FileNameAndPath[:-3]+'.png')
    else:
        return fig, InjParam_InjVals, SampleModes, X, lambda_bins, Y, beta_bins, Z

def PlotLikeRatioFoMFromJson(FoMJsonFileAndPath, BF_lim=20., SaveFig=False):
    """
    Function to eventually replace the base plot functions in other plot util.
    Return an axis and figure handle that has the base colour plot for the
    lambda and beta sky mode replection using log(bayes factor) calculations only.

    NB: This function has hard coded paths in it and will likely be removed wince we have moved on from needing it
    """

    with open(FoMJsonFileAndPath) as f:
        FoMOverRange = json.load(f)
    f.close()

    # Set the Bayes limit
    # For inc vs maxf: beta BF_lim=31, lambda BF_lim=15.5
    # For inc vs beta: beta BF_lim=15.5, lambda BF_lim= anything less than 11,000
    # For M vs z: beta BF_lim= between 10.5 and 31, lambda BF_lim= between 50 and several 1000

    # Get the variables
    Vars = FoMOverRange["LoopVariables"]

    # Check that you loaded a grided FoMOverRange file
    if not FoMOverRange["IsGrid"]:
        raise ValueError("Must be a full grid FoMOverRange file, otherwise this code does not work.")
    else:

        ################### PLOT THE MAIN DATA ###################

        X,Y = np.meshgrid(FoMOverRange[Vars[0]+"_xs"], FoMOverRange[Vars[1]+"_xs"])

        if Vars[0] == "T_obs_end_to_merger":
            X = X/(60.*60.)
        if Vars[1] == "T_obs_end_to_merger":
            Y = Y/(60.*60.)

        if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA" or FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
            Z = np.log10(abs(np.array(FoMOverRange["grid_ys"])+BF_lim))
        else:
            Z = np.log10(FoMOverRange["grid_ys"])

        # Master figure and axes
        fig, ax = plt.subplots(constrained_layout=True)
        im = ax.pcolormesh(X, Y, Z, shading='gouraud', vmin=Z.min(), vmax=Z.max())
        cbar = fig.colorbar(im, ax=ax)

        if Vars[0] == "M":
            ax.set_xscale('log')
            plt.xlabel(r'M$_{\mathrm{tot}} \; (M_{\odot})$')
        elif Vars[0] == "inc":
            ax.set_xlim([0.,np.pi])
            plt.xlabel(r'$\iota \; (\mathrm{rad.})$')
        elif Vars[0] == "maxf":
            ax.set_xscale('log')
            plt.xlabel(r'$\mathrm{f}_{\mathrm{max}} \; (\mathrm{Hz})$')
        elif Vars[0] == "beta":
            plt.xlabel(r'$\beta_{SSB} \; (\mathrm{rad.})$')
        elif Vars[0] == "z":
            plt.xlabel(r'z')
        elif Vars[0] == "T_obs_end_to_merger":
            ax.set_xscale('log')
            plt.xlabel(r'T$_{obstm} \; (\mathrm{hr})$')

        if Vars[1] == "M":
            ax.set_yscale('log')
            ax.set_ylabel(r'M$_{\mathrm{tot}} \; (M_{\odot})$')
        elif Vars[1] == "inc":
            ax.set_ylabel(r'$\iota \; (\mathrm{rad.})$')
        elif Vars[1] == "maxf":
            ax.set_yscale('log')
            ax.set_ylabel(r'$\mathrm{f}_{\mathrm{max}} \; (\mathrm{Hz})$')
        elif Vars[1] == "beta":
            ax.set_ylabel(r'$\beta_{SSB} \; (\mathrm{rad.})$')
        elif Vars[1] == "z":
            ax.set_ylabel(r'z')
        elif Vars[1] == "T_obs_end_to_merger":
            ax.set_yscale('log')
            plt.ylabel(r'T$_{obstm} \; (\mathrm{hr})$')

        if FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA":
            cbar.set_label(r'$\log_{10}(|\mathcal{B}_{\beta}|)$') # (r'$\log_{10}(|\Sigma_{b=0}^{3}\log_e(L_{(-1,b)}) - \log_e(L_{(1,b)})+20|)$')
        elif FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
            cbar.set_label(r'$\log_{10}(|\mathcal{B}_{\lambda}|)$') # r'$\log_{10}(|\Sigma_{a=1,-1}\log_e(L_{(a,1)}) + \log_e(L_{(a,2)}) + \log_e(L_{(a,3)}) - \log_e(L_{(a,0)})+20|)$')
        else:
            cbar.set_label(r'$\log_{10}(\Delta \Omega \; (\mathrm{sq. deg.}))$')

    # Save?
    if SaveFig:
        if FoMOverRange["FoM"] == "SkyModeLikelihoodsLAMBDA":
            plt.savefig("/Users/jonathonbaird/Documents/LabEx_PostDoc/Figs/lnL_skymodes_LambdaReflections_" + Vars[0] + "_" + Vars[1] + "_BF_" + str(BF_lim) + ".png", facecolor='w', transparent=False)
        elif FoMOverRange["FoM"] == "SkyModeLikelihoodsBETA":
            plt.savefig("/Users/jonathonbaird/Documents/LabEx_PostDoc/Figs/lnL_skymodes_BetaReflections_" + Vars[0] + "_" + Vars[1] + "_BF_" + str(BF_lim) + ".png", facecolor='w', transparent=False)
        plt.show()
    else:
        return fig, ax

def PlotOrbit(config_struct,MollProj=False,SaveFig=False):
    """
    Plot orbit. Mostly to check we have done things right for orbital
    calculations in segments_athena.py.

    SkyProj=True will use "mollweide" projection but this takes radians as input
    and will give give some glitchy things for patch radii due to their small
    angular extents. If you set to cartesian (mollweide=False) patches will have
    radii in degrees and not have the same problem.
    """

    orbit_dict=config_struct["orbit_dict"]

    from mpl_toolkits.mplot3d import Axes3D
    import astropy.constants as consts

    Moon_RaDecs=orbit_dict["Moon_From_Athena_radecs"]
    Earth_RaDecs=orbit_dict["Earth_From_Athena_radecs"]
    Sun_RaDecs=orbit_dict["Sun_From_Athena_radecs"]

    if MollProj:
        # Then np.arctan will give radians and projections will work
        Factor=np.pi/180.
    else:
        # Then we convert np.arctan to degrees for rectangular plot
        Factor=1.

    # Everything in degrees so 'Factor' handles conversion if we chose a sky projection
    moon_radii=np.arctan(1737400./np.linalg.norm(orbit_dict["Moon_From_Athena"],axis=0))*180./np.pi
    earth_radii=np.arctan(consts.R_earth.value/np.linalg.norm(orbit_dict["Earth_From_Athena"],axis=0))*180./np.pi
    sun_radii=np.arctan(consts.R_sun.value/np.linalg.norm(orbit_dict["Sun_From_Athena"],axis=0))*180./np.pi

    fig = plt.figure()
    if MollProj:
        ax = plt.subplot(111, projection="mollweide")
    else:
        ax = plt.gca()
        ax.set_aspect(1)

    for i in range(len(Moon_RaDecs[0,:])):
        if not SkyProj and i==0:
            m = plt.Circle((Moon_RaDecs[0,i]*Factor, Moon_RaDecs[1,i]*Factor), moon_radii[i]*Factor, color="cyan",label="Moon")
            e = plt.Circle((Earth_RaDecs[0,i]*Factor, Earth_RaDecs[1,i]*Factor), earth_radii[i]*Factor, color="blue",label="Earth")
            s = plt.Circle((Sun_RaDecs[0,i]*Factor, Sun_RaDecs[1,i]*Factor), sun_radii[i]*Factor, color="red",label="Sun")
        else:
            m = plt.Circle((Moon_RaDecs[0,i]*Factor, Moon_RaDecs[1,i]*Factor), moon_radii[i]*Factor, color="cyan")
            e = plt.Circle((Earth_RaDecs[0,i]*Factor, Earth_RaDecs[1,i]*Factor), earth_radii[i]*Factor, color="blue")
            s = plt.Circle((Sun_RaDecs[0,i]*Factor, Sun_RaDecs[1,i]*Factor), sun_radii[i]*Factor, color="red")
        ax.add_patch(m)
        ax.add_patch(e)
        ax.add_patch(s)
    if not MollProj:
        plt.xlim([-180.,180.])
        plt.ylim([-90.,90.])
        plt.legend(fontsize="x-small")

    # Save?
    if SaveFig:
        t0 = config_struct["gps_science_start"]
        t = Time(t0, format='gps', scale='utc').isot
        f=SYNEX_PATH+"/Plots/orbits/"
        pathlib.Path(f).mkdir(parents=True, exist_ok=True)
        if MollProj:
            strname="Athena_mollweide_" + "".join(t.split("T")[0].split("-")) + "_" + str(int((config_struct["mission_duration"]*364.25)//1)) + "d_inc"+str(int(config_struct["inc"]//1))+"_R"+str(int(config_struct["MeanRadius"]//1e6))+"Mkm_ecc"+str(int(config_struct["eccentricity"]//0.1))
        else:
            strname="Athena_cartesian_" + "".join(t.split("T")[0].split("-")) + "_" + str(int((config_struct["mission_duration"]*364.25)//1)) + "d_inc"+str(int(config_struct["inc"]//1))+"_R"+str(int(config_struct["MeanRadius"]//1e6))+"Mkm_ecc"+str(int(config_struct["eccentricity"]//0.1))
        strname+="_ArgPeri"+str(int(config_struct["ArgPeriapsis"]//1))+"_AscNode"+str(int(config_struct["AscendingNode"]//1))+"_phi0"+str(int(config_struct["ArgPeriapsis"]//1))
        strname+="_P"+str(int(config_struct["period"]//1))+"_frozen"+str(config_struct["frozenAthena"])+".pdf"
        plt.savefig(f+strname,dpi=200)

    plt.show()

def AnimateOrbit(config_struct,include_sun=False,SaveAnim=False):
    """
    Animation of orbit. Mostly just to be cool.

    To change from 'mp4' to 'swf' format (eg for inclusion in latex slides), use command line:
    "ffmpeg -i in.mp4 out.swf"
    """
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import animation
    import astropy.constants as consts

    orbit_dict=config_struct["orbit_dict"]

    obj_keys = ["Moon_From_Athena","Earth_From_Athena"]
    if include_sun:
        obj_keys.append("Sun_From_Athena")
        Scale = consts.au.value # config_struct["MeanRadius"] #
        lim = 1.1 # 1.+consts.au.value/Scale #
    else:
        Scale = config_struct["MeanRadius"]
        L2_from_Earth = config_struct["L2_from_Earth"]
        lim = 1.+L2_from_Earth/Scale

    labels={"Moon_From_Athena":"Moon","Earth_From_Athena":"Earth","Sun_From_Athena":"Sun"}
    colors={"Moon_From_Athena":"cyan","Earth_From_Athena":"blue","Sun_From_Athena":"red"}
    sizes ={"Moon_From_Athena":(4.*np.pi/3.)*(1737400./Scale)**2,"Earth_From_Athena":(4.*np.pi/3.)*(consts.R_earth.value/Scale)**2,"Sun_From_Athena":(4.*np.pi/3.)*(consts.R_sun.value/Scale)**2}

    # coordinate data into scaled vectors to plot
    data = np.array([orbit_dict[key] for key in obj_keys])/Scale
    data_labels = np.array([labels[key] for key in obj_keys])
    data_colors = np.array([colors[key] for key in obj_keys])
    data_sizes = np.array([5.]*len(obj_keys)) # sizes[key] for key in obj_keys])
    N=np.shape(data)[2]
    data=np.append(data,np.zeros((1,3,N)),axis=0)
    data_labels=np.append(data_labels,"Athena")
    data_colors=np.append(data_colors,"black")
    data_sizes=np.append(data_sizes,5.)

    def update(num, data, sc, data_sizes):
        sc.set_offsets(data[:,0:2,num])
        sc.set_3d_properties(data[:,2,num],zdir='z')
        # sc.set_offset_position("data") # Can we use this for sizes too? -- NB This threw an error that the class doesnt have an attribute 'set_offset_position'
        sc.set_sizes(data_sizes)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    sc = ax.scatter(data[:,0,0],data[:,1,0],data[:,2,0],c=data_colors,s=data_sizes) # ,s=sizes,label=data_labels,marker="o")
    ax.set_xlim3d([-lim,lim])
    ax.set_ylim3d([-lim,lim])
    ax.set_zlim3d([-1.1,1.1])
    # ax.set_xlabel('Earth Orb Tangent',fontsize="x-small")
    # ax.set_ylabel('Sun-Earth Axis',fontsize="x-small")
    # ax.set_zlabel('Normal to Ecliptic plane',fontsize="x-small")
    ani = animation.FuncAnimation(fig, update, N, fargs=(data, sc, data_sizes), interval=1, blit=False)
    # Save?
    if SaveAnim:
        t0 = config_struct["gps_science_start"]
        t = Time(t0, format='gps', scale='utc').isot
        f=SYNEX_PATH+"/Plots/OrbitAnimations/"
        pathlib.Path(f).mkdir(parents=True, exist_ok=True)
        strname="Athena_" + "".join(t.split("T")[0].split("-")) + "_" + str(int((config_struct["mission_duration"]*364.25)//1)) + "d_inc"+str(int(config_struct["inc"]//1))+"_R"+str(int(config_struct["MeanRadius"]//1e6))+"Mkm_ecc"+str(int(config_struct["eccentricity"]//0.1))
        strname+="_ArgPeri"+str(int(config_struct["ArgPeriapsis"]//1))+"_AscNode"+str(int(config_struct["AscendingNode"]//1))+"_phi0"+str(int(config_struct["ArgPeriapsis"]//1))
        strname+="_P"+str(int(config_struct["period"]//1))+"_frozen"+str(config_struct["frozenAthena"])+".mp4"
        ani.save(f+strname, writer=animation.FFMpegWriter(fps=60)) # writer=animation.PillowWriter(fps=1)) # writer='imagemagick')
    plt.show()

def AnimateSkyProjOrbit(config_struct,MollProj=False,SaveAnim=False):
    """
    Animation of orbit using skymap projection. Mostly just to be cool.

    SkyProj=True will use "mollweide" projection but this takes radians as input
    and will give give some glitchy things for patch radii due to their small
    angular extents. If you set to cartesian (mollweide=False) patches will have
    radii in degrees and not have the same problem.

    To change from 'mp4' to 'swf' format (eg for inclusion in latex slides), use command line:
    "ffmpeg -i in.mp4 out.swf"
    """

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import animation
    import astropy.constants as consts

    if MollProj:
        # Then np.arctan will give radians and projections will work
        Factor=np.pi/180.
    else:
        # Then we convert np.arctan to degrees for rectangular plot
        Factor=1.

    orbit_dict=config_struct["orbit_dict"]

    Moon_RaDecs=Factor*orbit_dict["Moon_From_Athena_radecs"]
    Earth_RaDecs=Factor*orbit_dict["Earth_From_Athena_radecs"]
    Sun_RaDecs=Factor*orbit_dict["Sun_From_Athena_radecs"]

    # Everything in degrees so 'Factor' handles conversion if we chose a sky projection
    moon_radii=Factor*np.arctan(1737400./np.linalg.norm(orbit_dict["Moon_From_Athena"],axis=0))*180./np.pi
    earth_radii=Factor*np.arctan(consts.R_earth.value/np.linalg.norm(orbit_dict["Earth_From_Athena"],axis=0))*180./np.pi
    sun_radii=Factor*np.arctan(consts.R_sun.value/np.linalg.norm(orbit_dict["Sun_From_Athena"],axis=0))*180./np.pi

    fig = plt.figure()
    if MollProj:
        ax = plt.subplot(111, projection="mollweide")
    else:
        ax = plt.gca()
        ax.set_aspect(1)

    if not MollProj:
        m = plt.Circle((Moon_RaDecs[0,0], Moon_RaDecs[1,0]), moon_radii[0], color="cyan",label="Moon")
        e = plt.Circle((Earth_RaDecs[0,0], Earth_RaDecs[1,0]), earth_radii[0], color="blue",label="Earth")
        s = plt.Circle((Sun_RaDecs[0,0], Sun_RaDecs[1,0]), sun_radii[0], color="red",label="Sun")
    else:
        m = plt.Circle((Moon_RaDecs[0,0], Moon_RaDecs[1,0]), moon_radii[0], color="cyan")
        e = plt.Circle((Earth_RaDecs[0,0], Earth_RaDecs[1,0]), earth_radii[0], color="blue")
        s = plt.Circle((Sun_RaDecs[0,0], Sun_RaDecs[1,0]), sun_radii[0], color="red")
    ax.add_patch(m)
    ax.add_patch(e)
    ax.add_patch(s)

    if not MollProj:
        plt.xlim([-180.,180.])
        plt.ylim([-90.,90.])
        plt.legend(fontsize="x-small")

    def update(num, ax, Moon_RaDecs, Earth_RaDecs, Sun_RaDecs, moon_radii, earth_radii, sun_radii):
        ax.patches[0].set(center=(Moon_RaDecs[0,num], Moon_RaDecs[1,num]), radius=moon_radii[num])
        ax.patches[1].set(center=(Earth_RaDecs[0,num], Earth_RaDecs[1,num]), radius=earth_radii[num])
        ax.patches[2].set(center=(Sun_RaDecs[0,num], Sun_RaDecs[1,num]), radius=sun_radii[num])

    ax.set_xlabel('Athena RA',fontsize="x-small")
    ax.set_ylabel('Athena Dec',fontsize="x-small")
    N = len(Moon_RaDecs[0,:])
    ani = animation.FuncAnimation(fig, update, N, fargs=(ax, Moon_RaDecs, Earth_RaDecs, Sun_RaDecs, moon_radii, earth_radii, sun_radii), interval=2, blit=False)
    # Save?
    if SaveAnim:
        t0 = config_struct["gps_science_start"]
        t = Time(t0, format='gps', scale='utc').isot
        f=SYNEX_PATH+"/Plots/OrbitAnimations/"
        pathlib.Path(f).mkdir(parents=True, exist_ok=True)
        if MollProj:
            strname="Athena_SkyProjOrbit_mollweide_" + "".join(t.split("T")[0].split("-")) + "_" + str(int((config_struct["mission_duration"]*364.25)//1)) + "d_inc"+str(int(config_struct["inc"]//1))+"_R"+str(int(config_struct["MeanRadius"]//1e6))+"Mkm_ecc"+str(int(config_struct["eccentricity"]//0.1))
        else:
            strname="Athena_SkyProjOrbit_cartesian_" + "".join(t.split("T")[0].split("-")) + "_" + str(int((config_struct["mission_duration"]*364.25)//1)) + "d_inc"+str(int(config_struct["inc"]//1))+"_R"+str(int(config_struct["MeanRadius"]//1e6))+"Mkm_ecc"+str(int(config_struct["eccentricity"]//0.1))
        strname+="_ArgPeri"+str(int(config_struct["ArgPeriapsis"]//1))+"_AscNode"+str(int(config_struct["AscendingNode"]//1))+"_phi0"+str(int(config_struct["ArgPeriapsis"]//1))
        strname+="_P"+str(int(config_struct["period"]//1))+"_frozen"+str(config_struct["frozenAthena"])+".mp4"
        ani.save(f+strname, writer=animation.FFMpegWriter(fps=60)) # writer=animation.PillowWriter(fps=1)) # writer='imagemagick')
    plt.show()

def PlotEMFlux(source, SaveFig=False):
    xray_time=source.EM_Flux_Data["xray_time"]
    xray_flux=source.EM_Flux_Data["xray_flux"]

    plt.plot([t/(24.*60.*60.) for t in xray_time],xray_flux)
    ax=plt.gca()
    ax.set_yscale('log')
    plt.ylabel(r"Flux [erg s$^{-1}$ cm$^{-2}$]", fontsize="xx-small")
    plt.xlabel(r"Time [d]", fontsize="xx-small")
    plt.grid()

    # Save?
    if SaveFig:
        strname=SYNEX_PATH+"/Plots/"+".".join(source.ExistentialFileName.split("/Saved_Source_Dicts/")[-1].split(".")[:-1])+".pdf"
        f="/".join(strname.split("/")[:-1])
        pathlib.Path(f).mkdir(parents=True, exist_ok=True)
        plt.savefig(strname)

    plt.show()

def PlotCTR(source, SaveFig=False):
    xray_time=source.EM_Flux_Data["xray_time"]
    CTRs=source.CTR_Data["CTR"]

    plt.plot([t/(24.*60.*60.) for t in xray_time],CTRs)
    ax=plt.gca()
    ax.set_yscale('log')
    plt.ylabel(r"CTR [ ]", fontsize="xx-small")
    plt.xlabel(r"Time [d]", fontsize="xx-small")
    plt.grid()

    # Save?
    if SaveFig:
        strname=SYNEX_PATH+"/Plots/"+".".join(source.ExistentialFileName.split("/Saved_Source_Dicts/")[-1].split(".")[:-1])+".pdf"
        f="/".join(strname.split("/")[:-1])
        pathlib.Path(f).mkdir(parents=True, exist_ok=True)
        plt.savefig(strname)

    plt.show()

def PlotPhotonAccumulation(detectors, SaveFig=False, SaveFileName=None):
    """
    Plot accumulated photons from source only after tiling is run for a list of a list, a list or a single detector.
    """
    if not isinstance(detectors,list):
        detectors=[detectors]

    # Check if it's a list of lists
    if any([isinstance(detector,list) for detector in detectors]):
        detectors=[[detectors_el] if not isinstance(detectors[0],list) else detectors_el for detectors_el in detectors]
    else:
        detectors=[detectors]

    # Create list of detectors if dicts given or mix of dicts and detector objects given
    # detectors=[detector if not isinstance(detector,dict) else SYDs.Athena(**detector) for detector in detectors]

    # Create list of sources in case we have a mix of stuff like DeltatL_cut
    sources=[[GetSourceFromLisabetaData(detector.detector_source_coverage["source H5File"],**{"ExistentialFileName":detector.detector_source_coverage["source save file"]}) for detector in detectors_el] for detectors_el in detectors]
    SourcePreMergerCuts=[[source.DeltatL_cut for source in source_el] for source_el in sources]
    for ii in range(len(detectors)):
        if len(np.unique(SourcePreMergerCuts[ii]))>1:
            for jj in range(len(detectors[ii])):
                detectors[ii][jj].detector_config_struct.update({"cloning key":"DeltatL_cut","cloning value":-sources[ii][jj].DeltatL_cut/86400.})
    t_mergers=[[-source.DeltatL_cut/86400. for source in source_el] for source_el in sources]

    # Find the max height so all bars are even
    for ii in range(len(detectors)):
        height=max([sum(detector.detector_source_coverage["Source photon counts"]) for detector in detectors[ii]])

        for detector,t_merger in zip(detectors[ii],t_mergers[ii]):
            T0_mjd=detector.detector_source_coverage["Start time (mjd)"]
            Xs=[ExT0s/86400. for ExT0s in detector.detector_source_coverage["Source tile start times (s)"]]
            Xs_Widths=[ExTs/86400. for ExTs in detector.detector_source_coverage["Source tile exposuretimes (s)"]]
            Ys=list(np.cumsum(detector.detector_source_coverage["Source photon counts"]))
            
            # Plot accumulated photons is source was captured
            if len(detector.detector_source_coverage["Source tile start times (s)"])>0:
                if detector.detector_config_struct["cloning value"]!="DeltatL_cut" and detector.detector_config_struct["cloning value"]!=None:
                    label=detector.detector_config_struct["cloning key"]+r"="+str(int(detector.detector_config_struct["cloning value"]//1)) if isinstance(detector.detector_config_struct["cloning value"],(float,int)) else detector.detector_config_struct["cloning value"]
                elif detector.detector_config_struct["cloning value"]=="DeltatL_cut":
                    label=r"System " + str(ii) + r", Source time to merger="+str(int(t_merger//1))+r"\,d"
                else:
                    label=detector.detector_config_struct["telescope"]
                step_plot=plt.step([t-t_merger+t_exp for t,t_exp in zip(Xs,Xs_Widths)], Ys, where='post', label=label)
                plt.bar([t-t_merger for t in Xs],height=height,width=Xs_Widths,align="edge",color=step_plot[-1].get_color(),alpha=0.3)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.xlabel(r"Time to Merger [d]",fontsize="xx-small")
    plt.ylabel(r"Accumulated Photons",fontsize="xx-small")
    plt.xlim([-max([max(ts) for ts in t_mergers]),0.])
    plt.legend(fontsize="xx-small")

    # Save?
    if SaveFig:
        if SaveFileName==None:
            strname=SYNEX_PATH+"/Plots/"+".".join(source.ExistentialFileName.split("/Saved_Source_Dicts/")[-1].split(".")[:-1])+".pdf"
        else:
            SaveFileName=SYNEX_PATH+"/"+SaveFileName.split("/SYNEX/")[-1]
            SaveFileName=".".join(SaveFileName.split(".")[:-1])+".pdf"
        f="/".join(strname.split("/")[:-1])
        pathlib.Path(f).mkdir(parents=True, exist_ok=True)
        plt.savefig(strname)

    plt.show()

def PlotSourcePhotons(detectors, labels=None, SaveFig=False, SaveFileName=None):
    """
    Plot total accumulated photons for each list of detectors (in case we
    have several versions of the same cloned parameter but using eg different
    gwemopt flags...)
    """
    # Check if it's a list
    if not isinstance(detectors,list):
        detectors=[detectors]

    # Check if it's a list of lists
    if any([isinstance(detector,list) for detector in detectors]):
        detectors=[[detectors_el] if not isinstance(detectors[0],list) else detectors_el for detectors_el in detectors]
    else:
        detectors=[detectors]

    # Handle cases where some/all entries are dictionaries instead of detector classes
    # This does not preserve the shape of input detectors list...
    # detectors=[detector if not isinstance(detector,dict) else SYDs.Athena(**detector) for detectors_el in detectors for detector in detectors]

    # Create list of sources in case we have a mix of stuff like DeltatL_cut
    sources=[[GetSourceFromLisabetaData(detector.detector_source_coverage["source H5File"],**{"ExistentialFileName":detector.detector_source_coverage["source save file"]}) for detector in detectors_el] for detectors_el in detectors]
    SourcePreMergerCuts=[[source.DeltatL_cut for source in source_el] for source_el in sources]
    for ii in range(len(detectors)):
        if len(np.unique(SourcePreMergerCuts[ii]))>1:
            for jj in range(len(detectors[ii])):
                detectors[ii][jj].detector_config_struct.update({"cloning key":"DeltatL_cut","cloning value":-sources[ii][jj].DeltatL_cut/86400.})

    # Handle labels now if not given
    if labels==None: labels=[None]*len(detectors)

    # Init global figure
    fig=plt.figure()

    # Get the figure and axes for each list within detectors list
    Lines2D=[(PlotSourcePhotons_SingleDetList(detectors_el,sources=sources,fig=fig,label=label)) for detectors_el,label in zip(detectors,labels)]

    # Formatting some stuff globally
    if labels!=None:
        plt.legend(fontsize="x-small")

    # Save?
    if SaveFig:
        if SaveFileName==None:
            strname=SYNEX_PATH+"/Plots/"+".".join(source.ExistentialFileName.split("/Saved_Source_Dicts/")[-1].split(".")[:-1])+".pdf"
        else:
            SaveFileName=SYNEX_PATH+"/"+SaveFileName.split("/SYNEX/")[-1]
            SaveFileName=".".join(SaveFileName.split(".")[:-1])+".pdf"
        f="/".join(strname.split("/")[:-1])
        pathlib.Path(f).mkdir(parents=True, exist_ok=True)
        plt.savefig(strname)

    plt.show()

def PlotSourcePhotons_SingleDetList(detectors,sources=None,fig=None,label=None):
    """
    Plot total accumulated photons from source versus cloned parameter.
    """
    if not isinstance(detectors,list):
        detectors=[detectors]

    # Create list of detectors if dicts given or mix of dicts and detector objects given
    detectors=[detector if not isinstance(detector,dict) else SYDs.Athena(**detector) for detector in detectors]

    # Create sources if we don't have any
    if sources==None: sources=[GetSourceFromLisabetaData(detector.detector_source_coverage["source H5File"],**{"ExistentialFileName":detector.detector_source_coverage["source save file"]}) for detector in detectors]

    # Get cloned values and photon counts
    Xs = [detector.detector_config_struct["cloning value"] for detector in detectors] # This will store values cloned for sources too...
    Ys = [sum(detector.detector_source_coverage["Source photon counts"]) for detector in detectors]

    # Do we need logarithmic axes?
    ax = plt.gca()
    if min(Xs)!=0. and np.log10(max(Xs)/min(Xs))>1.5:
        ax.set_xscale('log')
    if min(Ys)!=0. and np.log10(max(Ys)/min(Ys))>2.5:
        ax.set_yscale('log')

    # Label axes
    plt.xlabel(detectors[0].detector_config_struct["cloning key"],fontsize="x-small")
    plt.ylabel(r"Accumulated Photons",fontsize="x-small")

    # Init figure if it isn't already
    if fig==None: fig=plt.figure()
    if label==None: label=''

    # Plot
    ErrBarCont=plt.errorbar(Xs,Ys,yerr=np.sqrt(Ys), marker='+', linestyle='', markersize=4, figure=fig, label=label)

    # Return
    return ErrBarCont, ax






#              ###############################################
#              #                                             #
#              #   SYNEX -- Optimization Utility functions   #
#              #                                             #
#              ###############################################

### This can be used if we have a saved skymap and Athena tiling strategy, then optimize the given schedule using latency time / jump time / other params...
def RunPSO(source, detector, fitness_function=None, N=50, w=0.8, c_1=1, c_2=1, auto_coef=True, max_iter=100, NumSwarms=1, priors=None, **kwargs):
    print("Initializing Swarm(s)...")
    PSO_Classes = [PSO(fitness_function, priors=priors,N=N, w=w, c_1=c_1, c_2=c_2,
                       max_iter=max_iter, auto_coef=True) for iSwarm in range(NumSwarms)]
    PSO_IsRunning = [PSO_Class.is_running for PSO_Class in PSO_Classes]

    # Pre-run output
    print("Running PSO... ")
    ProgressBarText = "PSO progress"
    StartProgress(ProgressBarText)

    # Run
    while all(PSO_IsRunning):
        # Step each swarm forward 1 generation at a time
        for PSO_Class in PSO_Classes:
            PSO_Class.next()
        if PSO_Classes[0].iter % 10 == 0:
            # Plot snapshot
            PlotPSOSnapshot(PSO_Classes)
            # Output progress bar
            progress = 100.*PSO_Classes[0].iter/PSO_Classes[0].max_iter
            UpdateProgress(progress)

        # Update condition list for outer while loop
        PSO_IsRunning = [PSO_Class.is_running for PSO_Class in PSO_Classes]

    # Close progress bar
    EndProgress()

    # print the final variable values
    iSwarm = 1
    for PSO_Class in PSO_Classes:
        OutStr = "Swarm " + str(iSwarm) + "/" + str(NumSwarms) + " optimal loc: " + str(PSO_Class.g_best)
        OutStr += ", " + PSO_Class.__str__()
        print(OutStr)
        iSwarm += 1

def PlotPSOSnapshot(PSO_Classes):
    surf_xs = np.linspace(-5.,5.,20)
    surf_ys = np.linspace(-5.,5.,20)
    surf_X, surf_Y = np.meshgrid(surf_xs, surf_ys)
    particles_surf = [[x,y] for x,y in zip(np.ndarray.flatten(surf_X), np.ndarray.flatten(surf_Y))]
    surf_Z = PSO_Classes[0].fitness_function(particles_surf)
    N_surf = len(particles_surf)
    surf_Z = np.reshape(surf_Z,[int(np.sqrt(N_surf)),int(np.sqrt(N_surf))])

    # Plot the surface.
    fig, ax = plt.subplots() # subplot_kw={"projection": "3d"})
    plt.contour(surf_X, surf_Y, surf_Z, 20, cmap='RdGy');
    # surf = ax.plot_surface(surf_X, surf_Y, surf_Z, linewidth=0, antialiased=False, cmap=mpl.cm.coolwarm,)

    # Extract particle positions for each swarm
    for PSO_Class in PSO_Classes:
        parts = PSO_Class.particles
        Vels = PSO_Class.velocities
        FuncVals = PSO_Class.fitness_function(PSO_Class.particles)

        # Plot particles
        ax.scatter(parts[:,0], parts[:,1], s=2.5) # s=FuncVals,
        # ax.scatter(parts[:,0], parts[:,1], FuncVals, s=2.5) # c="black",

    # Show
    plt.show()

class PSO:
    def __init__(self, fitness_function, particles=None, velocities=None, priors=None,
                 N=100,w=0.8, c_1=1, c_2=1, max_iter=100, auto_coef=True):

        # Set the class variables
        self.N = N # Number of particles in swarm
        self.ndim = len(priors["infer_params"]) # Number of params defining search volume
        self.w = w # Particle inertia
        self.c_1 = c_1 # Cognitive coefficient
        self.c_2 = c_2 # Social coefficient
        self.auto_coef = auto_coef # Adapt coefficient at each generation
        self.max_iter = max_iter # Max number of generations (iterations of search)
        self.priors = priors # Dictionary of prior ranges and priors matching prior dict for lisabeta
        if isinstance(self.priors, dict):
            self.params_range = np.array(self.priors["params_range"]) # ndim by 2 np.array
            self.SearchWidths = np.array([u-l for l,u in self.params_range]) # ndim np.array

        # Initiate particles and velocities if not already done
        if all([particles==None, velocities==None, priors==None]):
            raise ValueError("To initiate PSO you need to specify either particles and velocities, or priors")
        elif all([particles==None, velocities==None, isinstance(self.priors, dict)]):
            # Stochastic particle initial positions and velocities (within prior bounds for each dimension)
            particles = []
            velocities = []
            for idim in range(self.ndim):
                partiles_dim = self.SearchWidths[idim]*(np.random.random(self.N)-0.5)*0.9 # factor 0.9 to initialize far enough from the boundary to stop stray particles
                particles.append(partiles_dim)
                vel_xs = self.SearchWidths[idim]*0.000001*(np.random.random(self.N)-0.5) # Initial velocity doesnt seem to change much the test runs
                velocities.append(vel_xs)
            particles = np.reshape(np.array(particles),(self.N,self.ndim))
            velocities = np.reshape(np.array(velocities),(self.N,self.ndim))

        # Store as class variables
        self.particles = particles
        self.velocities = velocities
        self.fitness_function = fitness_function

        # Initiate local and global best estimates
        self.p_bests = self.particles
        self.p_bests_values = self.fitness_function(self.particles)
        self.g_best = self.p_bests[0]
        self.g_best_value = self.p_bests_values[0]
        self.update_bests()

        # Set some runtime variables
        self.iter = 0
        self.is_running = True
        self.update_coef()

    def __str__(self):
        return f'[{self.iter}/{self.max_iter}] $w$:{self.w:.3f} - $c_1$:{self.c_1:.3f} - $c_2$:{self.c_2:.3f}'

    def next(self):
        if self.iter > 0:
            self.move_particles()
            self.update_bests()
            self.update_coef()

        self.iter += 1
        self.is_running = self.is_running and self.iter < self.max_iter

        return self.is_running

    def update_coef(self):
        if self.auto_coef:
            t = self.iter
            n = self.max_iter
            self.w = (0.4/n**2) * (t - n) ** 2 + 0.4
            self.c_1 = -3 * t / n + 3.5
            self.c_2 =  3 * t / n + 0.5

    def move_particles(self):

        # add inertia
        new_velocities = self.w * self.velocities

        # add cognitive component
        r_1 = np.random.random(self.N)
        r_1 = np.tile(r_1[:, None], (1, self.ndim))

        new_velocities += self.c_1 * r_1 * (self.p_bests - self.particles)/self.SearchWidths

        # add social component
        r_2 = np.random.random(self.N)
        r_2 = np.tile(r_2[:, None], (1, self.ndim))
        g_best = np.tile(self.g_best[None], (self.N, 1))
        new_velocities += self.c_2 * r_2 * (g_best  - self.particles)/self.SearchWidths

        self.is_running = np.sum(self.velocities - new_velocities) != 0

        # Scale velocity components so it does not exceed some upper limit,
        # and  heck particles do not exit the allowed explorable volume (prior ranges)
        if isinstance(self.priors, dict)  and isinstance(self.priors, dict):
            for ivel in range(self.N):
                for idim in range(self.ndim):
                    # Check if velocity is too large - I think this might be obsoolete after adding factor 1/SearchWidths tp social and cognitive velocity components above...
                    if new_velocities[ivel,idim]/self.SearchWidths[idim] >= 0.2:
                        new_velocities[ivel,idim] = 0.2*self.SearchWidths[idim]
                    # Check if velocity will put particles out of region, and reverse vel component if true
                    # Velocity should always be small enough that if the above conditions are true, reversing
                    # the vel component should not place the particle out of the other side of the search volume
                    if self.particles[ivel,idim] + new_velocities[ivel,idim] >= self.params_range[idim,1]:
                        new_velocities[ivel,idim] = -new_velocities[ivel,idim]
                    if self.particles[ivel,idim] + new_velocities[ivel,idim] <= self.params_range[idim,0]:
                        new_velocities[ivel,idim] = -new_velocities[ivel,idim]

        # Update particle and particle velocities
        self.velocities = new_velocities
        self.particles = self.particles + new_velocities

    def update_bests(self):
        fits = self.fitness_function(self.particles)

        for i in range(len(self.particles)):
            # update best personnal value (cognitive)
            if fits[i] < self.p_bests_values[i]:
                self.p_bests_values[i] = fits[i]
                self.p_bests[i] = self.particles[i]
                # update best global value (social)
                if fits[i] < self.g_best_value:
                    self.g_best_value = fits[i]
                    self.g_best = self.particles[i]

# Fancy progress bar functions
def StartProgress(title):
    global progress_old
    sys.stdout.write(title + ": [" + "-"*40 + "]" + chr(8)*41)
    sys.stdout.flush()
    progress_old = 0
def UpdateProgress(x):
    global progress_old
    x = int(x * 40 // 100)
    sys.stdout.write("#" * (x - progress_old))
    sys.stdout.flush()
    progress_old = x
def EndProgress():
    sys.stdout.write("#" * (40 - progress_old) + "]\n")
    sys.stdout.flush()
