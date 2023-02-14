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
import glob
# Set the path to root SYNEX directory so file reading and writing is easier -- NOTE: NO TRAILING SLASH
SYNEX_PATH=str(pathlib.Path(__file__).parent.resolve()).split("SYNEX")[0]+"SYNEX"
import astropy.constants as const
import astropy.units as u
from astropy.time import Time
import numpy as np
import healpy as hp
import pandas as pd
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
import itertools
from scipy import stats
import time
import json, h5py

# Import lisebeta stuff
import lisabeta
import lisabeta.lisa.lisa_fisher as lisa_fisher
import lisabeta.lisa.snrtools as lisa_snr
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
import SYNEX.SYNEX_Telescopes as SYTs
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

# global bins variable for histogram plots and mode calculations.
hist_n_bins=1000

# List of main source parameters according to lisabeta input dictionaries.
list_params = [
    "M",
    "q",
    "chi1",
    "chi2",
    "Deltat",
    "dist",
    "inc",
    "phi",
    "lambda",
    "beta",
    "psi"]






#                     #################################
#                     #                               #
#                     #   General Utility functions   #
#                     #                               #
#                     #################################

def ClassesToParams(source, detector=None, CollectionMethod="Inference",**kwargs):
    """
    Function to change parameters embedded in the classes and handed to functions
    in dictionaries, to dictionaries and parameter lists for function calls in lisabeta.

    This function is clunky and could probably use some cleaning up...

    PARAMS
    ------
        - Source : SYNEX source object.
        - Detector : SYNEX GW detector object

    NOTE:
    -----
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
              "Lframe": source.Lframe,        # Params always in Lframe - this is easier for the sampler.
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
    LNKs=waveform_params["LISAnoise"].keys()
    for key,val in kwargs.items():
        if key in waveform_params: waveform_params[key]=val
        if key in LNKs: waveform_params["LISAnoise"][key]=val

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
    """
    Function to return SYNEX classes from dictionary of parameters.

    PARAMS
    ------
        - input_params : Dictionary
            Params to give to lisabeta inference methods. Usually saved by SYNEX
            in Json format in file SYNEX/inference_param_files/...

        - CollectionMethod : String
            Method of formatting dictionaries depending on what 'input_params' was
            used for. Options are 'Inference' or 'Fisher' depending on if the
            output from the lisabeta launch was a Fisher localization estimate
            or a full MCMC inference (with ptemcee inside lisabeta).

        - kwargs : Dictionary
            Optional dictionary of values to adjust the input params. This can be
            useful if we want many copies of the same objects with a subset of
            parameters changed. In this case we would call this function many times
            with the same input parameters but with **kwargs containing the adjusted
            param value at each call iteration.


    OUTPUT
    ------
        - Source : SYNEX source class
            Source class object created wth attributed given in input_params and
            modified by anything in kwargs, using formatting according to 'CollectionMethod'.

    TO DO:
    ------
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
    source_kwargs = {**param_dict,
                     "lisabetaFile":run_params["out_dir"] +
                     run_params["out_name"] + ".h5", **waveform_params}

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

# def CompleteLisabetaDataAndJsonFileNames(FileName,
#                                          rep_param="inference_param_files",
#                                          rep_inf="inference_data"):
#     """
#     Function to give back congruent json and h5 filenames, complete with paths and
#     proper file extensions. This is particularly useful when transfering files
#     from cluster to local workspaces where locations of files may be different
#     (e.g. SYNEX_PATH is different, AND we wish to group cluster runs into sub directories
#     according to test generation etc).

#     PARAMS
#     ------
#         - FileName :: String
#             Filename, with or without '.json' or '.h5' extensions, including any
#             desired architecture existing below the locations './SYNEX/inference_data'
#             or './SYNEX/inference_param_files'.
#     """

#     # Create paths to data and json folders
#     LisabetaJsonPath = SYNEX_PATH + "/inference_param_files"
#     # LisabetaJsonPath = os.path.join(SYNEX_PATH, rep_param)
#     LisabetaDataPath = SYNEX_PATH + "/inference_data"
#     # LisabetaDataPath = os.path.join(SYNEX_PATH, rep_inf)

#     print(FileName)
#     if FileName[0]!="/":
#         LisabetaJsonPath+="/"
#         LisabetaDataPath+="/"

#     # Figure out if its a json or h5 filename
#     if FileName[-3]==".":
#         JsonFileLocAndName = FileName[:-3] + '.json'
#         H5FileLocAndName = FileName
#     elif FileName[-5]==".":
#         JsonFileLocAndName = FileName
#         H5FileLocAndName = FileName[:-5] + '.h5'
#     else:
#         JsonFileLocAndName = FileName + '.json'
#         H5FileLocAndName = FileName + '.h5'

#     # Add file path if only the filenames specified
#     bool_vec = [len(FileName.split("inference_param_files"))==1,
#                 len(FileName.split("inference_data"))==1]
#     if bool_vec[0] and bool_vec[1]:
#         H5FileLocAndName = LisabetaDataPath + H5FileLocAndName
#         JsonFileLocAndName = LisabetaJsonPath + JsonFileLocAndName
#     elif bool_vec[0] and not bool_vec[1]:
#         JsonFileLocAndName = LisabetaJsonPath + H5FileLocAndName.split("inference_data")[-1]
#         JsonFileLocAndName = ".".join(JsonFileLocAndName.split(".")[:-1])+".json"
#     elif not bool_vec[0] and bool_vec[1]:
#         H5FileLocAndName = LisabetaDataPath + JsonFileLocAndName.split("inference_param_files")[-1]
#         H5FileLocAndName = ".".join(H5FileLocAndName.split(".")[:-1])+".h5"

#     # Check if the subdirectories exist for both data and json files
#     JsonPathOnly="/".join(JsonFileLocAndName.split("/")[:-1])
#     DataPathOnly="/".join(H5FileLocAndName.split("/")[:-1])
#     try:
#         # See if the directories exist in case we load ource from savefile or something
#         pathlib.Path(JsonPathOnly).mkdir(parents=True, exist_ok=True)
#         pathlib.Path(DataPathOnly).mkdir(parents=True, exist_ok=True)
#     except:
#         JsonPathOnly=LisabetaJsonPath
#         DataPathOnly=LisabetaDataPath
#         JsonFileLocAndName=JsonPathOnly+JsonFileLocAndName.split("inference_param_files")[-1]
#         H5FileLocAndName=DataPathOnly+H5FileLocAndName.split("inference_data")[-1]
#         print(JsonPathOnly)
#         pathlib.Path(JsonPathOnly).mkdir(parents=True, exist_ok=True)
#         pathlib.Path(DataPathOnly).mkdir(parents=True, exist_ok=True)

#     # Return full paths and names all harmonious and what not
#     return JsonFileLocAndName,H5FileLocAndName

def CompleteLisabetaDataAndJsonFileNames(FileName, path=SYNEX_PATH):
    """
    Function to give back congruent json and h5 filenames with their complete
    paths. This is particularly useful when transfering files
    from cluster to local workspaces where locations of files may be different
    (e.g. SYNEX_PATH is different, AND we wish to group cluster runs into
    sub-directories according to test generation etc).

    PARAMS
    ------
        - FileName :: String
            Filename, with or without '.json' or '.h5' extensions
        - path :: String
            Path where the subdirectories with the json and h5 files are stored
            Default: SYNEX path
    """

    if ("inference_param_files" in FileName): # check if the dirs are already OK
        splited = FileName.split("inference_param_files")
        path = splited[0]
        FileName = splited[1].split('.')[0] # Remove extension if necessary

    elif ("inference_data" in FileName):
        splited = FileName.split("inference_data")
        path = splited[0]
        FileName = splited[1].split('.')[0].split('/')[1] # Remove extension if necessary
    else:
        # Find the filename for the H5 and json files
        if ("json" in FileName) or ("h5" in FileName):
            FileName = FileName.split('.')[0].split('/')[1] # Remove extension if necessary
    # Setup the repertory path containing the parameter and data files.``
    LisabetaJsonPath = os.path.join(path, "inference_param_files")
    LisabetaDataPath = os.path.join(path,  "inference_data")

    # Join path, filenames and extension
    JsonFileLocAndName = os.path.join(LisabetaJsonPath, FileName + '.json')
    H5FileLocAndName = os.path.join(LisabetaDataPath, FileName + '.h5')

    # Create all the relevant repertories if they do not exist yet
    if not os.path.exists(LisabetaJsonPath):
        os.makedirs(LisabetaJsonPath)
    if not os.path.exists(LisabetaDataPath):
        os.makedirs(LisabetaDataPath)

    return JsonFileLocAndName, H5FileLocAndName

def GWEMOPTPathChecks(go_params, config_struct):
    """
    Function to check formatting of file paths specified at initiation of gwemopt
    params dictionary. This will check the locations as well as file names are
    coherent for all file names.

    NOTE: We assume that all files will be within the allocated SYNEX filespaces
    to remove possabilities of problems at run time after transfering classes over
    from cluster.

    PARAMS
    ------
        - go_params :: Dict
            GWEMOpt 'params' dict, based on the telescope.telescope_go_params base dictionary,
            but passed through the 'PrepareGwemoptDicts()' function with a source
            to prepare for running through GWEMOPT.
        - config_struct :: Dict
            GWEMOPT 'config_struct' dictionary containing relevant information
            describing a telescope object. Based on the telescope.telescope_config_struct
            base dictionary, but passed through the 'PrepareGwemoptDicts()' function with
            a source to prepare dictionaries for running through GWEOPT.

    OUTPUT
    ------
        - go_params
            Same as input but with filename fields congruent.
        - config_struct
            Same as input but with filename fields congruent.
    """

    PATH_VARS = ["outputDir", "tilingDir", "catalogDir"] # , "configDirectory"]
    FILE_VARS = ["coverageFiles", "lightcurveFiles"] ### lightcurveFiles will be used when we eventually interface with EM stuff
    CONFIG_FILE_VARS = ["tesselationFile","orbitFile"]

    # Now specified files
    for FileVar in FILE_VARS+PATH_VARS:
        # Strip out some stuff that can change when transfering files between systems
        if "/SYNEX/" in go_params[FileVar]: go_params[FileVar]=go_params[FileVar].split("/SYNEX/")[-1]
        if len(go_params[FileVar].split("."))>1: go_params[FileVar] = ".".join(go_params[FileVar].split(".")[:-1])
        go_params[FileVar].strip("/") # No trailing slashes either side of path / file

        # Add file ends
        if FileVar in FILE_VARS: go_params[FileVar] += ".dat"

        # Add correct directories
        if FileVar=="coverageFiles":
            if "gwemopt_cover_files/" not in go_params[FileVar]: go_params[FileVar]="gwemopt_cover_files/"+go_params[FileVar]
        elif FileVar=="lightcurveFiles":
            if "gwemopt/lightcurves/" not in go_params[FileVar]: go_params[FileVar]="gwemopt/lightcurves/"+go_params[FileVar]
        elif FileVar=="outputDir":
            if "gwemopt_output" not in go_params[FileVar]: go_params[FileVar]="gwemopt_output/"+go_params[FileVar]
        elif FileVar=="tilingDir":
            if "Tile_files" not in go_params[FileVar]: go_params[FileVar]="Tile_files/"+go_params[FileVar]
        elif FileVar=="catalogDir":
            if "gwemopt_catalogs" not in go_params[FileVar]: go_params[FileVar]="gwemopt_catalogs/"+go_params[FileVar]

        # Add correct SYNEX_PATH for current system
        go_params[FileVar] = SYNEX_PATH + "/" + go_params[FileVar]

        # Check if paths exist yet
        PathOnly = go_params[FileVar]
        if FileVar in FILE_VARS:
            PathOnly = "/".join(go_params[FileVar].split("/")[:-1])
        pathlib.Path(PathOnly).mkdir(parents=True, exist_ok=True)

    # Now config_struct paths
    for FileVar in CONFIG_FILE_VARS:
        # Strip out some stuff that can change when transfering files between systems
        if "/SYNEX/" in config_struct[FileVar]: config_struct[FileVar]=config_struct[FileVar].split("/SYNEX/")[-1]
        if len(config_struct[FileVar].split("."))>1: config_struct[FileVar] = ".".join(config_struct[FileVar].split(".")[:-1])
        config_struct[FileVar].strip("/") # No trailing slashes either side of path / file

        # Add file ends
        if FileVar in ["tesselationFile"]:
            config_struct[FileVar] += ".tess"
        elif FileVar in ["orbitFile"]:
            config_struct[FileVar] += ".dat"

        # Add correct directories
        if FileVar=="tesselationFile":
            if "gwemopt_tess_files/" not in config_struct[FileVar]: config_struct[FileVar]="gwemopt_tess_files/"+config_struct[FileVar]
        elif FileVar=="orbitFile":
            if "orbit_files/" not in config_struct[FileVar]: config_struct[FileVar]="orbit_files/"+config_struct[FileVar]

        # Add correct SYNEX_PATH for current system
        config_struct[FileVar] = SYNEX_PATH + "/" + config_struct[FileVar]

        # Check if paths exist yet
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

    ################# This function is now redundant #################
    This function was used before 'post tiling analysis' scripts were made to
    extract information from classes into dataframes that are much faster and
    reliable for calculating statistics.
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

def GetTotSkyAreaFromPostData_OLD(FileName,ConfLevel=0.9):
    """
    Counting function to give the total sky area for a confidence level given eiher
    a Json or h5 filename for lisabeta data.

    A 2D histogram is created to bin the drawn lambda and beta values saved.
    We do not use the log posterior probability here and we probably should.

    REQUIRES lisabeta full inference to have been run to completion.


    ################# This function is now redundant #################
    """
    # Check what filename we were given
    JsonFileLocAndName,H5FileLocAndName=CompleteLisabetaDataAndJsonFileNames(FileName)

    # Unpack data
    [infer_params, inj_param_vals, static_params, meta_data] = read_h5py_file(H5FileLocAndName)
    if not inj_param_vals: # Raw files at first didn't record this, so make sure it's there...
        # Get the inj values from the processed data file instead
        DataFileLocAndName_NotRaw = DataFileLocAndName[:-7] + ".h5"
        [_,inj_param_vals,_,_] = read_h5py_file(DataFileLocAndName_NotRaw)
    labels = ["lambda","beta"]
    if np.size(infer_params[labels[0]][0])>1:
        nsamples = len(infer_params["lambda"][0])
    else:
        nsamples = len(infer_params["lambda"])
    # print("Posterior sample length: " + str(nsamples))
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

    return TotArea

def hist_lam_bet(data,lambda_bins,beta_bins):
    """
    Helper function used in `GetTotSkyAreaFromPostData_OLD()' because I didn't
    trust the 2D histogram functions I found online.

    It moves through lambsa bins, taking out all rows in data where the lambdas values
    are within the current lambda bin, and then 1D histograms the beta data to get the
    entries for that caolumn.

    Note that the function runs with uneven bins also and returns both the bin population
    and the area of the [ii,jj] bin. This might be helpful if ever a skymap is required
    outside of healpy or of unequal binnings (?).

    PARAMS
    ------
        - data :: np.array [nsamples by n_inferred_params]
            This can be extracted from lisabeta output h5 files and handed to this function.
            Note though that we assume here only that lambda is at index 0 and
            beta is at index 1, e.g. posterior_lambdas = data[:,0].
        - lambda_bins :: np array or list of bins for lambda to be counted over
        - beta_bins :: np array or list of bins for beta to be counted over

    OUTPUT
    ------
        - hist2D_pops == np.array [nbins_lambda by nbins_beta]
            Bin populations according only to the lambda and beta bins input at
            start of function.
        - areas == np.array [nbins_lambda by nbins_beta]
            Bin areas in case we have a list of arrat of uneven bins.
    """
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

def GetOctantLikeRatioAndPostProb(FileName):
    """
    ################# This function is no longer frequently used #################

    Function used in skymode project where the likelihood function only was used
    to estimate how many modes would be present in a posterior. This project was
    explored briefly, but we decided the correlations between number of posterior
    modes and estimated number of modes from likelihood function was not
    strong enough to warrent further exploration.

    This function is left here as it demonstrates how to access the h5 files and
    extract the likelihoods and information on skymodes.

    PARAMS
    ------
        - FileName :: string
            File name for either Json OR h5 file where the two have congruent names and
            subarchitecture in the relevant folders for Json and H5 files
            ('inference_param_files' and 'inference_data' respectively).
    OUTPUT
    ------
        - lnLikeRatios :: dictionary with tuple keys corresponding to skymodes
            Log likelihood of each skymode. Can use this to subtract other skymodes
            as in commented out line in order to estimate the ratios between relevant
            skymodes. See manual section on preliminary studies for skymodes for
            details of which likelihood ratios were tested.
        - OctantPostProbs
            Posterior probability for each octant labeled by the skymodes. This is
            what each log likelihood ratio is compared to in order to test
            validity and strength of correlations.

    ################# This function is no longer frequently used #################
    """
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

def read_h5py_file(FileName, path=SYNEX_PATH):
    """
    Function to read a lisabeta output H5file.

    PARAMS
    ------
        - FileName :: String
            Name and sub-architecture of either JsonFile or H5File located in `inference_param_files'
            or `inference_data' respectively.
    OUTPUT
    ------
        - infer_params :: Dict
            Dictionary of posteriors with key for each parameter inferred by lisabeta
        - inj_param_vals :: Dict
            Dictionary of ALL source parameters injected into run. Should have keys
            corresponding to injected values for LISA frame and SSB frame.
        - static_params :: Dict
            Dictionary of values of all parameters NOT included in inference run.
        - meta_data :: Dict
            Dictionary of extra information for the run, including acceptance ratios
            etc. that can be used for e.g. convergence tests.
    """

    # Check filenames
    JsonFileLocAndName, H5FileLocAndName = CompleteLisabetaDataAndJsonFileNames(FileName,
                                                                                path)

    # Load data from h5 file
    f = h5py.File(H5FileLocAndName,'r')

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
            elif key in list_params:
                infer_params[key] = param_vals
            else:
                meta_data[key] = param_vals
    # f.close()
    return [infer_params, inj_param_vals, static_params, meta_data]

def RunInference(source_or_kwargs, detector, inference_params, PlotInference=False,**RunTimekwargs):
    """
    Function to handle lisabeta inference, automatically detecting and handling
    parallel processing capabilities (e.g. installation and availability of mpi4py).

    The required Json fiel to run the lisabeta binary file is created automatically
    here using helper functions to convert classes to the lisabeta-compatible
    dictionaries and then write those dictionaries to Json file. The lisabeta binary is
    then called.

    PARAMS
    ------
        - source_or_kwargs :: SYNEX source class *OR* kwargs Dictionary
            We can initiate with either a SYNEX source class, or a dictionary to
            initiate SYNEX source class with.
        - detector         :: SYNEX detector class
        - inference_params :: Dict
            Lisabeta ready dictionary of inference parameters along with their
            prior ranges, prior distribution types (e.g. uniform), the range to sample
            from, and an extra legacy parameter no longer used. See './test_functions/installation_tests'
            notebook for an example of this dictionary.
        - PlotInference    :: Bool [Default False]
            Flag to request a corner plot of the inference results upon completion.
            Note this should not be used on the cluster unless user is sure to not
            exceed cluster space limits.
        - RunTimekwargs    :: Dict
            Extra parameters to replace within the source dictionary or inference
            params dictionary. This option is no longer used but retained in case
            futher developements need this option.

    NOTE: This function used to live in a seperate file ./SYNEX/SYNEX_PTMC.py, but
          was removed to here as something much more simple in order to remove
          problems keeping up with lisabeta development. Now we simply call the (new)
          binary `/lisabeta/lisabeta/inference/ptemcee_smbh.py' since this is maintained
          better and we will no longer require updating a second script in order to
          maintain congruence with lisabeta updates.
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

        # Check times match and adjust times if they don't match
        det_mission_times_gps = [detector.gps_science_start, detector.gps_science_start+detector.mission_duration*365.25]
        source_times_gps = [source.gps_timetomerger_max, source.gps_timetomerger_max+source.timetomerger_max*365.25+source.DeltatL_cut/(60.*60.*24.)]
        if source_times_gps[0]<det_mission_times_gps[0]:
            source.gps_timetomerger_max = detector.gps_science_start
            source.timetomerger_max = source.timetomerger_max - (det_mission_times_gps[0] - source_times_gps[0])/365.25
        if source_times_gps[1]>det_mission_times_gps[1]:
            source.DeltatL_cut = source.DeltatL_cut - (source_times_gps[1] - det_mission_times_gps[1])*24.*60.*60.

        # Write params to json file
        print("Creating json file...")
        WriteParamsToJson(source, detector, inference_params,
                          is_master, **RunTimekwargs)
        sourceJsonFile=source.JsonFile
        print("Done.")
    else:
        sourceJsonFile=None # So we only have master node running later if one source and/or detector given

    # Send source to workers
    if use_mpi:
        comm.Barrier()
        sourceJsonFile = comm.bcast(sourceJsonFile, root=0)

    # Start the run. Data will be saved to the 'inference_data' folder by default
    # All processes must execute the run together. mapper (inside ptemcee) will handle coordination between p's.
    # SYP.RunPTEMCEE(source.JsonFile)
    if is_master: print(" --------------- START INFERENCE --------------- ")
    print("Does Json exist?",MPI_rank,os.path.isfile(sourceJsonFile))
    t1 = time.time()
    command = "python3 " + SYNEX_PATH + "/lisabeta/lisabeta/inference/ptemcee_smbh.py " + sourceJsonFile
    os.system(command)
    t2 = time.time()
    if is_master: print(" ---------------- END INFERENCE ---------------- ")
    if is_master: print("Time to execute ptemcee: ", round((t2-t1)*10.)/10., "s")
    if is_master: print(" --------------------------------------------- ")

    # Create Skymap struct here. Commented out so that we don't automatically do this
    # on the cluster, for example, where space is limited.
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

def WriteParamsToJson(source, detector, inference_params,
                      IsMaster=True, **RunTimekwargs):
    """
    Function to save source and GW detector params to json file for
    lisabeta inference runs, using the saved json file name in source class.
    """
    # Only verbose if both classes are verbose
    verbose=source.verbose and detector.verbose

    # Chose the repertory where the files are saved. Use SYNEX path b y default``
    #  if nothing provided in the RunTimekwargs parameters.
    if "out_rep" in RunTimekwargs:
        out_rep = RunTimekwargs["out_rep"]
    else:
        out_rep = SYNEX_PATH

    # See if we have output filenames to asign
    # if "out_file" in RunTimekwargs:
    #     out_file = RunTimekwargs["out_file"]
    #     del RunTimekwargs["out_file"]

    #     JsonFile, H5File = CompleteLisabetaDataAndJsonFileNames(out_file,
    #                                                             out_rep)
    #     source.JsonFile = JsonFile
    #     source.H5File = H5File


    # double check names and paths are ok
    # if source.JsonFile==None and source.H5File==None:
    # No filename provided at all
    if (not source.JsonFile) and (not source.H5File):
        # Prepare the default filename
        from datetime import date
        today = date.today()
        d = today.strftime("%d_%m_%Y")
        JsonOrOutFile = (d + "_InferParams_"
                         + '_'.join(inference_params["infer_params"]))

        # Create the repertories
        JsonFile, H5File = CompleteLisabetaDataAndJsonFileNames(JsonOrOutFile,
                                                                out_rep)

        if verbose: print("Creating default json name: {}".format(JsonFile))

    # Case where only the json file is provided
    elif (source.JsonFile) and (not source.H5File):
        JsonFile, H5File = CompleteLisabetaDataAndJsonFileNames(source.JsonFile,
                                                                out_rep)
    # Case where only the H5 file is provided and not the json
    elif (not source.JsonFile) and (source.H5File):
        JsonFile, H5File = CompleteLisabetaDataAndJsonFileNames(source.H5File,
                                                                out_rep)
    # Case where both json and h5 filename are given
    elif source.JsonFile!=None and source.H5File!=None:
        JsonFile, H5File = CompleteLisabetaDataAndJsonFileNames(source.JsonFile,
                                                                out_rep)
        if source.H5File!=H5File and verbose:
            print("Json and H5 filenames are not matched... \n"
                  "There is currently no check that the paths are ok here.\n"
                  "Make sure to pass the entire path if doing this.")

    # source.JsonFile = Jsons

    # import some default parameters defined the ptemcee handler script
    # from lisabeta.inference.ptemcee_smbh import run_params_default # , waveform_params_default
    run_params_default = {
        "sampler": "ptemcee",
        "sample_Lframe": True,
        "multimodal": True,
        "multimodal_pattern": "8modes",
        "p_jump": 0.5,
        "ensemble_proposal": "ptemcee",
        "likelihood_method": "fresnel",
        "likelihood_residuals_ngrid": None,
        "skip_fisher": False,
        "init_method": "fisher",
        "init_file": None,
        "init_scale_cov": 1.,
        "zerolike": False,
        "n_temps": 10,
        "temp_max": None,
        "adaptive_temp": False,
        "n_walkers": 96,
        "n_iter": 8000,
        "burn_in": 5000,
        "autocor_method": "autocor_new",
        "thin_samples": True,
        "seed": None,
        "print_info": True,
        "n_iter_info": 10,
        "output": True,
        "output_raw": False,
        "upsample": 1, ## Extra parameter not sure what its for...
        "params_map": None ## Extra parameter not sure what its for...
    }

    # Create the default json file parameters
    json_default_dict = {}
    json_default_dict["run_params"] = run_params_default

    # Get additional parameters and add to json param dict
    # these dictionaries are 1. Base source params 2. waveform params 3. extras for fisher function
    # Therefore can rename the dictionaries for clarity, and then delete the key,value unwrap next
    # for the run_params update, since the run_params are handed to the function at call time.
    # But double check that the run_params are not not in the returned waveform params list etc.
    [param_dict, waveform_params, _ ] = ClassesToParams(source,
                                                        detector,
                                                        "Inference")
    json_default_dict["source_params"] = param_dict
    json_default_dict["waveform_params"] = waveform_params
    json_default_dict["prior_params"] = inference_params

    # Add some missing kwargs in the LISAnoise subdictionary in waveform params.
    # Not really sure what this does but it is needed.
    # Need to fiure out if this is needed in theFisher stuff too - it will
    # change where we put the defaults etc at detector initialization.
    json_default_dict["waveform_params"]["LISAnoise"]["lowf_add_pm_noise_f0"] = 0.0
    json_default_dict["waveform_params"]["LISAnoise"]["lowf_add_pm_noise_alpha"] = 2.0

    # Update json dict field defaults where needed
    for key,value in waveform_params.items():
        if key in json_default_dict["run_params"]:
            json_default_dict["run_params"][key] = value # these are not in the
        # elif key in json_default_dict["waveform_params"]:
        #     json_default_dict["waveform_params"][key] = value

    # Change now any keys set in the run time dictionary of kwargs (this could plot flags, run params values, no. of walkers etc)
    # This is NOT meant for binary params or prior params - these need to be specified at the highest script level
    for key, value in RunTimekwargs.items():
        if key in json_default_dict["run_params"]:
            json_default_dict["run_params"][key] = value
        elif key in json_default_dict["waveform_params"]:
            json_default_dict["waveform_params"][key] = value
        elif key in json_default_dict["source_params"]:
            json_default_dict["source_params"][key] = value

    # Set the output file and directory location in the json param list
    H5FilePath = "/".join(source.H5File.split("/")[:-1])+"/"
    H5FileName = source.H5File.split("/")[-1]
    H5FileName=".".join(H5FileName.split(".")[:-1])
    json_default_dict["run_params"]["out_dir"] = H5FilePath
    json_default_dict["run_params"]["out_name"] =  H5FileName # this needs to not have the '.h5' added at the end to work

    # Write the json file only if master node or not mpi
    if IsMaster:
        if verbose: print("Writting params to", source.JsonFile)
        with open(source.JsonFile, 'w') as f:
            json.dump(json_default_dict, f, indent=2)
        f.close()
    else:
        time.sleep(10) # MPI.COMM_WORLD.Barrier()???

def RunFoMOverRange(source,detector,ParamDict,FigureOfMerit='SNR',RunGrid=False,inference_params=None,**InferenceTechkwargs):
    """
    ################# This function is no longer frequently used #################

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

    ################# This function is no longer frequently used #################
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
    """
    Helper function used pricipally in `RunFoMOverRange()' (which is no longer
    frequenctly used).

    Given a source and GW detector, log likelihood and parameters for each each
    skymode is returned.

    PARAMS
    ------
        - source           :: SYNEX source class
        - detector         :: SYNEX detector class
    OUTPUT
    ------
        - lnL_skymodes    :: Dict
            Dictionary with key for each skymode, giving the loglikelihood value
            for each skymode octant.
        - params_skymode  :: Dict
            Dictionary with key for each skymode, giving the parameters describing
            the skymode (e.g. the central lambda and beta values).

    NOTE: Although the main function `RunFoMOverRange()' is no longer frequently used,
          this function could be useful in the future depending on project ideas.
          Hence it is left here. It is also useful for information on contents of
          lisabeta files.
    """
    # Get the parameters out of the classes and assign if given in kwargs dict
    [param_dict, waveform_params, extra_params] = ClassesToParams(source,detector,"Inference",**kwargs)

    # Define the likelihood class if not already existant
    likelihoodClass = lisa.LikelihoodLISASMBH(param_dict, **waveform_params)
    return lisatools.func_loglikelihood_skymodes(likelihoodClass)

def GetSkyMultiModeProbFromJson(FileName):
    """
    Same purpose as `GetSkyMultiModeProbFromClasses()', but starting from either
    Json or H5 file (including sub-architecture under main folder `inference_param_files'
    or `inference_data' respectively).

    PARAMS
    ------
        - FileName :: String
            Name and sub-architecture of either JsonFile or H5File located in `inference_param_files'
            or `inference_data' respectively.
    OUTPUT
    ------
        - lnL_skymodes    :: Dict
            Dictionary with key for each skymode, giving the loglikelihood value
            for each skymode octant.
        - params_skymode  :: Dict
            Dictionary with key for each skymode, giving the parameters describing
            the skymode (e.g. the central lambda and beta values).

    NOTE: Although not (currently) used in any other SYNEX function, this function
          could be useful in the future depending on project ideas. Hence it is left
          here. It is also useful for information on contents of lisabeta files.
    """

    # Check filenames
    JsonFileLocAndName,H5FileLocAndName = CompleteLisabetaDataAndJsonFileNames(FileName)

    # Read contents of file
    with open(JsonFileAndPath, 'r') as input_file: input_params = json.load(input_file)

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
    """
    Function to grab the fisher cov matrix from lisabeta.

    PARAMS
    ------
        - Source : SYNEX source object
        - Detector : SYNEX detector object
        - kwargs : Dict
            Dictionary of values to replace source or detector values quickly in
            calculations of fisher matrix. Note that if a parameter is given here,
            the function will load the source and detector params to lisabeta dictionaries,
            and then replace the corresponding value within those dictionaries. It will
            *NOT* replace the values within the source and or detector classes that will remain
            unchanged.
    """
    # Get the parameters out of the classes and assign if given in kwargs dict
    [param_dict, waveform_params, extra_params] = ClassesToParams(source,detector,"Fisher",**kwargs)

    # Call the fisher function
    return lisa_fisher.fisher_covariance_smbh(param_dict, **waveform_params)

def GetSMBHGWDetection(source, detector, **kwargs):
    """
    Wrapper function to get the measured (TDI) waveform from a SYNEX source class
    and SYNEX GW detector class.

    PARAMS
    ------
        - Source : SYNEX source object
        - Detector : SYNEX detector object
        - kwargs : Dict
            Dictionary of values to replace source or detector values quickly.
            Note that if a parameter is given here, the function will load the
            source and detector params to lisabeta dictionaries, and then replace
            the corresponding values within those dictionaries. It will *NOT*
            replace the values within the source and or detector classes that
            will remain unchanged.
    OUTPUT
    ------
        - wftdi : Dict
          Dictionary with a key for each GW harmonic with a sub-dictionary
          containing 'freq' (frequencies of the waveform), 'amp' (amplitude), and
          'phase'.  'modes' (a list of all the modes) is also contained in 'wftdi'
    """
    # Get the parameters out of the classes and assign if given in kwargs dict
    [param_dict, waveform_params, extra_params] = ClassesToParams(source, detector, "Fisher", **kwargs)

    # This is a workaround for flexibility with the input dicts. Not sure how fast it is- to be improved later!
    s = ' '
    for key,value in waveform_params.items():
        if isinstance(value, str):
            s = s + ', ' +  key + '="' + str(value) + '"'
        else:
            s = s + ', ' +  key + '=' + str(value)
    return eval('lisa.GenerateLISATDI_SMBH(param_dict' + s + ',extra_params=extra_params)')

def ComputeSNR(source, detector, ReturnAllVariable=False, **kwargs):
    """
    Wrapper function to get the SNR of a source as detected by a GW detector.

    Old version of this function is kept for the time being but will be removed
    soon. Older version had extra parameters 'freqs'=None, 'Lframe'=False, and
    'ReturnAllVariable'=False.

    PARAMS
    ------
        - Source :: SYNEX source object
        - Detector :: SYNEX detector object
        - kwargs :: Dict
          Extra dictionary to alter any parameters inside either source or detector
          class. This makes it easier to quickly get the SNR of a large set of
          objects without having to manually create all the versions of each object.
    OUTPUT
    ------
        - SNR :: float
    """
    # # Sort into dictionaries params, freqs
    # [params, waveform_params, extra_params] = ClassesToParams(source,detector,"Fisher",**kwargs)
    #
    # if Lframe:
    #     raise ValueError('Lframe not implemented yet, sorry.')
    #
    # # Parameters
    # m1 = params['m1']
    # m2 = params['m2']
    # M = m1 + m2
    # q = m1 / m2
    # params['M'] = M
    # params['q'] = q
    #
    # LISAnoise = waveform_params.pop('LISAnoise', pyLISAnoise.LISAnoiseSciRDv1)
    # TDI = waveform_params.get('TDI', 'TDIAET')
    # TDIrescaled = waveform_params.get('TDIrescaled', True)
    # LISAconst = waveform_params.get('LISAconst', pyresponse.LISAconstProposal)
    #
    # # Default frequencies (safe and slow):
    # # linear spacing, adjust deltaf=1/(2T) with T approximate duration
    # if freqs is None:
    #     freqs = ['linear', None]
    # if not isinstance(freqs, np.ndarray):
    #     freqs = GenerateFreqs(freqs, params, **waveform_params)
    #
    # # Generate tdi freqseries
    # tdifreqseries_base = lisa.GenerateLISATDIFreqseries_SMBH(params, freqs, **waveform_params)
    #
    # # Compute noises
    # noise_evaluator = pyLISAnoise.initialize_noise(LISAnoise,
    #                                     TDI=TDI, TDIrescaled=TDIrescaled,
    #                                     LISAconst=LISAconst)
    # Sn1_vals, Sn2_vals, Sn3_vals = pyLISAnoise.evaluate_noise(
    #                       LISAnoise, noise_evaluator, freqs,
    #                       TDI=TDI, TDIrescaled=TDIrescaled, LISAconst=LISAconst)
    #
    # # Modes and channels
    # modes = tdifreqseries_base['modes']
    # channels = ['chan1', 'chan2', 'chan3']
    # Snvals = {}
    # Snvals['chan1'] = Sn1_vals
    # Snvals['chan2'] = Sn2_vals
    # Snvals['chan3'] = Sn3_vals
    # h_full = {}
    # loopCount=0
    # for chan in channels:
    #     h_full[chan] = np.zeros_like(freqs, dtype=complex)
    #     for lm in modes:
    #         h_full[chan] += tdifreqseries_base[lm][chan]
    #         loopCount+=1
    #
    # # Compute SNR
    # # We do not assume that deltaf is constant
    # df = np.diff(freqs)
    # SNR = 0.
    # for chan in channels:
    #     SNR += 4*np.sum(df * np.real(h_full[chan] * np.conj(h_full[chan]) / Snvals[chan])[:-1])
    #
    # if ReturnAllVariable:
    #     return np.sqrt(SNR), freqs, h_full, Snvals
    # else:
    #     # wftdi = GetSMBHGWDetection(source, detector, **kwargs)
    #     # return lisa.computeSNR(wftdi) # ,LISAconstellation=pyresponse.LISAconstProposal,LISAnoise=pyLISAnoise.LISAnoiseProposal,LDCnoise=None,TDIrescaled=False,Nf=None)
    #     return np.sqrt(SNR)

    # Get the parameters out of the classes and assign if given in kwargs dict
    [param_dict, waveform_params, _] = ClassesToParams(source,detector,"Fisher",**kwargs)

    # Call new snr tools wrapper in lisabeta
    SNR = lisa_snr.lisa_mbhb_snr(param_dict, time_to_merger_cut=None, **waveform_params)

    # Should we include option to return all variables as in old function here? then we have to pull out the inner workings of this function...
    # In this state the input variable `ReturnAllVariable' is redundant.
    return SNR

def GenerateFreqs(freqs, params, **waveform_params):
    """
    Helper function to automatically generate the frequencies needed for caluclations
    like SNR and detector noise (see for example ComputeDetectorNoise()).

    Params
    ------
        - freqs : list
          Object to specify what kind of frequencies is required. For example:
          freqs = ['linear', None] requests linear spacing, with no pre-required
          bounds. The generator will find the total time of the signal from the
          params and waveform_params dictionaries, and then generate a frequeny
          range that covers the full waveform.

        - params : dict
          Dictionary of parameters describving the source (mass, spin etc).

        - waveform_params : dict
          Dictionary of parameters needed for waveform generation (approximant,
          modes, etc.)
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

def ComputeDetectorNoise(source, detector, freqs=None, Lframe=False, **kwargs):
    """
    Function adapted from lisabeta lisa.py (I think) to compute the sensitivity
    curve of a GW detector object. A source object must be given so that the
    noise generator knows which frequencies to calculate the noise curves over.
    NOTE: LISA noise has THREE channels corresponding to the three time delay
    interferometry (TDI) responses to a signal, labelled 'chanx' for x={1,2,3}.

    PARAMS
    ------
        - Source :: SYNEX source object
        - Detector :: SYNEX detector object
        - freqs :: list, ndarray or None [default None]
            Frequencies over which to calculate noise. Can be one of following:
            1. ['log', npt]
                Log-spaced bins with npt bins. Set npt=None to autospace with T and
                source params.
            2. ['linear', deltaf]
                Linear-spaced bins with binwidth deltaf. deltaf=None for autospacing
                with T and source params.
            3. ndarray
                Skip generation of freqs and use array of frequencies given instead.
                Note that this option is not well testing in SYNEX.
            4. None
                The default option. Then a linear binning is generated with default
                deltaf generated based on T and source params.
        - Lframe :: Bool [default False]
            Option to change source params to or from Lframe.
        - kwargs :: Dict
          Extra dictionary to alter any parameters inside either source or detector
          class. This makes it easier to quickly get detector noise over a large
          set of similar objects without having to manually create all the versions
          of each object.
    """
    # Grab the variables from classes
    [params, waveform_params, extra_params] = ClassesToParams(source,detector,"Inference",**kwargs)

    LISAnoise = waveform_params.get('LISAnoise', pyLISAnoise.LISAnoiseSciRDv1)
    TDI = waveform_params.get('TDI', 'TDIAET')
    TDIrescaled = waveform_params.get('TDIrescaled', True)
    LISAconst = waveform_params.get('LISAconst', pyresponse.LISAconstProposal)

    # Convert input parameters to either Lframe or SSBframe
    if not waveform_params.get('Lframe', False) and Lframe:
        params_base = lisatools.convert_SSBframe_to_Lframe(
                            params,
                            t0=waveform_params['t0'],
                            frozenLISA=waveform_params['frozenLISA'])
    elif waveform_params.get('Lframe', False) and not Lframe:
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
    Snvals['freqs'] = freqs

    detector.Snvals = Snvals

def fmaxFromTimeToMerger(source, detector, T_obs_end_to_merger=4.*60.*60., ReturnVal=False, **kwargs):
    """
    ################# This function is no longer frequently used #################

    Legacy function from previous versions of GW sensitivity calculations. This
    is left though to show how to calculate the maximum frequency of a GW waveform
    and a detector noise curve based on source masses.

    ################# This function is no longer frequently used #################
    """
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
    """
    Function to turn lisabeta file into a source class. Lisabeta inference
    does * NOT * need to have been run for this function to work. All that needs
    to exist is a Json file containing all of the source parameters. If the H5File
    exists, then the class init function will generate the skymap file etc.

    PARAMS
    ------
        - FileName
            Either Json filename or H5 filename, with or without file extensions,
            to load data from. Only substructure below either 'inference_param_files'
            or 'inference_data' main folders is needed and the function will assume
            local SYNEX architecture.
        - kwargs : Dict
            Dictionary of values to replace source or detector values once classes
            have been confiured from base data files. Not sure that this is really
            useful if we don't then change the data itself too...
    """

    # Check filenames
    JsonFileLocAndName,H5FileLocAndName=CompleteLisabetaDataAndJsonFileNames(FileName)

    # Extract data
    with open(JsonFileLocAndName, 'r') as f: input_params = json.load(f)

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


def TileSkyArea(CloningTrackFile=None,sources=None,telescopes=None,base_telescope_params=None,cloning_params=None,SaveInSubFile=None,SaveFileCommonStart=None,SourceIndexingByString=None,verbose=True):
    """
    MASTER function to tile skyarea. This will handle all extra tiling techniques
    we might add in the future outside of GWEMOpt. This function is stupid large
    and needs to be broken into helper functions to help legibility and debugging
    in future developments.

    A list of sources or a single source is given, along either EITHER a list of
    telescopes (or single telescope) OR a dictionary of base_telescope_params and a
    dictionary of cloning_params to spawn many identical telescopes with a single
    parameter changed at a time.

    GWEMOpt is then run on all possible combinations of source and telescope, either
    returning a list of all output telescopes if no parallel processing is used
    (sensed via presence of mpi4py), or returning nothing but saving all output
    telescopes to their specified savefiles.

    Note that if we have specified a CloningTrackFile then this is preferred over
    all other possible initiation methods (i.e. if cloning_params is also specified,
    it will be ignored in favour of loading contents of CloningTrackFile).

    PARAMS
    ------
        - CloningTrackFile :: string [default None]
            Complete file path to file (with extension) that keeps track of all
            source and telescope combinations. This is useful if you are running
            a very large list of combinations that cannot all be run at the same
            time on the cluster (because space is restricted). All of the combinations
            are saved to a file that includes a line for each combination (comb)
            that includes the source save file, the telescope save file, and the
            values of each parameter that was changed between each combination.

            NOTE: When cloning_params is used to clone a base telescope, the
            combinations are stored in the file "./SYNEX/TileTrackFiles/ProgressTrackFile.txt".

            *** FUTURE DEV: CHANGE this argument from a string to a bool that tells the
            *** function to prefer the savefile or not, and force the savefile to
            *** be "./SYNEX/TileTrackFiles/ProgressTrackFile.txt" every time.

        - sources :: list, SYNEX source class, or None [default None]
            A list of source classes, a list of source save files, a single source
            class, or a single source save file. File names should be complete with
            extension, and saved in the proper location: './SYNEX/Saved_Source_Dicts/'.

            NOTE: if a single source save file is given, this can take the form of
            a wild card search string too:
            './SYNEX/Saved_Source_Dicts/Test_System_*' will run GWEMOpt over a
            list of sources that are located within the saved source files that
            satisfy the search string.
        - telescopes :: list, SYNEX telescope class, or None [default None]
            A list of telescope classes, a list of telescope save files, a single telescope
            class, or a single telescope save file. File names should be complete with
            extension, and saved in the proper location: './SYNEX/Saved_Telescope_Dicts/'.

            NOTE: if a single telescope save file is given, this can take the form of
            a wild card search string too:
            './SYNEX/Saved_Telescope_Dicts/Test_System_*' will run GWEMOpt over a
            list of telescopes that are located within the saved telescope files that
            satisfy the search string.
        - base_telescope_params :: dict [default None]
            TO BE USED WITH 'cloning_params'

            When initiating tiling using the cloning params options, we start with
            a base telescope and clone it, changing only one
            parameter at a time.
            NOTE: the base telescope class or savefile could also be handed to
            this function via the 'telescopes' argument while setting
            base_telescope_params=None.
        - cloning_params :: dict [default None]
            TO BE USED WITH 'base_telescope_params'

            A dictionary with a list of values to implement on each key, where
            each key corresponds to an attribute in the telescope class object.
            A list of input telescopes is then created from the cloned telescopes
            classes, after which the function runs as usual by calculating all
            possible combinations of input sources with cloned telescopes and passing
            the combinations, one by one, to GWEMOpt.
            NOTE: the base telescope class or savefile could also be handed to
            this function via the 'telescopes' argument while setting
            base_telescope_params=None.
        - SaveInSubFile :: string [default None]
            Sub-structure under './SYNEX/Saved_Telescope_Dicts/' to save all telescope
            files. This is useful if we clone a large number of telescopes and want
            to group them into a file structure together.
        - SaveFileCommonStart :: string [default None]
            The starting string for all telescope savefiles (e.g. the name of the
            current experiment).
        - SourceIndexingByString :: string [default None]
            String to use when labelling sources in telescope save files. For example
            if we have a list of randomized sources with savefiles
            './SYNEX/Saved_Source_Dicts/System_XX.dat', where 'XX' is a number
            unique to each source, then we will wish to include the labelling in
            the telescope savefile names to keep track of which source was used
            where. To do this, we set 'SourceIndexingByString:"System_"', and the
            telescope savefiles will then include 'XX' in each savefile name.

            If this is set to None (as is the default) then the systems are labelled
            by a numbering from one to the number of input sources.

            NOTE: when GWEMOpt tiling is run, the resulting telescope class has a
            record of which source was used, e.g. :

            telescopes[0].telescope_source_coverage["source save file"]=source.ExistentialFileName
            telescopes[0].telescope_source_coverage["source JsonFile"]=source.JsonFile
            telescopes[0].telescope_source_coverage["source H5File"]=source.H5File

            However, this options means the labeling is in the telescope savefile name
            so that specific sources can be located easily without having to load
            load all telescopes at once and then remove those that aren't associated
            with the required source.
        - verbose :: bool [default True]
            Option to include updates through the tiling procedure.

    OUTPUT
    ------
        - telescopes_out :: list ----- ONLY IF WE HAVE ONE PROCESSOR E.G. NO MPI
            If there is one processor active, e.g. running locally without mpi4py
            (we are not running in parallel) then we return a list of telescope
            classes with the completed tiling output stored as class attributes for
            each telescope in the list of output telescopes. The number of output
            telescopes corresponds to the total number of combinations between input
            telescopes and sources. E.g. if we input 5 telescopes as a list, and
            2 sources in a list, with all other inputs set to None, we will have
            an output list of 10 telescopes (each input source run once with each
            input telescope).
            If there is parallel processing, then no output is provided. Instead,
            each completed telescope is saved to it's savefile and NO OTHER INFORMATION
            (e.g. tesselation file) IS SAVED. This is to save workspace and memory
            when running on cluster. Everything else is the same as a single processor
            output case: if we begin with 5 telescopes and 2 sources at input, we
            expect 10 total savefiles in './SYNEX/Saved_Telescope_Dicts/' corresponding
            to the total number of combinations between input sources and input
            telescopes.

    TO DO / OPEN QUESTIONS :
    ------------------------
    1. Should we allow case where 'DeltatL_cut' can be in cloning_params? Then we gotta check
       filenames exist etc... this would allow us to move source creations to loop one by one too... PRIORITIZE THIS!
    2. Reduce memory usage by loading base telecope by master only and transfering?
    3. Change construction and handling of tile tracking function to dataframe so we conserve organized data
       and de-risk mixing combinations with save files. In it's current form (brute force saving to txt file) the permutations
       and savefiles for all objects align again within the save file but we don't know if this is always the case across all
       versions of mpi (because we use mpi_comm.gather to collect permuation info)
    4. Add changes to default go_params flags (Ncores etc) for case of single detector and single
       source input with MPI_size>1. Can we go one step further and calculate a share of cores
       to go into parallelizing gwemopt with remaining cores used to share input objects across.
    """
    # Using MPI or not
    if MPI is not None:
        MPI_size = MPI.COMM_WORLD.Get_size()
        MPI_rank = MPI.COMM_WORLD.Get_rank()
        comm = MPI.COMM_WORLD
    else:
        MPI_size=1
        MPI_rank=0

    # See if we have a savefile for progress to pick up from... -- Change this later to dataframe which would be easier using indexing for combination labeling
    # Otheriwise when we create the combinations file by gathering each node's set of cobs we risk de-prganizing things..
    # check for a default trackfile ? Would reduce input params...
    if CloningTrackFile:
        if MPI_rank==0:
            with open(CloningTrackFile, 'r') as f: lines=f.readlines()
            SaveInSubFile=lines[0]
            CloningKeys=lines[1].split(",")
            BaseTelescopeName=lines[2]
            BaseTelescopeExFileName=lines[3]
            data_all = lines[4:]
            Nvals=len(data_all)
        else:
            SaveInSubFile=None
            BaseTelescopeName=None
            BaseTelescopeExFileName=None
            Nvals=None

        # Bcast static data across cores                    ##### CHECK EACH CORE HAS THE RIGHT VALS !!!
        if MPI_size>1:
            SaveInSubFile = comm.bcast(SaveInSubFile,root=0)
            if isinstance(SaveInSubFile,str):
                SaveInSubFile=SaveInSubFile.strip('\n')
            if SaveInSubFile=="None":
                SaveInSubFile=None
            BaseTelescopeName = comm.bcast(BaseTelescopeName,root=0).strip('\n')
            BaseTelescopeExFileName = comm.bcast(BaseTelescopeExFileName,root=0).strip('\n')
            Nvals = comm.bcast(Nvals,root=0)

        # Scatter data across cores                         ##### CHECK EACH CORE HAS THE RIGHT VALS !!!
        if MPI_size>1:
            nPerCore=int(Nvals//MPI_size)
            remainder=Nvals-nPerCore*MPI_size
            if MPI_rank==0:
                data_ii_start=0
                data_ii_end=data_ii_start+nPerCore
                if remainder>0: data_ii_end+=1 ## Split any extras left from rounding errors evenly
                data = data_all[data_ii_start:data_ii_end]
                for core_ii in range(1,MPI_size):
                    data_ii_start=data_ii_end
                    data_ii_end=data_ii_start+nPerCore
                    if core_ii<remainder: data_ii_end+=1 ## Split any extras left from rounding errors evenly
                    comm.send(CloningKeys, dest=core_ii)
                    comm.send(data_all[data_ii_start:data_ii_end], dest=core_ii)
            else:
                CloningKeys = comm.recv(source=0)
                CloningKeys = [k.strip('\n') for k in CloningKeys]
                data = comm.recv(source=0) ### strip '\n' from these and bcast vals too.
        else:
            data = data_all
        CoreLenVals = len(data)                           ##### CHECK EACH CORE HAS THE RIGHT VALS !!!

        # Make base detector params dict
        if base_telescope_params==None: base_telescope_params=SYDs.Athena(**{"ExistentialFileName":BaseTelescopeExFileName,"verbose":verbose}).__dict__

        # Configure data
        TelescopeNewExNames=[]
        SourceNewExNames=[]
        CloningCombs=[]
        for line in data:
            linelist=line.split(":")
            TelescopeNewExNames.append(linelist[0])
            SourceNewExNames.append(linelist[1])
            Comb=[float(el) if el!="None" else None for el in linelist[2].split(",")]
            if "Tobs" in CloningKeys:
                Comb=[v if k!="Tobs" else np.array([0.,v]) for k,v in zip(CloningKeys,Comb)]
            CloningCombs.append(Comb) # What if we have a string? Like a flag or 'None' value?

        # Create sources
        sources=[SYSs.SMBH_Merger(**{"ExistentialFileName":ExName,"verbose":verbose}) for ExName in SourceNewExNames]

    elif cloning_params!=None:
        # Create base_telescope_params from whatever was handed to us, otherwise obtain from default vals
        if base_telescope_params==None:
            if telescopes==None:
                base_telescope_params=SYDs.Athena().__dict__
            elif isinstance(telescopes,list):
                if len(telescopes)>1 and verbose:
                    print("Cloning first telescope in list only...")
                base_telescope_params=telescopes[0].__dict__
            else:
                base_telescope_params=telescopes.__dict__

        # Base detector savefile and name
        BaseTelescopeExFileName=base_telescope_params["ExistentialFileName"] if "ExistentialFileName" in base_telescope_params else SYNEX_PATH+"/Saved_Telescope_Dicts/Athena_Base.dat"
        if 'detector_config_struct' in base_telescope_params:
            BaseTelescopeName=base_telescope_params["detector_config_struct"]["telescope"]
        elif "telescope" in base_telescope_params:
            BaseTelescopeName=base_telescope_params["telescope"]
        else:
            BaseTelescopeName="Athena"

        # Check if we have a 'sources' list as input
        if isinstance(sources,dict) and "ExistentialFileName" in sources:
            if sources["ExistentialFileName"][-1]=="/": source+="*"
            sources=glob.glob(sources["ExistentialFileName"]) if "*" in sources["ExistentialFileName"] else [sources["ExistentialFileName"]]
            cloning_params["SourceExName"]=sources
        elif isinstance(sources,dict) and "H5File" in sources:
            if sources["H5File"][-1]=="/": source+="*"
            sources=glob.glob(sources["H5File"]) if "*" in sources["H5File"] else [sources["H5File"]]
            cloning_params["H5File"]=sources
        elif isinstance(sources,str):
            if sources[-1]=="/": sources+="*"
            sources=glob.glob(sources) if "*" in sources else [sources]
            if sources[0][-3:]==".h5":
                cloning_params["H5File"]=sources
            elif sources[0][-3:]=="dat":
                cloning_params["SourceExName"]=sources #### This has potential to be broken easily...
        elif isinstance(sources,list) and isinstance(sources[0],str): ### OR PATH?
            if sources[0][-3:]==".h5":
                cloning_params["H5File"]=sources
            elif sources[0][-3:]=="dat":
                cloning_params["SourceExName"]=sources #### This has potential to be broken easily...
        else: # handed list of sources or just one source
            cloning_params["SourceExName"]=[source.ExistentialFileName for source in sources] if isinstance(sources,list) else [sources.ExistentialFileName] ## should we delete the list here or find a way to check later if it doesn't already exist?

        # Remove source related params from cloning dict since we replaced with source ExNames instead
        SourceCloneParams = ["Tcut","DeltatL_cut"] ## In case we add more later
        DO_TCUT=True if any([key in cloning_params for key in SourceCloneParams]) else False
        sourceTcut = cloning_params.pop("Tcut",None) ### Keep these seperate in case we need to access dictionary values in future developements
        sourceDeltatL_cut = cloning_params.pop("DeltatL_cut",None)

        # Need to know how many objects we will need --- later will include some source params in cloning params dict bith check that sources must == None to use cloning stuff
        CloningVals = list(cloning_params.values())
        CloningKeys = list(cloning_params.keys())
        CloningCombs = [list(Comb) for Comb in list(itertools.product(*CloningVals))] # returns list of tuples of all possible combinations
        Nvals = len(CloningCombs)

        # Work out how many items per cpu to reduce data usage asap
        NValsPerCore=int(Nvals//MPI_size)
        CoreLenVals=[NValsPerCore+1 if ii<Nvals%MPI_size else NValsPerCore for ii in range(MPI_size)]
        CPU_ENDs=list(np.cumsum(CoreLenVals))
        CPU_STARTs=[0]+CPU_ENDs[:-1]

        # Assign subset of combinations to each core
        CloningCombs = [CloningCombs[ii] for ii in range(CPU_STARTs[MPI_rank],CPU_ENDs[MPI_rank])]

        # Create list of sources
        if "H5File" in cloning_params:
            ExNames = [ParamComb[-1].replace("/inference_data/", "/Saved_Source_Dicts/").replace(".h5", ".dat") for ParamComb in CloningCombs]
            sources=[GetSourceFromLisabetaData(ParamComb[CloningKeys.index("H5File")],**{"ExistentialFileName":ExName,"verbose":verbose}) for ParamComb,ExName in zip(CloningCombs,ExNames)]
        elif "SourceExName" in cloning_params:
            sources=[SYSs.SMBH_Merger(**{"ExistentialFileName":ParamComb[CloningKeys.index("SourceExName")],"verbose":verbose}) for ParamComb in CloningCombs]

        # Get values of cloned source params
        if DO_TCUT:
            if "H5File" in CloningKeys:
                CloningKeys[CloningKeys.index("H5File")]="Tcut"
            elif "SourceExName" in CloningKeys:
                CloningKeys[CloningKeys.index("SourceExName")]="Tcut"
            SourcePVals=[-s.DeltatL_cut/86400. for s in sources]
            CloningCombs=[[v if k!="Tcut" else SourcePVals[iparam] for v,k in zip(comb,CloningKeys)] for iparam,comb in enumerate(CloningCombs)] # Inner list comp here is wried. Redo using list[list.index(val)]=NewVal maybe? but I think we are also restructuring here... In which case does this carry over when we don't have a source param cloned???

        # Add structure to detector savefile location and name if requested
        FolderArch = "/".join(BaseTelescopeExFileName.split("/")[:-1])
        if SaveInSubFile and SaveInSubFile[0]!="/": SaveInSubFile="/"+SaveInSubFile
        if SaveInSubFile: FolderArch += SaveInSubFile
        if FolderArch[-1]=="/": FolderArch=FolderArch[:-1]
        pathlib.Path(FolderArch).mkdir(parents=True, exist_ok=True)
        SaveFileCommon=FolderArch+"/"+SaveFileCommonStart if SaveFileCommonStart else FolderArch+"/"+BaseTelescopeExFileName.split("/")[-1].strip(".dat")

        # Get system IDs if we want them i.e. if we randomized some params and want to organize accordingly
        if SourceIndexingByString:
            SourceIndices=[source.ExistentialFileName.split(SourceIndexingByString)[-1] for source in sources]
            SourceIndices=[s[:s.rindex("_")] for s in SourceIndices]
        else:
            SourceIndices=[str(i) for i in range(1,len(sources)+1)]

        # Create list of detector dicts NB:: exfile names depend on what was changed... Make sure to treat long numbers so we dont get stupid savefile names.
        ### I am sure this can be optimized... ### -- Probably better to turn CloningCombs into a numpy array...
        # Seee if we changed the sources first...
        TelescopeNewExNames=[]
        for ii,ParamComb in enumerate(CloningCombs):
            # New Existential filename
            KeyValStrings=[]
            for k,v in zip(CloningKeys,ParamComb):
                if isinstance(v,np.ndarray):
                    vs=str(round(v[-1],2)) ### What if we have gaps in Tobs???? should be count number of gaps and put that in too?
                elif k=="Tcut":
                    if v<1:
                        vs=v*24 # Turn into hours... Never go below 1 hour anyway
                        vs=str(round(vs,2))+"hr"
                    else:
                        vs=str(round(v,2))+"d"
                elif isinstance(v,(float,int)):
                    vs=str(round(v,2))
                elif isinstance(v,bool):
                    # incase we start playing with flags later?
                    vs=str(int(v))
                if k!="SourceExName": KeyValStrings.append("{key}_{val}".format(key=k,val=vs))

            # Use pairings with source IDs if given
            if SourceIndexingByString:
                TelescopeNewExNames.append(SaveFileCommon+"_SourceInd_"+SourceIndices[ii]+"__"+"_".join(KeyValStrings)+"."+BaseTelescopeExFileName.split(".")[-1])
            elif SaveFileCommonStart:
                TelescopeNewExNames.append(SaveFileCommon+"_"+SourceIndices[ii]+"__"+"_".join(KeyValStrings)+"."+BaseTelescopeExFileName.split(".")[-1])
            else:
                TelescopeNewExNames.append(SaveFileCommon+"CombinationNo_"+SourceIndices[ii]+"__"+"_".join(KeyValStrings)+"."+BaseTelescopeExFileName.split(".")[-1])

        # Gether all information from each node
        pathlib.Path(SYNEX_PATH+"/TileTrackFiles/").mkdir(parents=True, exist_ok=True)
        CloningTrackFile=SYNEX_PATH+"/TileTrackFiles/ProgressTrackFile.txt"
        SourceExNames=[s.ExistentialFileName for s in sources]
        if MPI_rank>0:
            comm.send(SourceExNames, dest=0)
            comm.send(TelescopeNewExNames, dest=0)
            comm.send(CloningCombs, dest=0)
        else:
            SourceExNamesAll=copy.deepcopy(SourceExNames)
            TelescopeNewExNamesAll=copy.deepcopy(TelescopeNewExNames)
            CloningCombsAll=copy.deepcopy(CloningCombs)
            for core_ii in range(1,MPI_size):
                SourceExNamesAll+=comm.recv(source=core_ii)
                TelescopeNewExNamesAll+=comm.recv(source=core_ii)
                CloningCombsAll+=comm.recv(source=core_ii)

            # Save everything to txt file for tracking long lists of stuff on cluster -- Maybe this would be better as a dataframe save ? Easier to load again and handle using indexing to combinations ?
            with open(CloningTrackFile, 'w') as f:
                f.write(str(SaveInSubFile)+'\n')
                f.write(','.join(CloningKeys)+'\n')
                f.write(BaseTelescopeName+'\n')
                f.write(BaseTelescopeExFileName+'\n')
                for TelEx,SouEx,Comb in zip(TelescopeNewExNamesAll,SourceExNamesAll,CloningCombsAll):
                    f.write(TelEx+":"+SouEx+":"+",".join([str(el) if not isinstance(el,(np.ndarray,list)) else str(el[-1]) for el in Comb])+'\n') ### Need a way to save arrays better for when we switch to gaps etc. Maybe then we will switch to a '.dat' savefile instead and just pickle everything.
    else:
        ###
        #
        # This section is not up to date for cluster usage !!!!!!!!!!
        # Need to divide between cores, maybe find way to reduce memory usage
        # too like in the previous section where we pass only ExNames and create the
        # objects one by one?
        #
        # Need to check case where sources and detectors lists are not equal...
        #
        ###

        # Set base Ex name for tiling call later -- maybe we can do better with the way things are structured here. Seems like there are two calls here...
        # Maybe set the cloning stuff in a seperate function to clean things up here.
        TelescopeNewExNames=None

        # No cloning params - create detectors list from inputs
        if base_telescope_params==None and telescopes==None:
            telescopes=[SYDs.Athena(**{"verbose":verbose})]
        elif base_telescope_params!=None and telescopes==None:
            telescopes=[SYDs.Athena(**dict(base_telescope_params,**{"verbose":verbose}))]
        elif telescopes!=None:
            if base_telescope_params!=None and verbose: print("WARNING :: Case of input telescope list and base params ambiguous. Setting base params to None.")
            base_telescope_params=None
            if isinstance(telescopes,list):
                telescopes=[SYDs.Athena(**dict(tel,**{"verbose":verbose})) if isinstance(tel,dict) else setattr(tel,"verbose",verbose) for tel in telescopes]
            elif isinstance(telescopes,str):
                telescopes=glob.glob(telescopes) if "*" in telescopes else [SYDs.Athena(**dict(telescopes,**{"verbose":verbose}))]
            else:
                telescopes=[SYDs.Athena(**telescopes)] if isinstance(telescopes,dict) else [telescopes]

        # No cloning params - create sources list from inputs
        if sources==None:
            sources=[SYSs.SMBH_Merger(**{"verbose",verbose})]
        elif isinstance(sources,dict):
            sources=[SYSs.SMBH_Merger(**dict(sources,**{"verbose",verbose}))]
        elif isinstance(sources,list):
            sources=[SYSs.SMBH_Merger(**dict(s,**{"verbose",verbose})) if isinstance(s,dict) else setattr(s,"verbose",verbose) for s in sources]
        else:
            setattr(sources,"verbose",verbose)
            sources=[sources] # if MPI_rank==0 else None # if only one detector given

        # Check lengths of lists of objects and decide how many total tilings to do -- NB: both detectors and sources should now be lists but we never specified if they are the same length or not... Need to work on this ambiguous case...
        Nsources,Nteles = len(sources),len(telescopes)
        if Nsources>1 and Nteles>1:
            # List of sources and list of detectors given... Ensure lengths are the same.
            if Nsources!=Nteles: raise ValueError("inputting a list of sources with a list of telescopes with lengths not equal is ambiguous... Try implementing cloning function instead to ensure all combinations of sources and telescopes are accounted for.")
            Nvals=Nsources
        elif Nsources==1 and Nteles>1:
            # List of detectors given
            Nvals=Nteles
        elif Nsources>1 and Nteles==1:
            # List of sources given
            Nvals=Nsources
        else:
            # One of each; special case
            Nvals=1

        # Divide objects evenly between cores available
        if Nvals>1 and MPI_size>1:
            NValsPerCore=int(Nvals//MPI_size)
            CoreLenVals=[NValsPerCore+1 if ii<Nvals%MPI_size else NValsPerCore for ii in range(MPI_size)]
            CPU_ENDs=list(np.cumsum(CoreLenVals))
            CPU_STARTs=[0]+CPU_ENDs[:-1]
            telescopes = [telescopes[ii] for ii in range(CPU_STARTs[MPI_rank],CPU_ENDs[MPI_rank])]
            sources = [sources[ii] for ii in range(CPU_STARTs[MPI_rank],CPU_ENDs[MPI_rank])]

        # Sort some additional stuff
        if Nvals==1 and MPI_rank>0:
            ### Special case of single detector and single source ###
            # Not sure how to treat special case of many cores but one det and one source... Does gwemopt handle cores instrinsically? Need to verify this...
            # For now just have the master do stuff.
            telescopes=None
            sources=None
        else:
            # Fill missing EM data
            for s,tel in zip(sources,telescopes):
                if not hasattr(s,"EM_Flux_Data"): s.GenerateEMFlux(fstart22=1e-4,TYPE="const",**{})
                # Calculate source CTR data (detector dependent)
                s.GenerateCTR(tel.ARF_file_loc_name,gamma=1.7)

    #####
    #
    # At this point we should have a list of sources and detector Ex names (or detector objects if we passed H5file names at input).
    # Both sources and detector-related lists have equal length:
    #              * Each object in one list corresponds to the equivalently placed object in the second list *
    #
    # EXCEPTION :: input is single source and single detector WHILE using multiple cores.
    # Then detectors and sources both = None for MPI_rank>1 -- assuming there is scope within GWEMOPT to
    # handle multiple cores. To do this though we need to change some default go_params flags (Ncores etc). This needs to be done still.
    #
    #####

    # Tile stuff
    if telescopes!=None and sources!=None:
        # Output a check
        print("Beginning tiling for",len(telescopes),"telescopes, and", len(sources),"sources...")

        # Output location -- I think gwemopt_output folder is redundant now. Keep in case of relics.
        OutPutArch=SaveInSubFile.strip("/") if SaveInSubFile else None

        # Tile one by one using matching objects in each list
        T_tiling_start=time.time()
        for ii,(s,tel) in enumerate(zip(sources,telescopes)):
            print("-"*100)
            T_tiling_start_ii=time.time()
            if tel.telescope_source_coverage==None: go_params, map_struct, tile_structs, coverage_struct, tel = TileWithGwemopt(s,tel,OutPutArch,verbose)
            T_tiling_end_ii=time.time()
            print("Time to tile object no.",ii+1,"out of",len(telescopes),":",T_tiling_end_ii-T_tiling_start_ii,"s")
        T_tiling_end=time.time()
        TotTime=T_tiling_end-T_tiling_start
        TotTime=str(TotTime)+" s" if TotTime<3600. else str(TotTime/3600.)+" hr"
        print("-"*100)
        print("Total time for",len(telescopes),"objects:",TotTime)
        print("-"*100)

        # Return detectors output list if we are not on cluster
        if MPI_size==1: return telescopes
    elif TelescopeNewExNames!=None and sources!=None: # Skip special case where sources==None for one input object each on cluster...
        # Output check that cluster is provisioning correctly
        print(MPI_rank+1,"/",MPI_size,"with",len(TelescopeNewExNames),"/",Nvals,"telescopes to tile, and ",len(sources),"/",Nvals,"sources to tile.")

        # Output location -- I think gwemopt_output folder is redundant now. Keep in case of relics.
        OutPutArch=SaveInSubFile.strip("/") if SaveInSubFile else None

        # Initiate empty detectors output list if we are not on cluster
        if MPI_size==1: telescopes_out=[]

        # Loop over lists of sources and telescope save files
        for i,ExName in enumerate(TelescopeNewExNames):
            # Detector object vals -- load if file exists but include clone vals in case we changed something
            # (then we recalculate everything)
            Dictii=dict(base_telescope_params,
                      **{"ExistentialFileName":BaseTelescopeExFileName,
                      "NewExistentialFileName":ExName,
                      "verbose":verbose,
                      "cloning keys":CloningKeys,
                      "cloning values":CloningCombs[i],
                      "telescope":BaseTelescopeName+"_"+"_".join(CloningKeys)+"_"+str(i+1)},
                      **{CloningKeys[jj]:CloningCombs[i][jj] for jj in range(len(CloningKeys)) if CloningKeys[jj] not in ["Tcut"]})

            # Detector object
            telescope=SYTs.Athena(**Dictii) ## IF a param changed from existing detect file, Athena innit function forces detector_source_coverage = None.

            # Tile
            if MPI_rank==0:
                t_tile0=time.time()
                TelescopeCovered=True if telescope.telescope_source_coverage==None else False
                PrintMsg=telescope.telescope_source_coverage["source save file"] if telescope.telescope_source_coverage!=None else "No source coverage..."
                print("Pre-tiling master node check:",PrintMsg)
            if telescope.telescope_source_coverage==None: go_params, map_struct, tile_structs, coverage_struct, telescope_out = TileWithGwemopt(sources[i],telescope,OutPutArch,verbose)
            if MPI_rank==0:
                t_tile1=time.time()
                PrintMsg=telescope_out.telescope_source_coverage["source save file"] if telescope_out.telescope_source_coverage!=None else "No source coverage..."
                print("Post-tiling master node check:",PrintMsg)
                print("Time for telescope",i+1,"/",len(TelescopeNewExNames),"by master rank:",t_tile1-t_tile0,"s, with source coverage tileranges:",telescope_out.telescope_source_coverage["Source tile timeranges (isot)"],"\n")

            # Add source IDs
            telescope_out.telescope_source_coverage["source ID"] = SourceIndices[i]

            # Add to detectors output list if we not on cluster
            if MPI_size==1: telescopes_out.append(telescope_out)

        # Return detectors output list if we are not on cluster
        if MPI_size==1: return telescopes_out

def TileWithGwemopt(source,telescope,outDirExtension=None,verbose=True):
    """
    Helper function specifically to handle tiling through gwemopt.

    PARAMS
    ------
        - source :: SYNEX source object
            Unlike 'TileSkyArea()' master function, this must be a SINGLE class
            object at a time, since the GWEMOpt dictionaries are a mix of source and
            telescope parameters.
        - telescope :: SYNEX telescope object
            Unlike 'TileSkyArea()' master function, this must be a SINGLE class
            object at a time, since the GWEMOpt dictionaries are a mix of source and
            telescope parameters.
        - outDirExtension :: string [default None]
            Sub-architecture to save output telesope savefiles to.
        - verbose :: bool
            Option to output progress.

    OUTPUT
    ------
        - go_params :: dict
            GWEMOpt campatible dictionary
        - map_struct :: dict
            GWEMOpt campatible dictionary
        - tile_structs :: dict
            GWEMOpt campatible dictionary
        - coverage_struct :: dict
            GWEMOpt campatible dictionary
        - telescope :: SYNEX telescope class
            Output telescope class with extra attributes containing coverage information
            from GWEMOpt tiling and schedulers.
    """

    # Get the right dicts to use
    t_tile_0=time.time()
    go_params,map_struct=PrepareGwemoptDicts(source,telescope,outDirExtension)

    # Get segments
    import SYNEX.segments_athena as segs_a
    go_params = segs_a.get_telescope_segments(go_params)

    # Get tile_structs
    if go_params["tilesType"]=="MaxProb":
        moc_structs = gwemopt.moc.create_moc(go_params, map_struct=map_struct)
        tile_structs = gwemopt.tiles.moc(go_params,map_struct,moc_structs)
    elif go_params["tilesType"]=="moc":
        moc_structs = gwemopt.moc.create_moc(go_params, map_struct=map_struct)
        tile_structs = gwemopt.tiles.moc(go_params,map_struct,moc_structs,doSegments=False) # doSegments=False ?? Only for 'moc'... Otherwise it calls gwemopt.segments.get_segments_tiles
        for tel in tile_structs.keys():
            tile_structs[tel] = segs_a.get_segments_tiles(go_params, go_params["config"][tel], tile_structs[tel], verbose)
    elif go_params["tilesType"]=="greedy":
        tile_structs = gwemopt.tiles.greedy(go_params,map_struct)
        go_params["Ntiles"] = []
        for tel in go_params["telescopes"]:
            tile_structs[tel] = segs_a.get_segments_tiles(go_params, go_params["config"][tel], tile_structs[tel], verbose) # replace segs with our own
            go_params["config"][tel]["tesselation"] = np.empty((0,3))
            tiles_struct = tile_structs[tel]
            for index in tiles_struct.keys():
                ra, dec = tiles_struct[index]["ra"], tiles_struct[index]["dec"]
                go_params["config"][tel]["tesselation"] = np.append(go_params["config"][tel]["tesselation"],[[index,ra,dec]],axis=0)
            go_params["Ntiles"].append(len(tiles_struct.keys()))
    elif go_params["tilesType"]=="hierarchical":
        tile_structs = gwemopt.tiles.hierarchical(go_params,map_struct) # ,doSegments=False)
        go_params["Ntiles"] = []
        for tel in go_params["telescopes"]:
            tile_structs[tel] = segs_a.get_segments_tiles(go_params, go_params["config"][tel], tile_structs[tel], verbose) # replace segs with our own
            go_params["config"][tel]["tesselation"] = np.empty((0,3))
            tiles_struct = tile_structs[tel]
            for index in tiles_struct.keys():
                ra, dec = tiles_struct[index]["ra"], tiles_struct[index]["dec"]
                go_params["config"][tel]["tesselation"] = np.append(go_params["config"][tel]["tesselation"],[[index,ra,dec]],axis=0)
            go_params["Ntiles"].append(len(tiles_struct.keys()))
    elif go_params["tilesType"]=="ranked":
        moc_structs = gwemopt.rankedTilesGenerator.create_ranked(go_params,map_struct)
        tile_structs = gwemopt.tiles.moc(go_params,map_struct,moc_structs,doSegments=False)
        for tel in tile_structs.keys():
            tile_structs[tel] = segs_a.get_segments_tiles(go_params, go_params["config"][tel], tile_structs[tel], verbose)
    elif go_params["tilesType"]=="galaxy":
        # Really not sure how this works... Use this method with care.
        map_struct, catalog_struct = gwemopt.catalog.get_catalog(go_params, map_struct)
        tile_structs = gwemopt.tiles.galaxy(go_params,map_struct,catalog_struct)
        for tel in go_params["telescopes"]:
            # tile_structs[telescope] = segs_a.get_segments_tiles(go_params, go_params["config"][telescope], tile_structs[telescope], verbose)
            go_params["config"][tel]["tesselation"] = np.empty((0,3))
            tiles_struct = tile_structs[tel]
            for index in tiles_struct.keys():
                ra, dec = tiles_struct[index]["ra"], tiles_struct[index]["dec"]
                go_params["config"][tel]["tesselation"] = np.append(go_params["config"][tel]["tesselation"],[[index,ra,dec]],axis=0)

    # Allocate times to tiles
    tile_structs, coverage_struct = gwemopt.coverage.timeallocation(go_params, map_struct, tile_structs)

    # Get info for run
    telescope = GetCoverageInfo(go_params, map_struct, tile_structs, coverage_struct, telescope, source=source, verbose=False)

    # Add info about source to coverage summary
    SourceInfo={
                "source save file":source.ExistentialFileName,
                "source gps_timetomerger_max":source.gps_timetomerger_max,
                "source sky_map":source.sky_map,
                "source fisher area":source.FisherSkyArea,
                "source post area":source.PostSkyArea,
                "source JsonFile":source.JsonFile,
                "source H5File":source.H5File,
                "source Is3D":source.do3D,
                "source true loc":[source.true_ra,source.true_dec]
                }
    if source.do3D:
        SourceInfo["source true loc"]+=[source.true_distance]
    telescope.telescope_source_coverage.update(SourceInfo)

    # Save detector
    telescope.ExistentialCrisis()

    return go_params, map_struct, tile_structs, coverage_struct, telescope

def PrepareGwemoptDicts(source,telescope,outDirExtension=None):
    """
    Helper function to prepare GWEMOPT dictionary from a source object and a
    detector objct.

    PARAMS
    ------
        - Source :: SYNEX source object
            Unlike 'TileSkyArea()' master function, this must be a SINGLE class
            object at a time, since the GWEMOpt dictionaries are a mix of source and
            telescope parameters.
        - Detector :: SYNEX detector object
            Unlike 'TileSkyArea()' master function, this must be a SINGLE class
            object at a time, since the GWEMOpt dictionaries are a mix of source and
            telescope parameters.
        - outDirExtension :: string [default None]
            Sub-architecture to save output telesope savefiles to.
    OUTPUT
    ------
        - go_params :: dict
            A GWEMOpt - compatible dictionary containing global parameters for the
            tiling simulation (that is a mix of source and telescope specific parameters...)
        - map_struct :: dict
            A GWEMOpt - compatible dictionary containing parameters describing a
            single telescope.
    """
    # Prepare go_params from detector first
    go_params = copy.deepcopy(telescope.telescope_go_params)
    telescope_names = [telescope.telescope_config_struct["telescope"]] # They want this to be a list
    go_params["telescopes"]=telescope_names
    go_params["exposuretimes"]=np.array([telescope.telescope_config_struct["exposuretime"]]) # This line is necessary and not really sure why...
    go_params["config"]={telescope.telescope_config_struct["telescope"] : telescope.telescope_config_struct}

    # Ideal pixelation for detector FoV has pixel area <= FOV/4 (idk this for sure, I just made it up...)
    nside_arr=np.array([2**i for i in range(11,3,-1)])
    area_arr=hp.nside2pixarea(nside_arr,degrees=True)
    nside=nside_arr[np.searchsorted(area_arr, telescope.telescope_config_struct["FOV"]/4., side='left')-1]

    # Update missing params in go_params contained in source
    go_params["nside"]=nside
    go_params["do3D"]=source.do3D
    go_params["true_ra"]=source.true_ra
    go_params["true_dec"]=source.true_dec
    go_params["true_distance"]=source.true_distance

    ################################### THIS BLOCK NEEDS CHECKING ###########################################################################################
    # Check times match up an correct if they don't
    telescope_times_gps = [telescope.telescope_config_struct["gps_science_start"], telescope.telescope_config_struct["gps_science_start"] + telescope.telescope_config_struct["mission_duration"]*365.25]
    source_times_gps = [source.gps_timetomerger_max + source.timetomerger_max*365.25 + source.DeltatL_cut/(60.*60.*24.), source.gps_timetomerger_max + source.timetomerger_max*365.25]
    if source_times_gps[0]>telescope_times_gps[0]:
        go_params["gpstime"]=source_times_gps[0]
    else:
        go_params["gpstime"]=telescope.telescope_config_struct["gps_science_start"]

    # Constrain to mission time
    if go_params["gpstime"]+go_params["Tobs"][-1] > telescope_times_gps[1]:
        go_params["Tobs"][-1] = telescope_times_gps[1]-go_params["gpstime"]

    # Constain to source timeline ???
    if go_params["gpstime"]+go_params["Tobs"][-1] > source_times_gps[1]:
        go_params["Tobs"][-1] = source_times_gps[1]-go_params["gpstime"]
    #########################################################################################################################################################

    # Update the gwemopt output dir to comply with lisabeta h5 and json filenames
    # This will overwrite existing files with the same event name... Maybe later include checks?
    if not go_params["tilesType"] in ["hierarchical","greedy"]:
        ### DONT add subfolder structure to 'outputDir' if sampler based tilesType is requested:
        ### The pathlength becomes too long for pymultinest so filenames get truncated, leading to 'file not found' errors...
        ### "greedy" acutally uses emcee -- does this have the same problem?
        go_params["outputDir"]+="/"+".".join(source.JsonFile.split("/")[-1].split(".")[:-1])
        if not outDirExtension==None:
            go_params["outputDir"]+="_"+outDirExtension
    if source.PermissionToWrite: pathlib.Path(go_params["outputDir"]).mkdir(parents=True, exist_ok=True)
    if os.path.isfile(go_params["outputDir"]):
        print("WARNING: outputDir",go_params["outputDir"].split("/SYNEX/")[-1],"already exists... Files within will be overwriten.")

    # Now get map struct
    map_struct=source.map_struct

    # Run through the GWEMOPT checker for uniformity between go_params and sky_map (rescale nside etc)
    # We can also apply rotations etc here
    map_struct=gou.read_skymap(go_params,source.do3D,map_struct=map_struct)

    return go_params,map_struct

def GetCoverageInfo(go_params, map_struct, tile_structs, coverage_struct, telescope, source, verbose=True):
    """
    Helper function to extract coverage information, including a count for tiles
    that cover source and photons captured origintaing at the source.

    PARAMS
    ------
        - go_params :: dict
            A GWEMOpt - compatible dictionary containing global parameters for the
            tiling simulation (that is a mix of source and telescope specific parameters...)
        - map_struct :: dict
            A GWEMOpt - compatible dictionary containing parameters describing a
            single telescope.
        - tile_structs
            A GWEMOpt - compatible dictionary containing parameters describing a
            single telescope's tiling of a single source's skymap.
        - coverage_struct
            A GWEMOpt - compatible dictionary containing information of a single
            telescope's coerage of a single source's skymap (obtained after GWEMOpt is run).
            this function effectively adds information to this in a new attribute
            in an output telescope class.
        - telescope :: SYNEX telescope class
            Telescope class output from a GWEMOpt run.
        - source :: SYNEX source class [default None]
            Source class needed for CTR data to caluclate captured photons and tiles
            that contained the source pixel location.
        - verbose :: bool [default True]
            Option to output progress.

    OUTPUT
    ------
        - telescope :: SYNEX telescope class
            SYNEX telescope class with updated information on source coverage. This
            is where we add information on e.g. number of times the telescope
            catches the source, how many photons are collected per exposure and in
            total, as well as tile information like time ranges (isot AND mjd) and
            days till first exposure from start of tiling procedure.

    TO DO:
    ------
        Include check that source exists here - if not: create it, or read it's
        position from json, or can we just use the detector_go_params ? Is it updated?
        Or can the position and dist be included in the detector_source_coverage dictionary?
    """

    # Get source location
    phiSSB = -source.true_lamdaSSB
    thetaSSB = np.pi/2-source.true_betaSSB
    source_pix = hp.ang2pix(go_params["nside"],thetaSSB,phiSSB,lonlat=False) # Can we set lonlat to True and input labmda beta directly? Have not checked this...

    # Extract some data structures
    cov_data=coverage_struct["data"] # Data is N*9 with columns ra,dec,mjd (exposure start),mag,exposureTime,field,tile_probs,airmass,program_id
    cov_ipix=coverage_struct["ipix"] # list because each sub-list is of variable length
    map_probs=map_struct["prob"]

    # Tiling info
    telescope_name=telescope.telescope_config_struct["telescope"]
    tile_struct=tile_structs[telescope_name]
    tile_keys_list=list(tile_struct.keys())
    pix_record=np.unique(np.concatenate([tile_struct[tile_id]["ipix"] for tile_id in tile_keys_list if len(tile_struct[tile_id]["ipix"])>0],axis=0))
    probs1=np.array([tile_struct[tile_id]["prob"] for tile_id in tile_keys_list])
    probs2=np.array([np.sum(map_probs[tile_struct[tile_id]["ipix"]]) for tile_id in tile_keys_list])
    TotProb1=np.sum(probs1) # Tile struct probabilities
    TotProb2=np.sum(probs2) # Mp struct probs summed using ipix in tile struct
    TotProb3=np.sum(map_probs[pix_record]) # Unique probabilities
    telescope.telescope_tile_struct=tile_struct
    if verbose:
        print("\n\n")
        print("#"*20+"  "+telescope_name+"  "+"#"*20+"\n")
        print("Total prob (tile_struct, map_struct, map_struct drop dups):",
              TotProb1, TotProb2, TotProb3)
        print("Params/Config checks:", telescope_name,
              go_params["config"][telescope]["tot_obs_time"],
              go_params["config"][telescope_name]["exposuretime"],
              go_params["config"][telescope_name]["tot_obs_time"]/go_params["config"][telescope_name]["exposuretime"],
              "Tot tiles avail:", len(tile_struct.keys()))
        print("\n")

    # Coverage info
    cov_ipix_array_unique=np.unique([pix for cov_ipix_el in cov_ipix for pix in cov_ipix_el])
    prob1=np.sum(cov_data[:,6])
    prob2=np.sum([map_probs[pix] for tile_pix in cov_ipix for pix in tile_pix]) # Probs according to coverage struct ipix
    prob3=np.sum([map_probs[pix] for pix in cov_ipix_array_unique]) # Same as prob2 but only unique pixels
    ex_Times = cov_data[:,4]
    ex_Times_unique = np.unique(ex_Times)
    telescope.telescope_coverage_struct=coverage_struct
    telescope.telescope_coverage_struct.update({"tile exposure times":ex_Times,"unique tile exposure times":ex_Times_unique,
                                              "normalised tile exposure times":[T_e/telescope.telescope_config_struct["exposuretime"] for T_e in ex_Times]})

    # Print summary if asked for
    if verbose:
        print(len(ex_Times),"/",len(cov_data[:,4]), "coverage tiles")
        print("Detector coverage prob (cov_struct data, map_struct, map_struct dup drop):",
              prob1, prob2, prob3)
        print(len(ex_Times_unique),"/",len(cov_data[:,4]),
              "unique coverage tiles")
        print(len(ex_Times_unique),"unique exposure times with (min, mean, max):",
              np.min(ex_Times_unique),
              np.mean(ex_Times_unique),
              np.max(ex_Times_unique))
        print("\n")

    ## Source tile info -- was it covered?
    UniqueSourceCoverageCount=0
    cov_source_tile=[]
    cov_source_pix=[]
    SourceTileTimeRanges=[]
    SourceTileTimeRanges_days=[]
    for tile_ii,tile_pix in enumerate(cov_ipix):
        TotExpTime_day=cov_data[tile_ii,2]-cov_data[0,2]
        # print(tile_ii, source_pix, tile_pix, cov_data[tile_ii,4], detector.detector_config_struct["exposuretime"], cov_data[tile_ii,4]/detector.detector_config_struct["exposuretime"])
        if source_pix in tile_pix:
            cov_source_tile += [tile_ii] # [tile_ii for tile_ii,tile_pix in enumerate(cov_ipix) if source_pix in tile_pix]
            cov_source_pix += list(tile_pix)
            SourceTileTimeRanges_days.append([TotExpTime_day,TotExpTime_day+cov_data[tile_ii,4]/86400.])
            SourceTileTimeRange=[Time(cov_data[tile_ii,2], format='mjd', scale='utc').isot,
                               Time(cov_data[tile_ii,2]+cov_data[tile_ii,4]/86400., format='mjd', scale='utc').isot]
            SourceTileTimeRanges.append(SourceTileTimeRange)
            if all([pix not in cov_source_pix for pix in tile_pix]):
                UniqueSourceCoverageCount+=1
    SourceTile_prob1=[cov_data[i,6] for i in cov_source_tile]
    SourceTile_prob2=[map_probs[pix] for pix in cov_source_pix]
    SourceTile_accum_prob1=np.sum(SourceTile_prob1)
    SourceTile_accum_prob2=np.sum(SourceTile_prob2)

    telescope.telescope_source_coverage={"telescope":telescope_name,"Source Tile Prob (by cov tiles)":SourceTile_prob1,
                        "Source Tile Prob (by cov pixels)":SourceTile_prob2,"tile ID containing source pixel":cov_source_tile,
                        "tile pixels containing source pixel":cov_source_pix, "No. source coverages":len(SourceTile_prob1),
                        "No. unique source coverages":UniqueSourceCoverageCount,
                        "Source tile timeranges (isot)": SourceTileTimeRanges,
                        "Source tile timeranges (days)": SourceTileTimeRanges_days}

    # Print summary if asked for
    if verbose:
        print("--- Source Coverage ---")
        print("\n")
        print("Source covered by",telescope_name, len(cov_source_tile),"times.\n")
        print(" "*10+"cov_struct p data:")
        print(SourceTile_prob1,"\n")
        print(" "*10+"map_struct pixel p's:")
        print(SourceTile_prob2,"\n")
        print("Source tile accumulated p:",SourceTile_accum_prob1,SourceTile_accum_prob2)

    # if len(cov_source_tile)>0: print("Coverage stuct exposures within range:",Time(cov_data[0,2],format='mjd', scale='utc').isot,Time(cov_data[-1,2],format='mjd', scale='utc').isot)

    # Now add cuts to coverage info depending on photon flux
    if source!=None and len(cov_source_tile)>0:
        # Make sure source has relevant data and create it if not
        # if not hasattr(source,"EM_Flux_Data"): source.GenerateEMFlux(fstart22=1e-4,TYPE="const",**{})
        # if not hasattr(source,"CTR_Data"): source.GenerateCTR(detector.ARF_file_loc_name,gamma=1.7) # Should be include gamma as a source param? Would we ever want to change this at run time?
        source.GenerateEMFlux(fstart22=1e-4,TYPE="const",**{})
        source.GenerateCTR(telescope.ARF_file_loc_name,gamma=1.7)

        # Calculate the exposuretimes for each tile that covers source
        SourceTileExpTimes=[cov_data[tile_ii,4] for tile_ii in cov_source_tile]
        SourceTileStartTimes=[86400.*(cov_data[tile_ii,2]-Time(go_params["gpstime"], format='gps', scale='utc').mjd) for tile_ii in cov_source_tile]

        # Get CTR data out from source and cut to Tobs -- times here are seconds to merger (<0)
        CTRs=source.CTR_Data["CTR"]
        CTR_times=list(86400.*(-np.array(source.EM_Flux_Data["xray_time"])/86400. - Time(source.gpstime_merger_max, format='gps', scale='utc').mjd + Time(go_params["gpstime"], format='gps', scale='utc').mjd)) # Just in case these are different -- i.e. if we have a run with multiple sources within Tobs
        CTRs=[CTR for CTR,t in zip(CTRs,CTR_times) if t>=(go_params["Tobs"][0]*86400.) and t<=(go_params["Tobs"][-1]*86400.)]
        CTR_times=[t for t in CTR_times if t>=(go_params["Tobs"][0]*86400.) and t<=(go_params["Tobs"][-1]*86400.)] # because times are seconds TO MERGER (-ve)

        # Integrate CTR for each tile covering source to check no. of source photons received
        TileListTimes=[[ts+i*dur/49 for i in range(50)] for ts,dur in zip(SourceTileStartTimes,SourceTileExpTimes)]
        SourceTilePhotonCounts=[np.trapz(np.interp(TileTs,CTR_times,CTRs),TileTs) for TileTs in TileListTimes]

        telescope.telescope_source_coverage.update({"Source tile exposuretimes (s)":SourceTileExpTimes,
                                                  "Source photon counts":SourceTilePhotonCounts,
                                                  "Source tile start times (s)":SourceTileStartTimes, ### Maybe change this and exp times to be a list of pairs like 'timeranges' but in gps days from start time? Then you can easily make the required histpgrams of interest.
                                                  "Start time (mjd)":Time(go_params["gpstime"], format='gps', scale='utc').mjd})
    else:
        telescope.telescope_source_coverage.update({"Source tile exposuretimes (s)":[],
                                                  "Source photon counts":[0],
                                                  "Source tile start times (s)":[],
                                                  "Start time (mjd)":Time(go_params["gpstime"], format='gps', scale='utc').mjd})

    # Return values
    return telescope

def WriteSkymapToFile(map_struct, SkyMapFileName,
                      go_params=None, PermissionToWrite=True):
    """
    Helper funciton to save a skymap to fits file.

    PARAMS
    ------
        - map_struct
            A GWEMOpt - compatible dictionary containing parameters describing a
            single source's sky map attribute (source.map_struct).
        - SkyMapFileName
            file name including full path and extension to save the skymap to. This
            can be easily source.sky_map, which is the checked filename constructed
            at source creation usually based on the source.JsonFile file name.
        - go_params :: dict [default None]
            A GWEMOpt - compatible dictionary containing parameters describing
            global run parameters of a GWEMOpt run. This the GWEMOpt ready dictionary,
            output from the function 'PrepareGwemoptDicts()'.
        - PermissionToWrite :: bool [default True]
            Flag to allow the function to write to file all the relevant information.
            This is automatically set to False when multiprocessing is sensed (via
            a successful import of mpi4py) to avoid overly using processor memory
            when on the cluster for example.

    OUTPUT
    ------
        1. IF 'go_params' is specified, then:
            - go_params
                Updated version of the input, GWEMOpt-compatible dictionary describing
                global function for the GWEMOpt run.
            - map_struct
                The input, GWEMOpt-compatible dictionary describing sky map.
        2. IF 'go_params' is NOT specified, then:
            - SkyMapFileName
                Complete filename of where the sky map was saved.
            - map_struct
                The input, GWEMOpt-compatible dictionary describing sky map.

    Note: filename is JUST the name without the path- the path is calculated in situ.
    """
    # from astropy.table import Table, Column
    # import os.path

    # get path and check filename...
    # TODO: add this to a helper function to check paths of data directories
    # are congruent....
    # Can also do this for tile functions in the same function but different
    # to lisabeta function, but this function reference lisabeta one
    # to harmonize the architectures..

    # ???
    if len(SkyMapFileName.split("Skymap_files"))==1 and len(SkyMapFileName.split("SYNEX"))>1:
        raise ValueError("Sorry but for continuity across functions please"
                          "direct skymap directories to be within './SYNEX/Skymap_files/...'")
    elif len(SkyMapFileName.split("Skymap_files"))>1 and len(SkyMapFileName.split("SYNEX"))==1:
        print("Directing given filename to './SYNEX/Skymap_files/...'")
        SkyMapFileName = SYNEX_PATH + "/Skymap_files" + SkyMapFileName.split("Skymap_files")[-1]
    elif len(SkyMapFileName.split("Skymap_files"))==1 and len(SkyMapFileName.split("SYNEX"))==1:
        print("WARNING: Adding ./SYNEX/Skymap_files/ to filename provided...")
        SkyMapFileName = SYNEX_PATH + "/Skymap_files/" + SkyMapFileName

    if '.fits' not in SkyMapFileName:
        SkyMapFileName = SkyMapFileName + '.fits'
    # Check if directory exists and create if it doesnt
    SkyMapFilePath = "/".join(SkyMapFileName.split("/")[:-1])
    pathlib.Path(SkyMapFilePath).mkdir(parents=True, exist_ok=True)

    # Write to fits file
    if "distmu" in map_struct:
        data_to_save = np.vstack((map_struct["prob"],
                                  map_struct["distmu"],
                                  map_struct["distsigma"],
                                  map_struct["distnorm"]))
        # np.vstack((map_struct["prob"], map_struct["cumprob"],
        # map_struct["ipix_keep"], map_struct["pixarea"],
        # map_struct["pixarea_deg2"],map_struct["distmu"],
        # map_struct["distsigma"],map_struct["distnorm"]))
        if PermissionToWrite:
            hp.write_map(SkyMapFileName,
                         data_to_save,
                         overwrite=True,
                         column_names=["prob", "distmu",
                                       "distsigma", "distnorm"])
        if go_params!=None:
            go_params["do3D"]=True
    elif PermissionToWrite:
        data_to_save = map_struct["prob"]
        # np.vstack((map_struct["prob"], map_struct["cumprob"],
        #            map_struct["ipix_keep"], map_struct["pixarea"],
        #            map_struct["pixarea_deg2"]))
        hp.write_map(SkyMapFileName, data_to_save, overwrite=True)
        # , column_names=["prob"])

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

# plotName = plotName[:-5]+"_prob.pdf"
# Check extension is there - move this to a general function please
# if plotName[-4:]!=".pdf":
#     plotName = ".".join(plotName.split(".")[:-1]) + ".pdf"
# Check if directory tree exists
# PlotPath="/".join(plotName.split("/")[:-1])
# pathlib.Path(PlotPath).mkdir(parents=True, exist_ok=True)

##### Need to organize these better and rename #####
# SaveFig=False, plotName=None):
def PlotSkyMapData(source, save_path=None, plot_name=None, extension=".png"):
    """
    Plotting tool to plot a source object's skymap.
    Adapted from 'gwemopt.plotting.skymap.py'
    """
    unit = 'Gravitational-Wave probability'
    cbar = False
    title = r"Lisabeta Localisation"
    xsize = 1500 # Image size

    if np.percentile(source.map_struct["prob"], 99) > 0:
        hp.mollview(source.map_struct["prob"], unit=unit, title=None,
                    cbar=cbar, min=np.percentile(source.map_struct["prob"], 1),
                    max=np.percentile(source.map_struct["prob"],99), cmap=cmap,
                    xsize=xsize)
    else:
        hp.mollview(source.map_struct["prob"], unit=unit, title=None,
                    cbar=cbar, min=np.percentile(source.map_struct["prob"], 1),
                    cmap=cmap, xsize=xsize)

    # Projplot has funky conventions... This was discovered playing around with
    # coordinates and seeing how they moved.
    # They inverted phi and this needs to be in radians...
    phi = -source.true_lamdaSSB
    theta = np.pi/2 - source.true_betaSSB
    hp.projplot(theta, phi, lonlat=False, coord=None, marker='*',
                markersize=10, markeredgewidth=1, c='black', linestyle='None',
                label='True Location')
    fig = plt.gcf()
    ax = plt.gca()
    plt.rcParams.update({'font.size':16})

    ax.set_title(title, fontsize=20)
    ax.legend(fontsize=18)

    add_edges()
    if save_path:
        # Default save name
        extensions = ["png", "pdf"]
        if plot_name:
            # add extension if necessary
            if plot_name.split(".")[-1] not in extensions:
                plot_name += extension
            plot_name = os.path.join(save_path, plot_name)
        else:
            plot_name = os.path.join(save_path, "skymap_proba_density.pdf")
        # Save
        fig.savefig(plot_name, dpi=200)
        plt.close(fig)
    else:
        plt.show()

    # Extra plotting for 3D skymaps
    if "distmu" in source.map_struct:
        fin = np.copy(source.map_struct["distmu"])
        fin[~np.isfinite(fin)] = np.nan
        hp.mollview(source.map_struct["distmu"], unit='Distance [Mpc]',
                    min=np.nanpercentile(fin,10), max=np.nanpercentile(fin,90))
        add_edges()
        if save_path:
            plot_name_dist = plot_name.split(".")[0] + '_distance.pdf'
            plt.savefig(plot_name_dist, dpi=200)
        else:
            plt.show()
        plt.close('all')

    if "distmed" in source.map_struct:
        fin = np.copy(source.map_struct["distmed"])
        fin[~np.isfinite(fin)] = np.nan
        hp.mollview(source.map_struct["distmed"], unit='Distance [Mpc]',
                    min=np.nanpercentile(fin,10), max=np.nanpercentile(fin,90))
        add_edges()
        if save_path:
            plot_name_median_dist = plot_name.split(".")[0] + '_dist_median.pdf'
            plt.savefig(plot_name_median_dist, dpi=200)
        else:
            plt.show()
        plt.close('all')

def PlotTilesArea_old(TileFileName,n_tiles=10):
    """
    ################# This function is no longer frequently used #################

    Function to plot a sample of tiles from a saved dictionary of tiles.

    It is no longer used because it was created to plot tile dictionary after tiling
    with the most basic possible tiling algorithm, which has since been replaced with
    GWEMOpt.

    ################# This function is no longer frequently used #################
    """
    # Load Tile dictionary
    with open(TileFileName, 'rb') as f:
        TileDict = pickle.load(f)
    from matplotlib.patches import Rectangle

    # Set the default color cycle
    TileColor=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
               "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
               "#bcbd22", "#17becf"]
    n_tile_colours = 10

    # Get background data for liklihoods
    fig, InjParam_InjVals, SampleModes, X, lambda_bins, Y, beta_bins,Z = PlotInferenceLambdaBeta(TileDict["LISA Data File"],
                                                                                                 bins=50, SkyProjection=False,
                                                                                                 SaveFig=False, return_data=True)
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

def PlotLikeRatioFoMFromJsonWithInferenceInlays(FoMJsonFileAndPath, BF_lim=20.,
                                                SaveFig=False,
                                                InlayType="histogram"):
    """
    ################# This function is no longer frequently used #################

    Just as with all util functions relating to FoM analyses (e.g. "RunFoMOverRange()")
    this is no longer used since we want to focus on randomized drawings rather than
    focussed lines of systems drawn from parameter space.

    NB: This function has hard coded paths in it and will likely be removed wince we have moved on from needing it

    ################# This function is no longer frequently used #################
    """

    with open(FoMJsonFileAndPath) as f:
        FoMOverRange = json.load(f)

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
            # print("Posterior sample length: " + str(nsamples) +
            #       ", number of infered parameters: " + str(ndim))

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
            # print("Posterior sample length: " + str(nsamples) +
            #       ", number of infered parameters: " + str(ndim))
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
    ################# This function is no longer frequently used #################

    Just as with all util functions relating to FoM analyses (e.g. "RunFoMOverRange()")
    this is no longer used since we want to focus on randomized drawings rather than
    focussed lines of systems drawn from parameter space.

    NB: This function has hard coded paths in it and will likely be removed wince we have moved on from needing it

    ################# This function is no longer frequently used #################
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
    ################# This function is no longer frequently used #################

    Plotting function to output histogram plots of lambda and beta from a lisabeta
    inference run.

    This funciton is not well maintained, since we focus now on healpy oriented
    plotting to standardize conventions with GWEMOpt and observationalists.

    ################# This function is no longer frequently used #################
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
    # print("Posterior sample length: " + str(nsamples) + ", number of infered parameters: " + str(ndim))
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

def set_nice_label_corner_plot(inferred_param_list):
    """
    Take the list of inferred parameters and returns a list with the
    corresponding label in latex format for a more aestetically pleasant
    corner plot.
    """
    nice_label_dic = {
        # Total *redshifted* mass M=m1+m2, solar masses
        "M": r'M$_{tot}$',
        # Chirp Mass
        "Mchirp": r"$\mathcal{M}_{c}$",
        # Mass ratio q=m1/m2
        "q": '$q$',
        # Dimensionless spin component 1 along orbital momentum
        "chi1": r"$\chi_{1}$",
        # Dimensionless spin component 2 along orbital momentum
        "chi2": r"$\chi_{2}$",
        # Dimensionless mass-weighted spin positive component
        "chim": r"$\chi_{-}$",
        # Dimensionless mass-weighted spin negative component
        "chip": r"$\chi_{+}$",
        # Luminosity distance, Mpc
        "dist": r"$D_L$",
        # Inclination, observer's colatitude in source-frame
        "inc": r"$\iota$",
        # Phase, observer's longitude in source-frame
        "phi": r"$\phi$",
        # Longitude in the sky
        "lambda": r"$\lambda$",
        # Latitude in the sky
        "beta": r"$\beta$",
        # Polarization angle
        "psi": r"$\Psi$",
        }

    nice_label_list = []
    for inferred_param_list_item in inferred_param_list:
        item_label = nice_label_dic[inferred_param_list_item]
        nice_label_list.append(item_label)

    return nice_label_list


def PlotInferenceData(FileName, extension='png'):
    """
    Plotting function to output corner plot of inference data after a lisabeta
    inference run.

    PARAMS
    ------
        - FileName :: string
            File name of lisabeta data. This only have to include sub-architecture
            of JsonFile or H5File location, and does not nedd to include file extensions.
        - SaveFig :: bool
            Flag for the function to save the plot. The plot is saved as png to
            "./SYNEX/Plots/lisabeta/" with a name corresponding to the H5File filename.

    NOTE:
    The savefile does NOT include sub-architecture in the savefile location. This makes
    it NOT congruent with JsonFile and H5File.
    """
    import corner

    # Check filenames
    JsonFileLocAndName,H5FileLocAndName = CompleteLisabetaDataAndJsonFileNames(FileName)

    # Get some useful info from Json first
    with open(JsonFileLocAndName) as f: json_data = json.load(f)
    labels = json_data["prior_params"]["infer_params"] # should be a list of inferred parameters
    ndim = len(labels)

    # Unpack the data
    [infer_params, inj_param_vals, static_params, meta_data] = read_h5py_file(H5FileLocAndName)
    if not inj_param_vals: # Raw files at first didn't record this, so make sure it's there...
        # Get the inj values from the processed data file instead
        DataFileLocAndName_NotRaw = DataFileLocAndName[:-7] + ".h5"
        [_,inj_param_vals,_,_] = read_h5py_file(DataFileLocAndName_NotRaw)
    if np.size(infer_params[labels[0]][0])>1:
        nsamples = len(infer_params[labels[0]][0])
    else:
        nsamples = len(infer_params[labels[0]])
    # print("Posterior sample length: " + str(nsamples) + ", number of infered parameters: " + str(ndim))
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
    for key in labels:
        if json_data["run_params"]["sample_Lframe"]:
            InjParam_InjVals.append(inj_param_vals["source_params_Lframe"][key][0])
        else:
            InjParam_InjVals.append(inj_param_vals["source_params_SSBframe"][key][0])

    labels = set_nice_label_corner_plot(labels)
    # Corner plot of posteriors
    figure = corner.corner(data, labels=labels,
                           quantiles=[0.16, 0.5, 0.84],
                           show_titles=True,
                           label_kwargs=dict(fontsize=16),
                           title_kwargs=dict(fontsize=16),)
    # Set the size a the Figure
    figure.set_size_inches(20., 20.)

    for ax in figure.get_axes():
        ax.tick_params(axis='both', labelsize=16)

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
            # ax.xaxis.set_tick_params(labelsize=20)
            # ax.yaxis.set_tick_params(labelsize=20)
            # ax.set_ylabel(fontsize=20)

    # save the figure
    # Put in folder for all lisabeta-related plots
    SaveFile = (os.path.basename(FileName).split('.')[0] +
                '_corner.{}'.format(extension))
    SaveRep = os.path.dirname(FileName) + '/../Plots'
    if not os.path.exists(SaveRep):
        os.makedirs(SaveRep)
    SaveFile = os.path.join(SaveRep, SaveFile)
    # pathlib.Path(SYNEX_PATH + "/Plots/lisabeta/").mkdir(parents=True, exist_ok=True)
    # SaveFile = SYNEX_PATH + "/Plots/lisabeta/" + SaveFile.split("/")[-1]
    plt.savefig(SaveFile)
    plt.close(figure)
    return 0

def PlotInferenceLambdaBeta(FileName, bins=50, SkyProjection=False,
                            SaveFig=False, return_data=False, extension='png'):
    """
    Plotting function to show just the localization parameters lambda and beta
    after a lisabeta inference run. This is effectively a two-dimensional grid
    binning of the sky parameter posteriors (because I have trust issues with
    2D histograms I found online...).

    PARAMS
    ------
        - FileName :: string
            filename of either JsonFile for H5File with any sub-architecture under
            'inference_param_files' or 'inference_data' respectively and with or without
            file extensions.
        - bins :: int [default 50]
            Number of bins in one direction. The total number of bins plotted is bins^2.
            Note that the histogram bins is refined according to the range of values in
            both lambda and beta posteriors to ensure that we get the right amount
            of resolution depending on the spread of posterior points across the sky.
        - SkyProjection :: bool [default False]
            Flag to plot in the "mollweide" projection. Default is cartesian.
        - SaveFig :: bool [default False]
            Flag to save the plot. The plot is saved as png to
            "./SYNEX/Plots/lisabeta/" with a name corresponding to the H5File filename.

            NOTE:
            The savefile does NOT include sub-architecture in the savefile location. This makes
            it NOT congruent with JsonFile and H5File.

        - return_data :: bool [default False]
            Flag to return the histogramed data and figure object.

    OUTPUT
    ------
        1. If we specify 'return_data=True', then:
            - fig
                Figure handle.
            - InjParam_InjVals
                Values of injected parameters.
            - SampleModes
                Modes of posterior samples.
            - X :: np.meshgrid
                Meshgrid of lambda and beta bins.
            - lambda_bins
                Bins used for lambda posteriors.
            - Y :: np.meshgrid
                Meshgrid of lambda and beta bins.
            - beta_bins
                Bins used for beta posteriors.
            - Z :: np.meshgrid
                Bin populations normalized by bin area and total population.
                NOTE: bin normalization needs to be checked before trusted...
        2. Otherwise nothing is returned and the plot is rendered.
    """
    # Update some general properties
    # GeneralPlotFormatting()

    # Check filenames
    JsonFileLocAndName, H5FileLocAndName = CompleteLisabetaDataAndJsonFileNames(FileName)

    # Get some useful info from Json first
    with open(JsonFileLocAndName) as f: Jsondata = json.load(f)
    labels = ["lambda","beta"]
    # Jsondata["prior_params"]["infer_params"]
    # should be a list of inferred parameters
    ndim = len(labels)
    print(labels)

    # Unpack the data
    [infer_params, inj_param_vals,
     static_params, meta_data] = read_h5py_file(H5FileLocAndName)
    if np.size(infer_params[labels[0]][0])>1:
        nsamples = len(infer_params[labels[0]][0])
    else:
        nsamples = len(infer_params[labels[0]])
    # print("Posterior sample length: " + str(nsamples) +
    #       ", number of infered parameters: " + str(ndim))

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
            mode = mode[1]
            # Does now take the first on the list,
            # but if there are several we need to work out what to do there...
        SampleModes.append(mode)
    data = np.transpose(np.array(data)) # should have shape [nsamples, ndim]

    # Get injected values
    InjParam_InjVals = []
    for key in labels:
        if Jsondata["run_params"]["sample_Lframe"]:
            InjParam_InjVals.append(inj_param_vals["source_params_Lframe"][key][0])
        else:
            InjParam_InjVals.append(inj_param_vals["source_params_SSBframe"][key][0])

    # Scatter plot of labmda and beta posterior chains, marginalized over all other params
    # fig,ax = plt.figure()
    if SkyProjection: plt.subplot(111, projection="mollweide") # To match gwemopt projections

    # Get 2D hitogram data
    levels = [1.-0.997, 1.-0.95, 1.-0.9, 1.-0.68]
    # Plot contours for 1 sigma, 90% confidence, 2 sigma, and 3 sigma
    levels_labels=[str(int((1.-l)*1000)/10) for l in levels]

    # Manually do 2D histogram because I don't trust the ones I found online
    # data to hist = [data[:,0], data[:,1]] = lambda, beta. Beta is y data since it is declination.
    hist2D_pops = np.empty([bins,bins])
    areas = np.empty([bins,bins])
    bin_max = max(max(data[:,0]), SampleModes[0]+np.pi/10000.)
    bin_min = min(min(data[:,0]), SampleModes[0]-np.pi/10000.)
    if bin_max>np.pi:
        bin_max=np.pi
    if bin_min<-np.pi:
        bin_max=-np.pi
    if isinstance(bin_min,np.ndarray):
        bin_min = bin_min[0]
    if isinstance(bin_max,np.ndarray):
        bin_max = bin_max[0]
    lambda_bins = np.linspace(bin_min, bin_max, bins+1) # np.linspace(np.min(data[:,5]), np.max(data[:,5]), bins+1)
    bin_max = max(max(data[:,1]), SampleModes[1]+np.pi/20000.)
    bin_min = min(min(data[:,1]), SampleModes[1]-np.pi/20000.)
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
        list_len = list(range(len(data[:,0])))
        if xii == 0:
            values_ii = [valii for valii in list_len if (lambda_bins[xii]<=data[valii,0]<=lambda_bins[xii+1])]
        else:
            values_ii = [valii for valii in list_len if (lambda_bins[xii]<=data[valii,0]<lambda_bins[xii+1])]
        beta_pops, beta_bins = np.histogram(data[values_ii,1], bins=beta_bins) # beta_bins
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
    for l,s in zip(levels, levels_labels):
        labels_dict[l]=s

    # Plot contour
    contour = plt.contour(X, Y, Z, levels)
    plt.clabel(contour, colors ='k', fmt=labels_dict, fontsize=20)
    fig, ax = plt.gcf(), plt.gca()

    # Set the size a the Figure
    fig.set_size_inches(15., 15.)

    # Add injected and mode vertical and horizontal lines
    if not SkyProjection:
        linewidth = 4
        ax.axhline(InjParam_InjVals[1], color="r", linestyle=":",
                   linewidth=linewidth)
        ax.axhline(SampleModes[1], color="b", linestyle=":",
                   linewidth=linewidth)
        ax.axvline(InjParam_InjVals[0], color="r", linestyle=":",
                   linewidth=linewidth)
        ax.axvline(SampleModes[0], color="b", linestyle=":",
                   linewidth=linewidth)

    # Labels
    if not SkyProjection:
        fontsize=25
        if Jsondata["run_params"]["sample_Lframe"]:
            plt.xlabel(r"$\lambda_{L}$", fontsize=fontsize) # Lambda
            plt.ylabel(r"$\beta_{L}$", fontsize=fontsize) # beta
            frame = "LISA"
        else:
            plt.xlabel(r"$\lambda_{SSB}$", fontsize=fontsize) # Lambda
            plt.ylabel(r"$\beta_{SSB}$", fontsize=fontsize) # beta
            frame = "SSB"

        # Add points at injected and mode values
        ax.plot(InjParam_InjVals[0], InjParam_InjVals[1], "*",color="r", ms=25,
                label="Injected Position in {} frame".format(frame))
        ax.plot(SampleModes[0], SampleModes[1], "o", color="b",  ms=18,
                label="Recovered Position in {} frame".format(frame))

        # Increase size of the tick in axes so that I can read
        label_size = 20
        ax.xaxis.set_tick_params(labelsize=label_size)
        ax.yaxis.set_tick_params(labelsize=label_size)
        ax.set_title(r'2D Localisation with Lisabeta', fontsize=30)
        # legend informing if we are Lframe or not
        ax.legend(fontsize=fontsize)

    # show now or return?
    if not return_data:
        # plt.grid()
        # plt.show()

        # save the figure if asked
        if SaveFig:
            # Put in folder for all lisabeta-related plots
            # SaveFile = H5FileLocAndName[:-3]+'.png'
            # pathlib.Path(SYNEX_PATH + "/Plots/lisabeta/").mkdir(parents=True, exist_ok=True)
            # SaveFile = SYNEX_PATH + "/Plots/lisabeta/" + SaveFile.split("/")[-1]
            # plt.savefig(SaveFile)
            # save the figure
            # Put in folder for all lisabeta-related plots
            SaveFile = (os.path.basename(FileName).split('.')[0] +
                        '_skymap.{}'.format(extension))
            SaveRep = os.path.dirname(FileName) + '/../Plots'
            if not os.path.exists(SaveRep):
                os.makedirs(SaveRep)
            SaveFile = os.path.join(SaveRep, SaveFile)
            # pathlib.Path(SYNEX_PATH + "/Plots/lisabeta/").mkdir(parents=True, exist_ok=True)
            # SaveFile = SYNEX_PATH + "/Plots/lisabeta/" + SaveFile.split("/")[-1]
            plt.savefig(SaveFile)
            plt.close(fig)
    else:
        return fig, InjParam_InjVals, SampleModes, X, lambda_bins, Y, beta_bins, Z

def PlotLikeRatioFoMFromJson(FoMJsonFileAndPath, BF_lim=20., SaveFig=False):
    """
    ################# This function is no longer frequently used #################

    Function to eventually replace the base plot functions in other plot util.
    Return an axis and figure handle that has the base colour plot for the
    lambda and beta sky mode replection using log(bayes factor) calculations only.

    NB: This function has hard coded paths in it and will likely be removed wince we have moved on from needing it

    This function is left over from sky mode studies that are no longer being pursued.

    ################# This function is no longer frequently used #################
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
    Function to plot the Earth, Sun and Moon as seen by a telescope in orbit
    given a set of parameters describing the telescope trajectory, mostly to
    check we have done things right for orbital calculations in segments_athena.py.

    Note: the plot contains circular patches representing the Earth, Sun and Moon
    whose size on the plot are accurate to the angular extent as viewed by the
    telescope along it's orbit.

    PARAMS
    ------
        - config_struct :: dict
            GWEMOPT 'config_struct' dictionary containing relevant information
            describing a telescope object. Based on the telescope.telescope_config_struct
            base dictionary, but passed through the 'PrepareGwemoptDicts()' function with
            a source, and then run through GWEOPT via the SYNEX master tiling functions.
            This ensures that the extra SYNEX module segments_athena.py has processed
            the dictionaries and produced orbit information.
        - MollProj bool [default False]
            Option to plot in a "MollProj" projection instead of cartesian projection.

            NOTE:
            MollProj=True will use "mollweide" projection but this takes radians as input
            and will give some glitchy things for patch radii due to their small
            angular extents. If you set to cartesian (mollweide=False) patches will have
            radii in degrees and not have the same problem.
        - SaveFig :: bool [default False]
            Option to save the plot as pdf in "./SYNEX/Plots/orbits/" under a long name
            corresponding to the orbit parameters.
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
        if not MollProj and i==0:
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
    Animation of the Earth, Sun and Moon in 3D cartesian space as seen by a
    telescope in orbit.
    Mostly just to be cool.

    To change from 'mp4' to 'swf' format (eg for inclusion in latex slides), use command line:
    "ffmpeg -i in.mp4 out.swf"

    Note: the animation contains circular patches representing the Earth, Sun and
    Moon whose time varying size on the plot are accurate to the time varying angular
    extent as viewed by the telescope along it's orbit.

    PARAMS
    ------
        - config_struct :: dict
            GWEMOPT 'config_struct' dictionary containing relevant information
            describing a telescope object. Based on the telescope.telescope_config_struct
            base dictionary, but passed through the 'PrepareGwemoptDicts()' function with
            a source, and then run through GWEOPT via the SYNEX master tiling functions.
            This ensures that the extra SYNEX module segments_athena.py has processed
            the dictionaries and produced orbit information.
        - include_sun :: bool [default False]
            Option to include the Sun in the animation. If True, then the animation
            of the Earth and Moon are virtually undetectible given their distance
            from the telescope relative to the distance between the Suna nd the telescope.

            NOTE:
            If we include the Sun then we normalize distance to 1=1AU. If we do not
            include Sun, then distances are normalized to 1=1 mean radius of telescope
            orbit.
        - SaveAnim :: bool [default False]
            Option to save the animation as mp4 in "./SYNEX/Plots/OrbitAnimations/"
            under a long name corresponding to the orbit parameters.
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
    Animation of the Earth, Sun and Moon in 2D space as seen by a telescope in orbit.
    Mostly just to be cool.

    SkyProj=True will use "mollweide" projection but this takes radians as input
    and will give give some glitchy things for patch radii due to their small
    angular extents. If you set to cartesian (mollweide=False) patches will have
    radii in degrees and not have the same problem.

    To change from 'mp4' to 'swf' format (eg for inclusion in latex slides), use command line:
    "ffmpeg -i in.mp4 out.swf"

    PARAMS
    ------
        - config_struct :: dict
            GWEMOPT 'config_struct' dictionary containing relevant information
            describing a telescope object. Based on the telescope.telescope_config_struct
            base dictionary, but passed through the 'PrepareGwemoptDicts()' function with
            a source, and then run through GWEOPT via the SYNEX master tiling functions.
            This ensures that the extra SYNEX module segments_athena.py has processed
            the dictionaries and produced orbit information.
        - MollProj bool [default False]
            Option to plot in a "MollProj" projection instead of cartesian projection.

            NOTE:
            MollProj=True will use "mollweide" projection but this takes radians as input
            and will give some glitchy things for patch radii due to their small
            angular extents. If you set to cartesian (mollweide=False) patches will have
            radii in degrees and not have the same problem.
        - SaveAnim :: bool [default False]
            Option to save the animation as mp4 in "./SYNEX/Plots/OrbitAnimations/"
            under a long name corresponding to the orbit parameters.
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
    """
    Simple loglog plot of EM flux as a function of time.

    PARAMS
    ------
        - source :: SYNEX source class
            Must have had 'source.GenerateEMFlux()' run so that it contains EM data
            attributes to plot.
        - SaveFig :: bool [default False]
            Option to save the plot as a pdf at "./SYNEX/Plots/" with filename
            constructed from "source.ExistentialFileName".
    """
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
    """
    Simple loglog plot of CTR as a function of time.

    PARAMS
    ------
        - source :: SYNEX source class
            Must have had 'source.GenerateCTR()' run with a specific telescope ARF
            file so that it contains CRT data attributes to plot.
        - SaveFig :: bool [default False]
            Option to save the plot as a pdf at "./SYNEX/Plots/" with filename
            constructed from "source.ExistentialFileName".
    """
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

def PlotPhotonAccumulation(telescopes, SaveFig=False, SaveFileName=None):
    """
    ################# This function is no longer frequently used #################

    Plot accumulated photons from source only after tiling is run through GWEMOPT.

    Function is less well maintained after implementing dataframe options instead.
    this function avoids this and grabs information directly from the telescope objects
    which can be cumbersome if we revisit a stuy later and want to plot something
    but we have to load a long list of telescopes instead of just a CSV-saved dataframe.

    ################# This function is no longer frequently used #################
    """
    if not isinstance(telescopes,list):
        telescopes=[telescopes]

    # Check if it's a list of lists
    if any([isinstance(telescope,list) for telescope in telescopes]):
        telescopes=[[telescopes_el] if not isinstance(telescopes[0],list) else telescopes_el for telescopes_el in telescopes] ### Should the else part be after the for blah in blahs part?
    else:
        telescopes=[telescopes]

    # Create list of detectors if dicts given or mix of dicts and detector objects given
    # detectors=[detector if not isinstance(detector,dict) else SYDs.Athena(**detector) for detector in detectors]
    print("Detectors shape check:",type(telescopes),len(telescopes),[len(telescopes_el) for telescope_el in telescopes])

    # Create list of sources in case we have a mix of stuff like DeltatL_cut
    sources=[[GetSourceFromLisabetaData(telescope.telescope_source_coverage["source H5File"],**{"ExistentialFileName":telescope.telescope_source_coverage["source save file"],"verbose":telescope.verbose}) for telescope in telescopes_el] for telescopes_el in telescopes]
    SourcePreMergerCuts=[[source.DeltatL_cut for source in source_el] for source_el in sources]
    for ii in range(len(telescopes)):
        if len(np.unique(SourcePreMergerCuts[ii]))>1:
            for jj in range(len(detectors[ii])):
                telescopes[ii][jj].telescope_config_struct.update({"cloning key":"DeltatL_cut","cloning value":-sources[ii][jj].DeltatL_cut/86400.})
    t_mergers=[[-source.DeltatL_cut/86400. for source in source_el] for source_el in sources]

    # Find the max height so all bars are even
    for ii in range(len(telescopes)):
        height=max([sum(telescope.telescope_source_coverage["Source photon counts"]) for telescope in telescopes[ii]])

        for telescope,t_merger in zip(telescopes[ii],t_mergers[ii]):
            T0_mjd=telescope.telescope_source_coverage["Start time (mjd)"]
            Xs=[ExT0s/86400. for ExT0s in telescope.telescope_source_coverage["Source tile start times (s)"]]
            Xs_Widths=[ExTs/86400. for ExTs in telescope.telescope_source_coverage["Source tile exposuretimes (s)"]]
            Ys=list(np.cumsum(telescope.telescope_source_coverage["Source photon counts"]))

            # Plot accumulated photons if source was captured
            if len(telescope.telescope_source_coverage["Source tile start times (s)"])>0:
                if telescope.telescope_config_struct["cloning value"]!="DeltatL_cut" and telescope.telescope_config_struct["cloning value"]!=None:
                    label=telescope.telescope_config_struct["cloning key"]+r"="+str(int(telescope.telescope_config_struct["cloning value"]//1)) if isinstance(telescope.telescope_config_struct["cloning value"],(float,int)) else telescope.telescope_config_struct["cloning value"]
                elif telescope.telescope_config_struct["cloning value"]=="DeltatL_cut":
                    label=r"System " + str(ii) + r", Source time to merger="+str(int(t_merger//1))+r"\,d"
                else:
                    label=telescope.telescope_config_struct["telescope"]
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

def PlotSourcePhotons(telescopes, labels=None, BoxPlot=True, SaveFig=False, SaveFileName=None):
    """
    ################# This function is no longer frequently used #################

    Plot total accumulated photons for each list of detectors (in case we
    have several versions of the same cloned parameter but using eg different
    gwemopt flags...)

    Function is less well maintained after implementing dataframe options instead.
    this function avoids this and grabs information directly from the telescope objects
    which can be cumbersome if we revisit a stuy later and want to plot something
    but we have to load a long list of telescopes instead of just a CSV-saved dataframe.

    ################# This function is no longer frequently used #################
    """
    # Check if it's a list
    if not isinstance(telescopes,list):
        telescopes=[telescope]

    # Check if it's a list of lists
    if any([isinstance(telescope,list) for telescope in telescopes]):
        telescopes=[[telescopes_el] if not isinstance(telescopes[0],list) else telescopes_el for telescopes_el in telescopes]
    else:
        telescopes=[telescopes]

    # Create list of sources in case we have a mix of stuff like DeltatL_cut
    sources=[[GetSourceFromLisabetaData(telescope.telescope_source_coverage["source H5File"],**{"ExistentialFileName":telescope.telescope_source_coverage["source save file"],"verbose":detector.verbose}) for telescope in telescopes_el] for telescopes_el in telescopes]



    ############# Did we change something? ------ First handle two special cases then gneral case (TO DO) #############
    # Add this to master tiling function.

    # Did we change DeltatL_cut?
    SourcePreMergerCuts=[[source.DeltatL_cut for source in source_el] for source_el in sources]
    for ii in range(len(detectors)):
        if len(np.unique(SourcePreMergerCuts[ii]))>1:
            for jj in range(len(detectors[ii])):
                telescopes[ii][jj].telescope_config_struct.update({"cloning keys":"DeltatL_cut","cloning values":-sources[ii][jj].DeltatL_cut/86400.})

    # Did we change Tobs?
    telescopeTobs=[[telescope.telescope_go_params["Tobs"] for telescope in telescopes_el] for telescopes_el in telescopes]
    for ii in range(len(telescopes)):
        if not np.array_equal(np.unique(telescopeTobs[ii]),telescopeTobs[ii][0]):
            for jj in range(len(telescopes[ii])):
                telescopes[ii][jj].telescope_config_struct.update({"cloning keys":"Tobs","cloning values":telescopeTobs[ii][jj][1]})

    ###################################################################################################################



    # Handle labels now if not given
    if not isinstance(labels,(list,int,str,float)): labels=[None]*len(telescopes)
    # elif isinstance(labels,(int,str,float)): labels=[labels]*len(detectors)   ##### Not sure what the default behaviour should be here...

    # Get the figure and axes for each list within detectors list
    if BoxPlot:
        Xs = [telescope.telescope_config_struct["cloning values"] for telescopes_el in telescopes for telescope in telescopes_el] # This will store values cloned for sources too...
        Ys = [sum(telescope.telescope_source_coverage["Source photon counts"]) for telescopes_el in telescopes for telescope in telescopes_el]
        UniqueXs = pd.unique(Xs)
        UniqueXs.sort()
        clone_key = telescopes[0][0].telescope_config_struct["cloning keys"]
        data_points = pd.DataFrame(data={clone_key:Xs, "photon_count":Ys})
        xdata_stats = data_points.groupby(by=clone_key).describe() # Default is only include numeric data
        print(xdata_stats)
        BoxPlots = data_points.boxplot(vert=True, by=clone_key, column=["photon_count"], return_type="dict")
        for k in BoxPlots["photon_count"]:
            for Line2D in BoxPlots["photon_count"][k]:
                xtmp = Line2D.get_xdata()
                if len(xtmp)>0:
                    ii=int(round(np.mean(xtmp)))-1
                    Line2D.set_xdata(xtmp - round(np.mean(xtmp)) + UniqueXs[ii])
    else: ### Can we just use dataframes to draw the plots below? Do we need the extra helper function?
        # Init global figure
        fig=plt.figure()

        Lines2D=[(PlotSourcePhotons_SingleDetList(telescopes_el,sources=sources,fig=fig,label=label,BoxPlot=BoxPlot)) for telescopes_el,label in zip(telescopes,labels)]

        # Plot the mean
        clone_key = telescopes[0][0].detector_config_struct["cloning keys"]
        data_points = pd.DataFrame(data={clone_key:[x for Line2D in Lines2D for x in Line2D[2]],
                                   "photon_count":[y for Line2D in Lines2D for y in Line2D[3]]})
        xdata_stats = data_points.groupby(by=clone_key).describe() # Default is only include numeric data
        print(xdata_stats)
        mean_xs=pd.unique(data_points[clone_key])
        mean_xs.sort()
        mean_ys=data_points.groupby(by=clone_key).mean().to_numpy().flatten()
        mean_stds=data_points.groupby(by=clone_key).std().to_numpy(na_value=0.).flatten() ### Make suree to treat case where we have only one data line to avoid Nan errorbars
        mean_ns=data_points[clone_key].value_counts().sort_index(inplace=False).to_numpy()
        errors=mean_stds/np.sqrt(mean_ns)
        print(" --- Means -- ")
        for i in range(len(mean_xs)):
            print(mean_xs[i],":",round(mean_ys[i]),r"\pm",round(errors[i]))
        print(" ------------ ")
        plt.errorbar(mean_xs,mean_ys,yerr=errors, marker='+', color='black', linestyle='', markersize=2, figure=fig, label='mean')

        # Formatting labels
        if labels!=None: plt.legend(fontsize="x-small",ncol=data_points[clone_key].nunique()//10+1) # Wrap columns when we get more than 10 points

    # Sort special case of premerger time cut
    if telescopes[0][0].telescope_config_struct["cloning keys"]=="DeltatL_cut":
        ax=plt.gca()
        ax.invert_xaxis()
        X_Datas_unique = pd.unique(data_points[clone_key])
        plt.xticks(X_Datas_unique, labels=['-{0:.1f}'.format(d) for d in X_Datas_unique]) # , rotation='vertical')
        plt.xlabel(r"Pre-merger cut [days to merger]",fontsize="small")

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

def PlotSourcePhotons_SingleDetList(telescopes,sources=None,fig=None,label=None,BoxPlot=True):
    """
    ################# This function is no longer frequently used #################

    Helper function for "PlotSourcePhotons()" function that takes a list of telescopes
    instead of the structured list (a list of lists of telescopes depending on what we tiled)
    that PlotSourcePhotons() takes.

    Function is less well maintained after implementing dataframe options instead.
    this function avoids this and grabs information directly from the telescope objects
    which can be cumbersome if we revisit a stuy later and want to plot something
    but we have to load a long list of telescopes instead of just a CSV-saved dataframe.

    ################# This function is no longer frequently used #################
    """
    if not isinstance(telescopes,list):
        telescopes=[telescopes]

    # Create list of detectors if dicts given or mix of dicts and detector objects given
    telescopes=[telescope if not isinstance(telescope,dict) else SYDs.Athena(**telescope) for telescope in telescopes]

    # Get cloned values and photon counts
    Xs = [telescope.telescope_config_struct["cloning value"] for telescope in telescopes] # This will store values cloned for sources too...
    Ys = [sum(telescope.telescope_source_coverage["Source photon counts"]) for telescope in telescopes]

    # Init figure if it isn't already
    if fig==None: fig=plt.figure()
    if label==None: label=''

    # Do we need logarithmic axes?
    ax = plt.gca()
    if min(Xs)!=0. and np.log10(max(Xs)/min(Xs))>1.5:
        ax.set_xscale('log')
    if min(Ys)!=0. and np.log10(max(Ys)/min(Ys))>2.5:
        ax.set_yscale('log')

    # Label axes
    plt.xlabel(telescopes[0].telescope_config_struct["cloning key"],fontsize="small")
    plt.ylabel(r"Accumulated Photons",fontsize="small")

    # Plot
    if BoxPlot:
        plt.boxplot(Ys,meanline=True)
        plt.show()
    else:
        ErrBarCont=plt.errorbar(Xs,Ys,yerr=np.sqrt(Ys), marker='+', linestyle='', markersize=2, figure=fig, label=label)
        # Return
        return ErrBarCont, ax, Xs, Ys

def PlotCoverage(source,telescope):
    """
    Plot all coverage information for a single telescope and source. This uses
    GWEMOpt built-in functions, and therefore requires conversion from SYNEX classes
    to GWEMOpt compatible dictionaries. We also hard set the options 'catalog_struct'
    and 'plot_sun_moon' to None and Fale respectively since these are not useful
    in SYNEX applications. The former is for catalogue information pertaining to
    lower mass binaries outside the SMBH applications of SYNEX, and the latter
    plots the sun and moon location under the assumption that the telescope is
    Earth bound. In order to get an equivalent plot for SYNEX purposes, we must
    access, for example, the orbit dictionary in the telescope class caluclated
    from orbit parameters stored in the telescope.telescope_config_struct attribute.

    PARAMS
    ------
        - source :: SYNEX source class
        - telescope :: SYNEX telescope class
    """
    # Get necessary gwemopt-compatible dictionaries
    go_params,map_struct=PrepareGwemoptDicts(source,telescope)
    coverage_struct=telescope.telescope_coverage_struct

    # Plot using gwemopt built-in stuff
    gwemopt.plotting.coverage(go_params, map_struct, coverage_struct, catalog_struct=None,plot_sun_moon=False)








#              #####################################
#              #                                   #
#              #   SYNEX -- Post tiling analysis   #
#              #                                   #
#              #####################################

def GetDataFrame(telescopes=None, SaveFile=None):
    """
    Function to get or create a dataframe from a list or array of
    telescopes with tiling information.

    Dataframes are preferred for post-tiling analyses since loading a h5-saved
    table is much more memory and time efficient rather than loading the
    list of telescoped each time we wish to start a new study.

    PARAMS
    ------
        - telescopes :: list/array of SYNEX telescope classes or savefile names.
            Must have gwemopt output information (e.g. detector.detector_source_coverage!=None).
        - SaveFile :: string.
            Filename with sub-architecture in which dataframe should be saved or loaded from.
            If the file already exists and a list of detectors given also, then the
            savefile will be overwritten. Otherwise the file will be loaded or
            created given input parameters.

    NOTE: The save location for a new dataframe is always "/./SYNEX/DataFrames/"
    to save any problems with architecture etc.

    TO DO:
    ------
    Need to add check that the location of DataFrames exists before we try to save,
    otherwise new users will get errors that the file location does not exist.

    """
    # Work out what we want to do based on inputs.
    if telescopes==None and SaveFile==None:
        raise ValueError("Need to give either a list of detectors, a savefile to load, or a list of detectors and savefile to save the dataframe to.")
    elif telescopes==None and SaveFile:
        store = pd.HDFStore(SaveFile)
        return store['data']
    else:
        return CreateDataFrameFromDetectorList(telescopes=telescopes, SaveFile=SaveFile)

def CreateDataFrameFromDetectorList(telescopes, SaveFile=None):
    """
    Create a pandas dataframe from a list of telescopes after running through the
    SYNEX master tiling function. It takes out all useful information from the
    telescope classes as well as the useful information from corresponding sources
    whose skymap was tiled over, and puts it all into a table. This is much more
    efficient when tiling is finished and analysis of results begins, since instead
    of loading all classes again from scratch, we can load instead a simple
    h5-saved table. This drastically reduces computation time and memory.

    PARAMS
    ------
        - Detectors :: list/array of SYNEX detectors or telescope savefile names.
            Must have gwemopt output information (e.g. detector.detector_source_coverage!=None).
        - SaveFile :: string.
            Savefile and directory sub-architecture in which dataframe should be saved.
            The base folder is forced to be "/./SYNEX/DataFrames/" in order to be
            congruent with other save locations like JsonFile. The default when this is None
            is to save at "./SYNEX/DataFrames/data.h5". Extension to file is not needed
            since ".h5" is manually added when dataframe is saved.

    OUTPUT
    ------
        - data :: pandas dataframe
            Information from all detectors on coverage formatted into a dataframe.
            this maked visualization and statistical manipulation much faster and compact
            to study correlations between coverage information and input parameters.

            For example, we can tile the same source many times with a number of
            "exposureTime", and then we can use the dataframe to see a positive
            correlation between all the telescope.exposureTime and
            telescope.telescope_source_coverage["Source photon counts"] without having
            to manually load all the telescopes each time we initiate the study notebook.
            Instead, dataframes are much faster to load, and can be saved easily
            as csv files for much smaller disk space.

    TO DO:
    ------
        - Include error bounds from lisabeta inferece? Might be interesting to
          compare error bounds for a static parameter when randomizing others ?
        - Optimise... Creating each source from posterior data files every time
          is ridiculously laboursome and I am sure there is a better way to do this.
    """
    ###
    # Variables of interest:
    # source -- {Mtot, q, chi1, chi2, lambda, beta, dist, DeltatL_cut}
    # detect -- {Tobs, Texp}
    # FoMs -- {n_photons, DaysToSourceExp}
    ###

    # Did we get a list?
    if not isinstance(telescopes,list): telescopes=[telescopes]

    # Did we get a set of savefiles?
    if isinstance(telescopes[0],str): telescopes=[SYDs.Athena(**{"ExistentialFileName":ExName}) for ExName in telescopes]

    # List of all params of interest... This should not be hardcoded moving forward as a priori we don't know what is interesting...
    SourceTestParamKeys = ["M","q","chi1","chi2","lambda","beta", "dist", "DeltatL_cut"]
    telescopeTestParamKeys = ["Tobs", "exposuretime","n_photons_per_tile","n_photons","DaysToSourceExp"]
    AllTestParamKeys = SourceTestParamKeys+telescopeTestParamKeys
    AllTestParams = {key:[] for key in AllTestParamKeys}

    # Load params
    for tel in telescopes:
        with open(tel.telescope_source_coverage["source JsonFile"]) as f:
            input_params = json.load(f)
        for key in AllTestParamKeys:
            if key in input_params['source_params'] or key in ["M","q"]:
                if key=="M":
                    AllTestParams[key].append(round(input_params['source_params']["m1"]+input_params['source_params']["m2"])) ## Otherwise we get fractional differences and this can screw up methods in dataframe
                elif key=="q":
                    AllTestParams[key].append(round(input_params['source_params']["m1"])/round(input_params['source_params']["m2"]))
                else:
                    AllTestParams[key].append(input_params['source_params'][key])
            elif key in input_params['waveform_params']:
                AllTestParams[key].append(input_params['waveform_params'][key])
            elif key in tel.telescope_go_params:
                if key=="Tobs":
                    AllTestParams[key].append(tel.telescope_go_params[key][1]) ### But what if we have gaps... How do we treat array of pairs? Seperate out photons per obs period maybe? Effectively 'exploding' the dataframe? Need to include this calc then inside 'GetCoverageInfo'
                else:
                    AllTestParams[key].append(tel.telescope_go_params[key])
            elif key in tel.telescope_config_struct:
                AllTestParams[key].append(tel.telescope_config_struct[key])

        # Add system ID if asked for
        # AllTestParams["source ID"] = d.detector_source_coverage["source ID"]

        # Add total source photons captured
        AllTestParams["n_photons_per_tile"].append(tel.telescope_source_coverage["Source photon counts"])
        AllTestParams["n_photons"].append(round(sum(tel.telescope_source_coverage["Source photon counts"])))

        # Add first day after detector start of science at which source is exposed -- append 'NAN' if source never exposed
        if len(tel.telescope_source_coverage["Source tile timeranges (days)"])>0:
            AllTestParams["DaysToSourceExp"].append(tel.telescope_source_coverage["Source tile timeranges (days)"])
        else:
            AllTestParams["DaysToSourceExp"].append([[np.nan]])

        # Add sky areas
        ##################################################### TMP CODE #####################################################
        if "source fisher area" in tel.telescope_source_coverage:
            AllTestParams["Fisher Area (sq deg)"]=tel.telescope_source_coverage["source fisher area"]
            AllTestParams["Posterior Area (sq deg)"]=tel.telescope_source_coverage["source post area"]
        else:
            source_kwargs={"ExistentialFile":tel.telescope_source_coverage["source save file"],
                           "verbose":False}
            source=GetSourceFromLisabetaData(tel.telescope_source_coverage["source JsonFile"],**source_kwargs) # hopefully ust gotta do this once...
            AllTestParams["Fisher Area (sq deg)"]=source.FisherSkyArea
            AllTestParams["Posterior Area (sq deg)"]=source.PostSkyArea
            # Re-write source and detector savefiles so everything is right
            tel.ExistentialCrisis()
            source.ExistentialCrisis()
        ##################################################### TMP CODE #####################################################

    # Create DataFrame
    data = pd.DataFrame(data=AllTestParams)

    # Savefile name and location (force placement in ../SYNEX/DataFrames/ directory)
    SaveFileDir=SYNEX_PATH + '/DataFrames/'
    if SaveFile:
        if SaveFile[-1]=="/": SaveFile+="data" # IF we are only given a directory
        SaveFile=".".join(SaveFile.split('.')[:-1]) + ".h5" if "." in SaveFile else SaveFile+".h5" # Make sure file type is right
        SaveFile=SaveFile.split("/SYNEX/")[-1] # Take out path to SYNEX in case this is changed between systems
        SaveFile=SaveFileDir+SaveFile.split("/DataFrames/")[-1] # Make sure we save in the right place
    else:
        SaveFile=SaveFileDir+'data.h5'

    # Check path exists
    pathlib.Path("/".join(SaveFile.split("/")[:-1])).mkdir(parents=True, exist_ok=True)

    # Save dataframe
    print("Saving dataframe to:",SaveFile)
    store = pd.HDFStore(SaveFile)
    store['data'] = data  # save it

    # Return dataframe
    return data






####################################################################
#                                                                  #
#    PSO -- NO LONGER USED IN SYNEX BUT MIGHT BE USEFUL LATER ?    #
#                                                                  #
####################################################################

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
