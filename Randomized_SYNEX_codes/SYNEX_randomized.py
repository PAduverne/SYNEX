import numpy as np
import glob,copy,sys

from SYNEX import SYNEX_Detectors as SYDs
from SYNEX import SYNEX_Sources as SYSs
from SYNEX import SYNEX_Utils as SYU
from SYNEX.SYNEX_Utils import SYNEX_PATH

try:
    from mpi4py import MPI
except ModuleNotFoundError:
    MPI = None

# Using MPI ?
use_mpi=False
is_master = True
mapper = map
if MPI is not None:
    MPI_size = MPI.COMM_WORLD.Get_size()
    MPI_rank = MPI.COMM_WORLD.Get_rank()
    comm_global = MPI.COMM_WORLD
    use_mpi = (MPI_size > 1)
    if use_mpi:
        pool = ptemcee.mpi_pool.MPIPool(debug=False)
        is_master = pool.is_master()
        mapper = pool.map
    else:
        is_master = True
        mapper = map




# Functions to randomize spins, angles, Mass rat.
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

# Check for existing json files (optional)
CHECK_FOR_JSONS=True
if is_master and CHECK_FOR_JSONS:
    JsonFiles = glob.glob(SYNEX_PATH+"/inference_param_files/Randomized_*.json")
elif is_master:
    JsonFiles=[]
else:
    JsonFiles=None

if is_master and len(JsonFiles)==0:
    # Draw random values
    n = 10
    rand_spins = draw_random_spins(size=n)
    rand_angles = draw_random_angles(size=n)
    rand_massratios = draw_random_massratio(size=n)

    # Detector object options
    LISA_base_kwargs = {"TDI":'TDIAET',"verbose":False}

    # Initialize detector object
    LISA = SYDs.LISA(**LISA_base_kwargs)

    # Set Tcut if looking at pre-merger - 'None' to ignore
    OneHour = 60.*60.
    OneDay = 24.*OneHour
    OneWeek = 7.*OneWeek
    T_obs_end_to_mergers = [-30.*OneDay,-3.*OneWeek,-2.OneWeek,-OneWeek,
                        -3.*OneDay,-OneDay,-10.*OneHour,-5.*OneHour,None]
    T_obs_labels = ["1mon","3wk","2wk","1wk","3d","1d",
                        "10hr","5hr","0cut"]

    # Source object options
    # NB Dec = beta and Ra = lambda
    Merger_base_kwargs = {"q": 1.1, "M": 8e6, "z": 1.3, "chi1": 0.9,
            "chi2": 0.95, "beta" : -3.*np.pi/8., "lambda" : np.pi/3.,
            "inc": np.pi/10., "psi": 0.4,  "approximant" : 'IMRPhenomHM',
            "Lframe":True, "DeltatL_cut":None,
            "verbose":False,
            "ExistentialFileName":
                SYNEX_PATH+"/Saved_Source_Dicts/TestSystem_9d_base.dat",
            "NewExistentialFileName":
                SYNEX_PATH+"/Saved_Source_Dicts/TestSystem_9d_tmp.dat"
    }

    # Params to be infered
    inference_params = {
      "infer_params": ["chi1", "chi2", "dist", "inc", "phi", "lambda", "beta",
                        "psi"],
      "params_range": [[-1., 1.],[-1., 1.],[5e3, 2e5],[0.,np.pi],
                        [-np.pi, np.pi],[-np.pi, np.pi], [-np.pi/2.,np.pi/2.],
                        [0.,np.pi]],
      "prior_type": ["uniform", "uniform", "uniform", "sin", "uniform",
                        "uniform", "cos", "uniform"],
      "wrap_params": None
    }

    # MCMC options
    RunTimekwargs = {"print_info": True,
                    "n_walkers":  96, # must be >= 2 * inference cube dimension
                    "n_iter": 8000,
                    "burn_in": 5000, # Throw away at least 2/3
                    "autocor_method": "autocor_new",
                    "thin_samples": True,
                    "TDI": "TDIAET",
                    "multimodal": True,
                    "multimodal_pattern": "8modes",
                    "p_jump": 0.5,
                    "init_method": "fisher",
                    "skip_fisher": False,
                    "n_temps": 10,
                    "output_raw":False # Avoid raw file save (740KB vs. 59 MB)
    }

    # Loop over variablez for each Json.
    JsonFiles = []
    for iiLoop in range(n):
        for iiCut in range(len(T_obs_end_to_mergers)):
            # Copy then update the source class kwargs
            Merger_kwargs = copy.deepcopy(Merger_base_kwargs)
            Merger_kwargs["DeltatL_cut"] = T_obs_end_to_mergers[iiCut]
            Merger_kwargs["q"] = 1./rand_massratios[iiLoop][0]
            Merger_kwargs["m2"] = Merger_kwargs["M"]/(1.+Merger_kwargs["q"])
            Merger_kwargs["m1"] = Merger_kwargs["M"]-Merger_kwargs["m2"]
            Merger_kwargs["chi1"] = rand_spins[iiLoop][0]
            Merger_kwargs["chi2"] = rand_spins[iiLoop][1]
            Merger_kwargs["inc"] = rand_angles[iiLoop][0]
            Merger_kwargs["phi"] = rand_angles[iiLoop][1]
            Merger_kwargs["lambda"] = rand_angles[iiLoop][2]
            Merger_kwargs["beta"] = rand_angles[iiLoop][3]
            Merger_kwargs["psi"] = rand_angles[iiLoop][4]

            # Set output file names
            OutFileName="Randomized_angles_spins_MRat_"
            OutFileName+=str(iiLoop+1)+"_"+T_obs_labels[iiCut]
            RunTimekwargs["out_file"]=OutFileName

            # Create source object
            source=SYSs.SMBH_Merger(**Merger_kwargs)

            # Write params to json file
            SYU.WriteParamsToJson(source,LISA,inference_params,
                                      is_master,**RunTimekwargs)
            JsonFiles+=[source.JsonFile]


# Return list of jsonfiles to shell script
if is_master:
    for JsonFile in JsonFiles: print(JsonFile)
if use_mpi: comm_global.Barrier()
if use_mpi: sys.exit(0)
