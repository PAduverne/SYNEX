#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 18:06:04 2023

@author: duverne
"""

import os, sys, json
import numpy as np
from astropy.table import Table


# Define SYNEX_PATH variable, add to path if not already there
try:
    from SYNEX.SYNEX_Utils import SYNEX_PATH
except:
    SYNEX_PATH=os.popen('pwd').read().split("SYNEX")[0]+"SYNEX"
    sys.path.insert(1, SYNEX_PATH)


# Reading the table containing the CBC population
parameters_cbc = Table.read(SYNEX_PATH + '/astrophysical_population/parameters_cbc.dat',
                            format='ascii.commented_header')

# Setting some parameters for the simulations.
# We use four values for the mass ratio, four for the total mass and
# four for the time before merger
# mass_ratio = [1., 3., 10., 100.]
mass_ratio = [1.1]
# log_M = [4, 5, 6, 7]
log_M = [4]
# index = 0
# DeltatL_cut = [30*24*60*60, # 1 month before merger
#                 7*24*60*60, # 1 week before merger
#                 3*24*60*60, # 3 day before merger
#                 24*60*60] # 1 day before merger
DeltatL_cut = [3*24*60*60]

# Number of simulations per {q, M_tot} couple
n_sim = 1

# Total number of simulations
N_tot_sim = len(mass_ratio) * len(log_M) * n_sim

for i, q in enumerate(mass_ratio):
    for j, M in enumerate(log_M):
        for k in range(n_sim):

            # Computing the index to chose the correct raw
            # in the injection table
            index = i * len(mass_ratio) * n_sim + j * n_sim + k
            print(index)

            # create the repertory and sub-repertories containing the files
            rep_name = 'astrophysical_population/q_{}_Mtot_{}/{}'.format(q, M, k)
            path = os.path.join(SYNEX_PATH, rep_name)
            if not os.path.exists(path):
                os.makedirs(path)

            # Get the parameter raw from table
            parameters = Table(parameters_cbc[index])

            # Replacing the mass with the right values
            # based on the mass ratio q and M_tot
            M_tot = float("1e{}".format(M))
            parameters["mass1"] = q * M_tot / (q + 1)
            parameters["mass2"] = M_tot / (q + 1)

            # Adding the chi_plus and chi_minus mass-weighted spin-parameters
            parameters["chi_p"] = (parameters["spin1z"] * parameters["mass1"] +
                                   parameters["spin2z"] * parameters["mass2"]) / M_tot
            parameters["chi_m"] = (parameters["spin1z"] * parameters["mass1"] -
                                   parameters["spin2z"] * parameters["mass2"]) / M_tot


            # Saving the parameters as an astropy Table for each injection
            parameters.write(path + "/parameters_SMBBH.dat",
                             format='ascii.commented_header',
                             overwrite=True)

            # Preparing the files for all the PE before the merger
            for cut in DeltatL_cut:

                # Name of the existential crisis file
                crisis_file = path + "/existential_crisis_file_{}.dat".format(cut)

                # Name of the files containing the parameters
                # and the inference results
                FileName = "merger_file_DeltatL_cut_{}".format(cut)

                # Parameters of the merger
                Merger_kwargs = {
                    # Total *redshifted* mass M=m1+m2, solar masses
                    "M": M_tot,
                    # Mass ratio q=m1/m2
                    "q": q,
                    # Dimensionless spin component 1 along orbital momentum
                    "chi1": parameters["spin1z"][0],
                    # Dimensionless spin component 2 along orbital momentum
                    "chi2": parameters["spin2z"][0],
                    # Dimensionless mass-weighted spin positive component
                    "chim": parameters["chi_p"][0],
                    # Dimensionless mass-weighted spin negative component
                    "chip": parameters["chi_m"][0],
                    # Time shift of coalescence, s -- coalescence is at t0*yr + Deltat*s,
                    # t0 in waveform_params
                    "Deltat": 0.0,
                    # t0 - coalescence time
                    "t0": parameters["t0"][0],
                    # Luminosity distance, Mpc
                    "dist": parameters["distance"][0],
                    # Inclination, observer's colatitude in source-frame
                    "inc": parameters["inclination"][0],
                    # Phase, observer's longitude in source-frame
                    "phi": parameters["coa_phase"][0],
                    # Longitude in the sky
                    "lambda": parameters["lambda"][0],
                    # Latitude in the sky
                    "beta": parameters["beta"][0],
                    # Polarization angle
                    "psi": parameters["polarization"][0],
                    # Flag indicating whether angles and Deltat pertain to the
                    # L-frame or SSB-frame
                    "Lframe": True,
                    # Stop time before merger event
                    "DeltatL_cut": - cut,
                    # Where to save the files
                    "rep_path": path,
                    "ExistentialFileName": crisis_file,
                    "JsonFile": FileName,
                    "H5File": FileName,
                    "MUTATED": False
                    }

                # Inference params definition with prior ranges and types
                inference_params = {
                    "list_params": ["Mchirp",
                                    "q",
                                    "Deltat",
                                    # "chi1", # dimensionless spin component 1
                                    # "chi2", # dimensionless spin component 2
                                    "chip", # dimensionless spin component 1
                                    "chim", # dimensionless spin component 2
                                    "dist", # Distance
                                    "inc", # inclination
                                    "phi", # Phase
                                    "lambda", # Longitude
                                    "beta", # Latitude
                                    "psi"], # Polarization

                    "infer_params": ["Mchirp",
                                     "q",
                                     "chim", # dimensionless spin component 1
                                     "chip", # dimensionless spin component 2
                                     "dist", # Distance
                                     "inc", # inclination
                                     "phi", # Phase
                                     "lambda", # Longitude
                                     "beta", # Latitude
                                     "psi"], # Polarization

                    "params_range": [[1e3, 8e7],
                                     [0., 15.],
                                     [-1., 1.],
                                     [-1., 1.],
                                     [1e2, 2e5],
                                     [0.,np.pi],
                                     [-np.pi, np.pi],
                                     [-np.pi, np.pi],
                                     [-np.pi/2.,np.pi/2.],
                                     [0.,np.pi]],

                    "prior_type": ["uniform",
                                   "uniform",
                                   "uniform",
                                   "uniform",
                                   "uniform",
                                   "sin",
                                   "uniform",
                                   "uniform",
                                   "cos",
                                   "uniform"],
                    "wrap_params": None
                    }

                # Parameters governing MCMC inference
                # NB: samples used are SMALL to reduce runtime.
                # Larger numbers commented out are recommended
                # for proper inference runs.
                # (n_iter, n_walkers, burn_in, n_temps) = (2000, 64, 500, 10)
                RunTimekwargs = {
                    "print_info": True,
                    "n_walkers": 64, # 20, # 96,
                    "n_iter": 60,# 30, # 8000, # 15000, #
                    "burn_in": 20, #20, # 4000, # 5000, #
                    "autocor_method": "autocor_new",
                    "likelihood_method": "fresnel", # "residuals",
                    "thin_samples": True, # False,
                    "multimodal": True,
                    "multimodal_pattern": "8modes",
                    "p_jump": 0.5,
                    "init_method": "fisher",
                    "skip_fisher": False,
                    "n_temps": 10,
                    "output_raw":False,
                    "out_rep": path
                    }

                # Create a json file with all the necessary parameters for SYNEX
                json_files_dict = {}

                # add the merger parameters to the dictionnary
                json_files_dict["Merger_kwargs"] = Merger_kwargs

                # add the inference parameters to the dictionnary
                json_files_dict["inference_params"] = inference_params

                # add the run time kwarg to the dictionnary
                json_files_dict["RunTimekwargs"] = RunTimekwargs

                # Prepare the repertory for saving the json file
                json_rep = os.path.join(path, "inference_param_files")
                if not os.path.exists(json_rep):
                    os.makedirs(json_rep)

                # setup the name of the json file
                name = "synex_file_cut_{}.json".format(cut)
                json_files_name = os.path.join(json_rep, name)

                # Removing previously existing files
                if os.path.exists(json_files_name):
                    os.remove(json_files_name)

                # save the dictionnary as a json file
                with open(json_files_name, 'w') as f:
                    json.dump(json_files_dict, f, indent=2)
                    f.close()
