#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 23:44:56 2023

@author: duverne
"""

import os, sys, argparse, json, glob, ptemcee, time
import numpy as np
from astropy.table import Table

# Define SYNEX_PATH variable, add to path if not already there
try:
    from SYNEX.SYNEX_Utils import SYNEX_PATH
except:
    SYNEX_PATH=os.popen('pwd').read().split("SYNEX")[0]+"SYNEX"
    sys.path.insert(1, SYNEX_PATH)

# synex
import SYNEX.SYNEX_Utils as SYU
import SYNEX.SYNEX_Sources as SYSs
import SYNEX.SYNEX_Detectors as SYDs
import SYNEX.SYNEX_Telescopes as SYTs
import SYNEX.segments_athena as segs_a

from astropy.time import Time




# Using MPI or not
try:
    from mpi4py import MPI
except ModuleNotFoundError:
    print('MPI not found')
    MPI = None

use_mpi=False
is_master = True
mapper = map
if MPI is not None:
    MPI_size = MPI.COMM_WORLD.Get_size()
    MPI_rank = MPI.COMM_WORLD.Get_rank()
    comm_global = MPI.COMM_WORLD
    use_mpi = (MPI_size > 1)
    if use_mpi:
        # print("MPI rank/size: %d / %d" % (MPI_rank, MPI_size), flush=True)
        pool = ptemcee.mpi_pool.MPIPool(debug=False)
        is_master = pool.is_master()
        mapper = pool.map
    else:
        # print("No MPI", flush=True)
        is_master = True
        mapper = map


# Instanciate the LISA objects for the inference
# LISA object params
LISA_base_kwargs = {"TDI":'TDIAET', "verbose":False}

# Create LISA
LISA = SYDs.LISA(**LISA_base_kwargs)


def main():

    """Run a lisabeta inference on a given injection."""

    parser = argparse.ArgumentParser(description="Run lisabeta for a given simulated SMBBH.")

    parser.add_argument("--path",
                        dest="path",
                        required=True,
                        type=str,
                        help="Path where the SYNEX file is stored.")

    parser.add_argument("--cut",
                        dest="cut",
                        default="259200",
                        type=str,
                        help="Time shift before merger in seconds.\n"
                             "Default: 259200s, ie: 3j")

    # parser.add_argument("--skip-file-generation",
    #                     required=False,
    #                     action="store_true",
    #                     help="Skip parameter file generation if it already \n"
    #                          "exists.")

    parser.add_argument("--skip-inference",
                        dest="skip",
                        required=False,
                        action="store_true",
                        help="Skip the inference if the result file already\n"
                             "exists")

    args = parser.parse_args()

    path = args.path

    # get all the synex files for the simulation
    synex_files = glob.glob(path + "/synex_file_cut_{}.json".format(args.cut))

    content_directory = os.path.dirname(path)
    parameters = Table.read(content_directory + "/parameters_SMBBH.dat",
                            format='ascii.commented_header')

    for synex_file in synex_files:
        with open(synex_file, 'r') as file:

            print("Reading the SYNEX parameter file.")
            synex_param = json.load(file)

            # get the merger parameters
            Merger_kwargs = synex_param["Merger_kwargs"]

            # get the inference parameters
            inference_params = synex_param["inference_params"]

            # get the run time kwarg paralmeters
            RunTimekwargs = synex_param["RunTimekwargs"]

            # Create source parameter files
            Merger = SYSs.SMBH_Merger(**Merger_kwargs)

            # Write params to json file for lisabeta
            print("Writing the parameter in Lisabeta files")
            SYU.WriteParamsToJson(Merger, LISA, inference_params,
                                  is_master, **RunTimekwargs)

            print("Estimating the optimal SNR.")
            SNR = SYU.ComputeSNR(Merger, LISA, freqs=None,
                                 Lframe=False, ReturnAllVariable=False)
            parameters["SNR"] = SNR
            print("OPTIMAL SNR = {}".format(SNR))

            # Begin inference if the test_file doesn't already exist (this takes time...)
            script_path = os.path.join(SYNEX_PATH,
                                       'lisabeta/lisabeta/inference/ptemcee_smbh.py')

            # print("python3 " + script_path + " " + Merger.JsonFile)
            if os.path.isfile(Merger.H5File) and not args.skip:
                print("Starting the inference with Lisabeta.")
                t0 = time.time()
                os.system(("python3 " + # Call the right python version
                           script_path + " " + # Path for the inference script
                           Merger.JsonFile)) # Inference parameters
                t1 = time.time()
                print("Time to complete the inference: ", (t1-t0)/60, "min")

            # See inference results
            print("Creating the corner plot.")
            SYU.PlotInferenceData(Merger.H5File, extension='pdf')

            # See just lambda and beta (equivalent of RA and Dec)
            # Note: red is true location, blue is best fit posterior.
            print("Plotting the skymap.")
            SYU.PlotInferenceLambdaBeta(Merger.H5File, bins=50,
                                        SkyProjection=False,
                                        SaveFig=True,
                                        return_data=False)

            #setup the skymap name
            print("Creating the fits file for the skymap.")
            Merger.sky_map = (Merger.H5File.split("inference_data")[0] +
                              'Skymap_files' +
                              Merger.H5File.split("inference_data")[-1].strip(".h5"))

            # Create the FITS file for the skymap
            Merger.CreateSkyMapStruct(SkyMapFileName=None)

            # See the skymap
            print("Plot the skymap as a Mollweid projection.")
            save_dir = os.path.join(content_directory, 'Plots')
            cut = os.path.basename(synex_file).split('.')[0].split('_')[-1]
            SYU.PlotSkyMapData(Merger, save_path=save_dir,
                               plot_name="skymap_proba_density_{}.pdf".format(cut))
            print('Whouhouuuu')

    return 0

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover