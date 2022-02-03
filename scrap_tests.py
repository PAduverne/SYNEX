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
import json
import healpy as hp
import gwemopt
from gwemopt import utils as gou
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pylab as pylab
from SYNEX.SYNEX_Utils import pylab_params
pylab.rcParams.update(pylab_params)
mpl.use('MacOSX')


########################### Example - Plotting Athena orbit ###########################

import ligo.segments as segments
from astropy.time import Time
import SYNEX.segments_athena as segs_a

t0 = '2034-06-01T00:00:00.00' # YYYY-MM-DDTHH:mm:SS.MS
t = Time(t0, format='isot', scale='utc').gps

# Athena_kwargs={"FOV":1.,"exposuretime":60.,"slew_rate":1., "telescope":"Athena_1",
#             "ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_1_base.dat"}
Athena_kwargs={
                "ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_1_base.dat",
                "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_1_dev.dat",
                "frozenAthena" : False, # False,
                "exposuretime" : 12*60*60,
                "inc" : 80., # 60., # In DEGREES, incline of orbital plane normal to Sun-Earth axis.
                "MeanRadius" : 750000000*1.5, # 750000000., # meters (from earth-orbit normal axis)
                "semi_maj" : 750000000.*1.5, # 750000000., # equivalent to MeanRadius axis ONLY IF we say we are really orbiting the centre and not the focal point
                "eccentricity" : 0.8, # 0.8
                "ArgPeriapsis" : 0., # In DEGREES, angle of point of closest approach to FOCAL POINT IN ORBIT PLANE
                "AscendingNode" : 0., # In DEGREES
                "phi_0" : 0., # in DEGREES, initial phase of Athena when measurments start
                "period" : 180., # 180., # In days, for one complete halo orbit about L2
                "gps_science_start" : t, # 1703721618.0, # 01/01/2034 00:00:00.000 UTC -- gps start time of science meaasurements
                "mission_duration" : 2.
               }

Athena_1=SYDs.Athena(**Athena_kwargs)

config_struct = segs_a.calc_telescope_orbit(Athena_1.detector_config_struct,SAVETOFILE=False)

# "orbitFile" default - build on data at start...
strname="Athena_" + "".join(t0.split("T")[0].split("-")) + "_" + str(int((Athena_kwargs["mission_duration"]*364.25)//1)) + "d_inc"+str(int(Athena_kwargs["inc"]//1))+"_R"+str(int(Athena_kwargs["MeanRadius"]//1e6))+"Mkm_ecc"+str(int(Athena_kwargs["eccentricity"]//0.1))
strname+="_ArgPeri"+str(int(Athena_kwargs["ArgPeriapsis"]//1))+"_AscNode"+str(int(Athena_kwargs["AscendingNode"]//1))+"_phi0"+str(int(Athena_kwargs["ArgPeriapsis"]//1))
strname+="_P"+str(int(Athena_kwargs["period"]//1))+"_frozen"+str(Athena_kwargs["frozenAthena"])+".png"
print(strname)

SYU.PlotOrbit(config_struct)




# ########################### Example - Get source from file and plot skymap ###########################
#
# # Source file we want to resurract - note we only specify everything afer .../inference_data or .../inference_param_files and don't need to give suffix either
# # FileName = "IdeaPaperSystem_9d"
#
# # Merger args
# Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_base.dat",
#                  "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_dev.dat"}
# # Merger_kwargs = {"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Source_Dicts/IdeaPaperSystem_9d_base.dat"}
#
# # Resurrect - either from lisabeta data or saved source file
# # Merger = SYU.GetSourceFromLisabetaData(FileName,**Merger_kwargs)
# Merger = SYSs.SMBH_Merger(**Merger_kwargs)
#
# # Plot the CL bounds
# # SYU.PlotInferenceLambdaBeta(Merger.H5File, bins=50, SkyProjection=True, SaveFig=False, return_data=False)
#
# # Plot
# # SYU.PlotSkyMapData(Merger,SaveFig=False,plotName=None)
#
# Make some test telescopes
# Athena_kwargs={"FOV":1.,"exposuretime":60.,"slew_rate":1., "telescope":"Athena_1",
#             "ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_1_base.dat"}
# Athena_kwargs={"ExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_1_base.dat",
#                "exposuretime" : 60.,
#                "elevation" : 1500000000., # In meters
#                "NewExistentialFileName":"/Users/baird/Documents/LabEx_PostDoc/SYNEX/Saved_Telescope_Dicts/Athena_1_dev.dat"}
# Athena_1=SYDs.Athena(**Athena_kwargs)
#
# # Test preparation of gwemopt objects
# # go_params,map_struct = SYU.PrepareGwemoptDicts(Merger,detectors=Athena_1)
#
# # Test tiling directly with gwemopt
# tiling_kwargs={}
# TileDict = SYU.TileWithGwemopt(Merger, detectors=Athena_1, **tiling_kwargs)






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
