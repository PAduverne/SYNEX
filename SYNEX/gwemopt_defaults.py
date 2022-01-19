"""


GWEMOPT parameters dictionary definition -- all of their flags are here for convenience...


"""

from SYNEX.SYNEX_Utils import SYNEX_PATH

# Define flags that are either not useful to Athena and LISA, or we want to default to False.
go_params_default={
"do3D":False,
"doEvent":False,
"doSuperSched":False,
"doMovie_supersched":False,
"doSamples":False,
"doPlots":False,
"doDatabase":False,
"doMovie":False,
"doIterativeTiling":False,
"doOverlappingScheduling":False,
"doPerturbativeTiling":False, # Not sure what this is yet...
"doOrderByObservability":False,
"doCatalog":False,
"doUseCatalog":False,
"doCatalogDatabase":False,
"doObservability":False,
"doSkybrightness":False,
"doEfficiency":False,
"doTransients":False,
"doSingleExposure":False, # I think this means single exposure time - i.e. one latency time per tile.
"doAlternatingFilters":False,
"doMaxTiles":False,
"doReferences":False,
"doChipGaps":False, # This is for a specific telescope that we don't care about
"doUsePrimary":False,
"doUseSecondary":False,
"doSplit":False,
"doParallel":False,
"writeCatalog":False,
"doFootprint":False,
"doBalanceExposure":False,
"doBlocks":False,
"doUpdateScheduler":False, # I guess if we have an improving skymap.... We can write functions to do this later.
"doTreasureMap":False,
"doRASlice":False,
"doRASlices":False,
"doRotate":False,
"doMindifFilt":False,
"doTrueLocation":False,
"doAvoidGalacticPlane":False
}



# Define flags that are useful to default to True.
go_params_default["doSkymap"] = True # We always want to do this
go_params_default["doCoverage"] = True # Not sure what it is yet because it's linked to scheduler and I'm not that far yet.
go_params_default["doSchedule"] = True # We always want to do this
go_params_default["doTiles"] = True # We always want to calculate these unless we already have them, in which case we can set this to False.
go_params_default["doMinimalTiling"] = True # I think this means minimize overlap of tiles when calculating tesselation...
go_params_default["doCalcTiles"] = True # We always want to calculate these unless we already have them, in which case we can set this to False.


# Define remaining fields that are dependent on out implementation (e.g. paths etc)
go_params_default["skymap"] = None
go_params_default["outputDir"] = SYNEX_PATH+"/gwemopt_output"
go_params_default["tilingDir"] = SYNEX_PATH+"/Tile_files"
go_params_default["catalogDir"] = SYNEX_PATH+"/gwemopt_catalogs"
go_params_default["event"] = "IdeaPaperSystem"
go_params_default["coverageFiles"] = SYNEX_PATH+"/gwemopt_cover_files/Athena_test.dat"
go_params_default["lightcurveFiles"] = SYNEX_PATH+"/gwemopt/lightcurves/Me2017_H4M050V20.dat" ### THIS NEEDS TO BE CHANGED LATER WHEN WE HAVE SOME LIGHTCURVES...
go_params_default["tilesType"] = "moc"
go_params_default["scheduleType"] = "greedy"
go_params_default["timeallocationType"] = "powerlaw"
go_params_default["configDirectory"] = SYNEX_PATH+"/gwemopt_conf_files"
go_params_default["gpstime"] = 1703721618.0 # 01/01/2034 00:00:00.000 UTC -- Athena launch
go_params_default["Ninj"] = 1000
go_params_default["Ndet"] = 1
go_params_default["Ntiles"] = 10
go_params_default["Ntiles_cr"] = 0.70
go_params_default["Dscale"] = 1.0
go_params_default["nside"] = 256
go_params_default["Tobs"] = np.array([0.0,1.0])
go_params_default["powerlaw_cl"] = 0.9
go_params_default["powerlaw_n"] = 1.0
go_params_default["powerlaw_dist_exp"] = 0
go_params_default["galaxies_FoV_sep"] = 1.0
go_params_default["footprint_ra"] = 30.0            ### What is this?
go_params_default["footprint_dec"] = 60.0           ### What is this?
go_params_default["footprint_radius"] = 10.0        ### What is this?
go_params_default["transientsFile"] = SYNEX_PATH+"/gwemopt/data/GW190425/transients.dat" # Although we probably wont use this...
go_params_default["transients_to_catalog"] = 0.8
go_params_default["dt"] = 14.0
go_params_default["galaxy_catalog"] = "GLADE"
go_params_default["filters"] = ['r','g','r']
go_params_default["exposuretimes"] = np.array([30.0,30.0,30.0])     ### Time per tile? Latency time?
go_params_default["max_nb_tiles"] = np.array([-1,-1,-1])            ### What is this?
go_params_default["mindiff"] = 0.0
go_params_default["airmass"] = 2.5
go_params_default["iterativeOverlap"] = 0.0
go_params_default["maximumOverlap"] = 1.0
go_params_default["catalog_n"] = 1.0
go_params_default["galaxy_grade"] = "S"
go_params_default["AGN_flag"] = False
go_params_default["splitType"] = "regional"
go_params_default["Nregions"] = 768
go_params_default["Ncores"] = 4
go_params_default["Nblocks"] = 4
go_params_default["unbalanced_tiles"] = None
go_params_default["treasuremap_token"] = ""
go_params_default["treasuremap_status"] = ["planned","completed"]
go_params_default["graceid"] = "S190426c"
go_params_default["raslice"] = [0.0,24.0]
go_params_default["nside_down"] = 2
go_params_default["max_filter_sets"] = 4
go_params_default["absmag"] = -15
go_params_default["phi"] = 0
go_params_default["theta"] = 0
go_params_default["program_id"] = -1
go_params_default["galactic_limit"] = 15.0
go_params_default["true_ra"] = None             # Source param to be defined when data file is provided
go_params_default["true_dec"] = None            # Source param to be defined when data file is provided
go_params_default["true_distance"] = None       # Source param to be defined when data file is provided
go_params_default["telescopes"] = "Athena_test"




# Create default config struct
# NB: ATHENA WILL FOLLOW L1 ORBIT (FOR THERMAL STABILITY AND VISIBILITY), WHICH COULD
# BE TURNED INTO A TIME VARYING LAT/LONG/ELEVATION IF THIS WORKS INSIDE GWEMOPT
# FUNCTIONS BUT NEED TO CHECK THIS...
config_struct_default = {
"telescope" : "Athena_test",
"filt" : "c",
"magnitude" : 18.7,
"exposuretime" : np.array([10000.]), # 10^4 s... Does gwemopt require this in seconds or hours or what?
"latitude" : 20.7204,       ### this could be a problem... Need to understand how gwemopt uses telesope location...
"longitude" : -156.1552,    ### this could be a problem... Need to understand how gwemopt uses telesope location...
"elevation" : 3055.0,       ### this could be a problem... Need to understand how gwemopt uses telesope location...
"FOV_coverage" : 1., # In deg^2
"FOV" : 1., # In deg^2
"FOV_coverage_type" : "square",
"FOV_type" : "square",
"tesselationFile" : SYNEX_PATH+"/gwemopt_tess_files/Athena_test.tess",
"slew_rate" : 1., # in s/deg
"readout" : 6,
"horizon" : 30,       ### this could be a problem... Need to understand how gwemopt uses horizon...
"overhead_per_exposure" : 10., # Settle time after each slew/per tile? in seconds or what?
"filt_change_time" : 60
}
