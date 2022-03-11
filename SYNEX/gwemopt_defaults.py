"""


GWEMOPT parameters dictionary definition -- all of their flags are here for convenience...

Definitions can be found here: https://iopscience.iop.org/article/10.3847/1538-4357/838/2/108/meta


"""

from SYNEX.SYNEX_Utils import SYNEX_PATH
import numpy as np

# Define flags that are either not useful to Athena and LISA, or we want to default to False.
go_params_default={
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
"doObservability":False, # Order tile priorities by airmass etc induced changes to observability for telescope... Not relevant for Athena
"doSkybrightness":False,
"doEfficiency":False,
"doTransients":False,
"doSingleExposure":False, # I think this means single exposure time - i.e. one latency time per tile. I False time times are weighted by their skymap probability (sum of map_struct pixel probs)
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
"doAvoidGalacticPlane":False,
"doTrueLocation":False # Not sure what this is used for...
}



# Define flags that are useful to default to True.
go_params_default["doSkymap"] = True # We always want to do this
go_params_default["doCoverage"] = True # Not sure what it is yet because it's linked to scheduler and I'm not that far yet.
go_params_default["doSchedule"] = True # We always want to do this
go_params_default["doTiles"] = True # We always want to calculate these unless we already have them, in which case we can set this to False.
go_params_default["doMinimalTiling"] = True # This is a flag to tile only 90% CL area of skymap instead of the whole sky - if False the runtime is huge and plotting tesselation shows the whole sky is tiled
go_params_default["doCalcTiles"] = True # We always want to calculate these unless we already have them, in which case we can set this to False.


# Define remaining fields that are dependent on out implementation (e.g. paths etc)
go_params_default["skymap"] = None
go_params_default["outputDir"] = SYNEX_PATH+"/gwemopt_output"
go_params_default["tilingDir"] = SYNEX_PATH+"/Tile_files"
go_params_default["catalogDir"] = SYNEX_PATH+"/gwemopt_catalogs"
go_params_default["event"] = "IdeaPaperSystem"
go_params_default["coverageFiles"] = SYNEX_PATH+"/gwemopt_cover_files/Athena_test.dat"
go_params_default["lightcurveFiles"] = SYNEX_PATH+"/gwemopt/lightcurves/Me2017_H4M050V20.dat" ### THIS NEEDS TO BE CHANGED LATER WHEN WE HAVE SOME LIGHTCURVES...
go_params_default["tilesType"] = "moc" #  Tiling options are moc/greedy/hierarchical/ranked/galaxy.
go_params_default["scheduleType"] = "greedy" # Scheduling options are greedy/sear/weighted/airmass_weighted, or with _slew.
go_params_default["timeallocationType"] = "powerlaw" # "absmag" / "powerlaw" / "waw" / "manual" / "pem"
go_params_default["configDirectory"] = None # SYNEX_PATH+"/gwemopt_conf_files" # Is this needed? I don't think so...
go_params_default["gpstime"] = None # 1703721618.0 # 01/01/2034 00:00:00.000 UTC -- time of event -- set by source and if not by detectors
go_params_default["Ninj"] = 1000
go_params_default["Ndet"] = 1
go_params_default["Ntiles"] = 50 # Speific to tilesType=hierarchical and greedy
go_params_default["Ntiles_cr"] = 0.70 # Speific to tilesType=hierarchical and greedy (I think)
go_params_default["Dscale"] = 1.0
go_params_default["nside"] = 128
go_params_default["Tobs"] = None # Source param to be defined when data file is provided -- e.g. np.array([0.0,1.0]) = [Tstart,Tend], for times in DAYS
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
go_params_default["max_nb_tiles"] = np.array([-1,-1,-1])            ### What is this?
go_params_default["mindiff"] = 0.0 # minimum time in seconds between the start of exposure segment and the telesscope's tile segment
go_params_default["airmass"] = 2.5
go_params_default["iterativeOverlap"] = 0.0 # Only used when 'slicing' tiles (eg in gou.slice_map_tiles())
go_params_default["maximumOverlap"] = 1.0
go_params_default["catalog_n"] = 1.0
go_params_default["galaxy_grade"] = "S"
go_params_default["AGN_flag"] = True    ######### What does this change in GWEMOPT?
go_params_default["splitType"] = "regional"
go_params_default["Nregions"] = 768
go_params_default["Ncores"] = None # parallelizing moc creation where FoV to MoC function is used for whole telescope tesselation.
go_params_default["Nblocks"] = 4
go_params_default["unbalanced_tiles"] = None
go_params_default["treasuremap_token"] = ""
go_params_default["treasuremap_status"] = ["planned","completed"]
go_params_default["graceid"] = "S190426c"
go_params_default["raslice"] = [0.0,24.0] #### Can we set this to None? or empty list?
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


##### Parameters required in go_params but SYNEX handles differently. These are calculated in SYNEX_Utils functions
##### when arrays of sources/detectors are configured into gwemopt compatible objects...
# go_params_default["exposuretimes"] = np.array([30.0,30.0,30.0]) ### Time per tile? Latency time? -- Calculated when go_params for GWEMOPT is configured using each detector(s) config
# go_params_default["telescopes"] = np.array(["Athena_test_1","Athena_test_2","Athena_test_3"])
# go_params_default["do3D"]=False



# Create default config struct
# NB: ATHENA WILL FOLLOW L1 ORBIT (FOR THERMAL STABILITY AND VISIBILITY), WHICH COULD
# BE TURNED INTO A TIME VARYING LAT/LONG/ELEVATION IF THIS WORKS INSIDE GWEMOPT
# FUNCTIONS BUT NEED TO CHECK THIS...
config_struct_default = {
"telescope" : "Athena_test",
"filt" : "c",
"magnitude" : 18.7,
"exposuretime" : 10000.,    ### IN SECONDS, INCLUDES SLEW AND READONOUT TIME!
"min_observability_duration" : 0., ### IN HOURS -- minimum time per time if we set doSingleExposure=False and weight tile times by their probability
"latitude" : 0., # None, # 20.7204,       ### None if we want a telesscopic orbit?
"longitude" : 0., # None, # -156.1552,    ### None if we want a telesscopic orbit?
"elevation" : 0., # None, # 3055.0,       ### None if we want a telesscopic orbit? GWEMOPT uses these for airmass calcs... Ask to raise flag for this? It's always written to file but only used if you choose "airmass_weighted" scheduleType
"FOV_coverage" : 1., # In deg^2
"FOV" : 1., # In deg^2 ### Is it deg^2? line 26 of samplers.py suggests it's the radius/side length of the observation window...
"FOV_coverage_type" : "square",
"FOV_type" : "square",
"tesselationFile" : None,
"slew_rate" : 1., # in s/deg
"readout" : 6,                      ### In seconds?
"horizon" : None, # 30.,            ### None if we want a telesscopic orbit
"overhead_per_exposure" : 10., # Settle time after each slew/per tile? in seconds or what?
"filt_change_time" : 0.,
"sat_sun_restriction" : 45.,
"sat_earth_constraint" : 30.,
"sat_moon_constraint" : 20.0,
"orbitFile" : None, # Where to save orbit to
"frozenAthena" : False, # Consider stuck at L2 exactly so we only take Earth, Sun and Moon orbital motions into account
"inc" : 60., # deg incline of orbital plane normal to Sun-Earth axis.
"MeanRadius" : 750000000., # meters (from earth-orbit normal axis)
"semi_maj" : 750000000., # equivalent to MeanRadius axis ONLY IF we say we are really orbiting the centre and not the focal point
"L2_from_Earth" : 1500000000., # meters (L2 from Earth)
"eccentricity" : 0.8,
"ArgPeriapsis" : 0., # In Degrees, angle of point of closest approach to FOCAL POINT IN ORBIT PLANE
"AscendingNode" : 0., # In DEGREES
"phi_0" : 0., # In DEGREeS, initial phase of Athena when measurments start
"period" : 180., # In days, for one complete halo orbit about L2
"gps_science_start" : 1703721618.0, # 01/01/2034 00:00:00.000 UTC -- gps start time of science meaasurements
"mission_duration" : 3. # in YEARS- total time from *start of science* to end of Athena mission
}
# "dec_constraint" : "-90,90",





###########
#
# exposuretimes if doSingleExposure=True, then it takes them one by one and calculates things
# exposuretime if doSingleExposure=False, then it takes it as a base and reduces/enlarges based on tile prob
#
#


#######
# Segments is a list of ligo time segments (like sets) for when the telescope is available to make a measurement.
# exposurelist is a similar list of segments calculated from segments that slices them into times when there is an exposure to the source location and reduces times by overhead for each orientation, time per tile, etc.
# These are first calculated for the detector based on location and time, and then it is calculated for each tile depending on time allocation and position of moon/sun.
#
# I think we have to write our own versions of segments.py and adjust the usage of observer objects. the telescope part is easily done, even
# if you modify the returned set of segments and then recalculate the exposures based on that... But doing this for tiles is going to be a pain.





# Coverage struct has "data" as numpy array with N*9 dimensions for N*[ra,dec,mjd,mag,exposureTime,field,prob,airmass,program_id]
# If more than one detector, then numpy array of (n_det*N)*9 dimensions for (n_det*N)*[ra,dec,mjd,mag,exposureTime,field,prob,airmass,program_id]

# ha_constraint?
# moon_constraint?
# sat_sun_restriction? -- angular distance between satelite and sun in degrees. Need to modify lines 263-283 of 'coverage.py' to account for this in the right frame of reference.
# Line 455 of 'scheduler.py' ::: segments.py: tile_struct[key]["segmentlist"] is a list of segments when the tile is available for observation
# dec_constraint?

### Use lisatools to convert SSB frame to LISA, and then adjust angles to be at L1 instead of traiing Earth.
### Adjust telescope segments manually to be just one segment including all observation time
### Adjust tile segments manually according to moon and sun location relative to telescope using astropy functions as in 'coverage.py'
### Check these are not then manually done anywhere else in gwemopt...

# Bottom of page 46 in https://arxiv.org/pdf/1903.04083.pdf states AXIS (NOT Athena...) needs 45 deg at least from Sun for resolution and thermal control...
# Page iv in https://arxiv.org/pdf/1207.2745.pdf -- L2 Halo orbit with 750,000 km amplitude.
# Page 64 in previous reference says 'orbital timescales of 1 to 300ks', but not sure if this is Athena or source related... (period = 16 mins - 3.47 days, so probs source...)
# Second lagrange point is approx 1.5 million km from Earth.
# If no proposed orbit files are available, can we use JW Telescope orbit files? It too is in L2 Halo orbit... Could be a good first estimate.
