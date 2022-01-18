import astropy.units as u
import astropy.stats as astat
import numpy as np
import SYNEX.SYNEX_Utils as SYU
import gwemopt
import json
import pickle
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 8, # 'x-large',
         'axes.labelsize': 8, # 'x-large',
         'xtick.labelsize': 4, # 'x-large',
         'ytick.labelsize': 4, # 'x-large'}
         'lines.markersize': 2}
pylab.rcParams.update(params)
import scipy

class LISA:
    """
    Class to make a Space Based Instrument using the Interferometer base class

    Parameters
    ----------



    """

    def __init__(self, **kwargs):
        for key,value in kwargs.items():
            if key == 'tmin':
                self.tmin = value
            elif key == 'tmax':
                self.tmax = value
            elif key == 'TDI':
                self.TDI = value
            elif key == 'order_fresnel_stencil':
                self.order_fresnel_stencil = value
            elif key == 'LISAconst':
                self.LISAconst = value
            elif key == 'responseapprox':
                self.responseapprox = value
            elif key == 'frozenLISA':
                self.frozenLISA = value
            elif key == 'TDIrescaled':
                self.TDIrescaled = value
            elif key == "LISAnoise":
                self.LISAnoise = value
            elif key == "InstrumentalNoise" and hasattr(self,"LISAnoise"):
                self.LISAnoise["InstrumentalNoise"] = value
            elif key == "InstrumentalNoise":
                self.LISAnoise = {}
                self.LISAnoise["InstrumentalNoise"] = value
            elif key == "WDbackground" and hasattr(self,"LISAnoise"):
                self.LISAnoise["WDbackground"] = value
            elif key == "WDbackground":
                self.LISAnoise = {}
                self.LISAnoise["WDbackground"] = value
            elif key == "WDduration" and hasattr(self,"LISAnoise"):
                self.LISAnoise["WDduration"] = value
            elif key == "WDduration":
                self.LISAnoise = {}
                self.LISAnoise["WDduration"] = value
        if not hasattr(self,"tmin"):
            self.tmin = None     # Included
        if not hasattr(self,"tmax"):
            self.tmax = None     # Included
        if not hasattr(self,"TDI"):
            self.TDI = "TDIAET"     # Included
        if not hasattr(self,"order_fresnel_stencil"):
            self.order_fresnel_stencil = 0     # Included
        if not hasattr(self,"LISAconst"):
            self.LISAconst = "Proposal"     # Included
        if not hasattr(self,"responseapprox"):
            self.responseapprox = "full"     # Included
        if not hasattr(self,"frozenLISA"):
            self.frozenLISA = False     # Included
        if not hasattr(self,"TDIrescaled"):
            self.TDIrescaled = True     # Included
        if not hasattr(self,"LISAnoise"):
                self.LISAnoise = {
                "InstrumentalNoise": "SciRDv1",
                "WDbackground": True,
                "WDduration" : 3.0,
                'lowf_add_pm_noise_f0': 0.0,
                'lowf_add_pm_noise_alpha': 2.0}     # Included
        if not hasattr(self.LISAnoise,"InstrumentalNoise"):
                self.LISAnoise["InstrumentalNoise"] = "SciRDv1"
        if not hasattr(self.LISAnoise,"WDbackground"):
                self.LISAnoise["WDbackground"] = True
        if not hasattr(self.LISAnoise,"WDduration"):
                self.LISAnoise["WDduration"] = 3.0








class Athena:
    """
    Class to make a Space Based Instrument using the Interferometer base class

    Parameters
    ----------
    FoView : Float
        Field of view of a hypothetical Athena

    T_lat : Float
        Latency time of the tile for each move

    T_repos : Float
        Time to re-orient between each latency time

    T_init : Float
        Time to initialize Athena. E.g. if prior to the GW alert it was pointing
        elsewhere in the sky or in the middle of commissioning/calibration runs
        then this time gives the delay before tiling can begin.
    """

    def __init__(self, **kwargs):
        # Import default gwemopt dict
        from SYNEX.gwemopt_defaults import go_params_default, config_struct_default as go_params, config_struct

        # Change all things in go_params that are specified in kwargs,
        # and warn the user if any of the things they set are not used...
        print_reminder = False
        for key,value in kwargs.items():
            if key in go_params:
                go_params[key]=value
            else:
                print("'",key,"' not contained in gwemopt 'params' dict...")
                print_reminder = True
        # Now initiate config struct for this telescope
        if "ConfigFileName" in kwargs.keys():
            config_struct = gwemopt.utils.readParamsFromFile(kwargs["ConfigFileName"])
            config_struct["telescope"] = go_params["telescopes"]
        else:
            for key,value in kwargs.items():
                if key in config_struct:
                    config_struct[key]=value
                else:
                    print("'",key,"' not contained in gwemopt 'config_struct' dict...")
                    print_reminder = True
        # go_params["config"][go_params["telescopes"]]=config_struct

        # Check that file names are all coherent with telescope name
        go_params, config_struct = SYU.GWEMOPTPathChecks(go_params,config_struct)

        # Now calculate tesselation for telescope
        if go_params["doSingleExposure"]:
            exposuretime = np.array(go_params["exposuretimes"].split(","),dtype=np.float)[0]

            nmag = -2.5*np.log10(np.sqrt(config_struct["exposuretime"]/exposuretime))
            config_struct["magnitude"] = config_struct["magnitude"] + nmag
            config_struct["exposuretime"] = exposuretime
        if "tesselationFile" in config_struct:
            if not os.path.isfile(config_struct["tesselationFile"]):
                if config_struct["FOV_type"] == "circle":
                    gwemopt.tiles.tesselation_spiral(config_struct)
                elif config_struct["FOV_type"] == "square":
                    gwemopt.tiles.tesselation_packing(config_struct)
            if go_params["tilesType"] == "galaxy":
                config_struct["tesselation"] = np.empty((3,))
            else:
                config_struct["tesselation"] = np.loadtxt(config_struct["tesselationFile"],usecols=(0,1,2),comments='%')

        if "referenceFile" in config_struct:
            from astropy import table
            refs = table.unique(table.Table.read(
                config_struct["referenceFile"],
                format='ascii', data_start=2, data_end=-1)['field', 'fid'])
            reference_images =\
                {group[0]['field']: group['fid'].astype(int).tolist()
                for group in refs.group_by('field').groups}
            reference_images_map = {1: 'g', 2: 'r', 3: 'i'}
            for key in reference_images:
                reference_images[key] = [reference_images_map.get(n, n)
                                         for n in reference_images[key]]
            config_struct["reference_images"] = reference_images

        location = astropy.coordinates.EarthLocation(config_struct["longitude"],config_struct["latitude"],config_struct["elevation"])
        observer = astroplan.Observer(location=location)
        config_struct["observer"] = observer

        # Issue reminder of where to find list of gwemopt variables and flags
        if print_reminder:
            print("Some keys given at initiatiation of Athena class are not contained in gwemop params - see 'SYNEX/gwemopt_defaults.py' for full list of possible field names.")

        # Now set class attributes as dicts - then we can integrate these with many classes in utils when running gwemopt
        self.go_params = go_params
        self.config_struct = go_params

    def TileSkyArea(self,LISAPosteriorDataFile,TilePickleName=None,TileStrat="MaxProb", overlap=0.5, SAVE_SOURCE_EM_PROPERTIES_IN_TILE_JSON=True, go_params=None):
        """
        Function to load and tile a posterior sky area, saving a file of sky positions
        in order of a specification for a given tiling strategy.

        TileStrat : STRING
            - "MaxProb"
                Most basic tiling algorithm that selects locations in descending posterior probability.
                No asumptions (for the time being ) are made about slew rates, dead time, reboot times, etc.

            - "greedy"
                GWEMOPT tiling strategy named 'greedy'.
        """

        # Security check on filename - in gwemopt tiling methods this will be checked again... Need to optimize layout of checks.
        JsonFileLocAndName,H5FileLocAndName=SYU.CompleteLisabetaDataAndJsonFileNames(LISAPosteriorDataFile)

        # Choose strategy for tiling
        if TileStrat=="MaxProbOld":
            """
            Most basic way to tile - find the max prob sky location, and start there.
            """
            # Params for run time
            extent = np.sqrt(self.FoView) # in radians to match lisabeta units - side length of FoV
            BreakFlag = False
            TileNo = 0 # start at 0 to comply with gwemopt convention

            # Get data
            _, _, _, X, _, Y, _, Z = SYU.PlotInferenceLambdaBeta(H5FileLocAndName, bins=50, SkyProjection=False, SaveFig=False, return_data=True)
            betas = Y.flatten()

            # Repeat binning to ensure you can get the right points for overlap
            if betas[-1]-betas[0]>extent*overlap:
                bins = 2*int((betas[-1]-betas[0])/(extent*overlap))
                _, _, _, X, _, Y, _, Z = SYU.PlotInferenceLambdaBeta(H5FileLocAndName, bins=bins, SkyProjection=False, SaveFig=False, return_data=True)
            Post_probs = Z.flatten()
            lambdas = X.flatten()
            betas = Y.flatten()
            plt.close('all')

            # Output dictionary of tile properties
            TileDict = {"Tile Strat": TileStrat, "overlap": overlap, "tile_structs":{self.telesope:{}}, "LISA Data File":H5FileLocAndName}

            # Run through the space of confidence area, extracting tile properties and reducing the list of total tiles
            while BreakFlag == False:
                tile_prob = np.max(Post_probs)
                arg_max = np.argmax(Post_probs)
                if not isinstance(arg_max,np.int64):
                    arg_max = arg_max[0]
                tile_beta = betas[arg_max]
                tile_lambda = lambdas[arg_max]

                # Put in the tile to the tile dictionary
                lambda_range=[tile_lambda-extent/2.,tile_lambda+extent/2.]
                beta_range=[tile_beta-extent/2.,tile_beta+extent/2.]
                ra_range=[np.rad2deg(lambda_range[0]+np.pi),np.rad2deg(lambda_range[1]+np.pi)]
                dec_range=[np.rad2deg(beta_range[0]),np.rad2deg(beta_range[1])]
                corners=np.array([[ra_range[0],dec_range[0]], ## Following convention of gwemopt ordering
                         [ra_range[1],dec_range[0]],
                         [ra_range[0],dec_range[1]],
                         [ra_range[1],dec_range[1]]])
                Tile = {"ra":np.rad2deg(tile_lambda+np.pi), ## gwemopt dict fields -- 'ra', 'dec', 'ipix', 'corners', 'patch', 'area', 'segmentlist'
                        "dec":np.rad2deg(tile_beta),
                        "ipix":[], ## this is a number but I'm not sure what it is yet...
                        "corners":corners,
                        "patch":[], ## this is a matplotlib 'path' thing but I'm not sure what it is yet...
                        "area":self.FoView*(180./np.pi)**2,
                        "segmentlist":[],
                        "lambda":tile_lambda, ## lisabeta dict fields
                        "beta":tile_beta,
                        "lambda_range":lambda_range,
                        "beta_range":beta_range,
                        "Center post prob": tile_prob}
                TileDict["tile_structs"][self.telescope][TileNo] = Tile

                # Delete coordinates of mode that's done minus some overlap
                mask = [np.all([betas[ii]<tile_beta+overlap*extent, betas[ii]>tile_beta-overlap*extent, lambdas[ii]<tile_lambda+overlap*extent, lambdas[ii]>tile_lambda-overlap*extent]) for ii in range(len(lambdas))]
                betas = np.delete(betas, mask, axis=0)
                lambdas = np.delete(lambdas, mask, axis=0)
                Post_probs = np.delete(Post_probs, mask, axis=0)

                # exit loop if the remaining data is
                if len(Post_probs) <= 100 or TileNo>100:
                    BreakFlag = True
                TileNo+=1

        elif any([TileStrat=="MaxProb", TileStrat=="moc" , TileStrat=="greedy" , TileStrat=="hierarchical" , TileStrat=="ranked" , TileStrat=="galaxy"]):
            """
            "MaxProb" is the same as "MaxProbOld", but using gwemopt features...
            """
            import gwemopt.tiles as gots
            from gwemopt import utils as gou

            # check if go_params given and create if not -- this will later be integrated into the athena init function...
            if go_params==None:
                go_params={
                            "nside":16,
                            "TileStrat":TileStrat, # Since we want to keep the SYNEX call and gwemopt calls different for added strats in SYNEX
                            "tilesType":TileStrat,
                            "exposuretimes":np.array([self.T_lat]), # This has the same length as number of filters.... Need to look at this for Athena.
                            "Ntiles":10, # set 200, 10 for debugging, total tiles for now and later we cut ones we dont need. this helps if we only change things like latency times later...
                            "iterativeOverlap":overlap,   #### WHAT IS THIS USED FOR?
                            "maximumOverlap":1.0          #### WHAT IS THIS USED FOR?
                            }
            if len(go_params["exposuretimes"])==1:
                go_params["doSingleExposure"]=True
            else:
                go_params["doSingleExposure"]=False
            go_params["iterativeOverlap"]=overlap

            # Make tile dict, which contains gwemopt 'tile_structs' object, plus some stuff for SYNEX and lisabeta combatability
            if TileStrat=="MaxProb":
                # Adjust gwemopt params for requested tiling strat -- this will later be integrated into the athena init function...
                go_params["tilesType"]="moc"
                go_params["timeallocationType"] = "powerlaw"
                go_params["powerlaw_n"]=1
                if "powerlaw_cl" not in go_params: # i.e. if this is not already specified at input..
                    go_params["powerlaw_cl"]=0.9
                go_params["powerlaw_dist_exp"]=0
                go_params["doMinimalTiling"]=True ## Can we ignore this later or somehow 'sense' if it's True or False?

                # Run checker to complete gwemopt fields
                go_params = SYU.gou_params_checker(go_params,detector=self)
                # Create sky_map dictionary
                go_params,map_struct = SYU.CreateSkyMapStruct(go_params,LisabetaFileName=H5FileLocAndName)

                # Get tiling structs
                moc_structs = gwemopt.moc.create_moc(go_params, map_struct=map_struct)
                tile_structs = gots.moc(go_params,map_struct,moc_structs)
            elif TileStrat=="moc":
                # Adjust gwemopt params for requested tiling strat -- this will later be integrated into the athena init function...
                go_params["timeallocationType"] = "powerlaw"

                # Run checker to ensure dict is complete according to gwemopt definitions
                go_params = SYU.gou_params_checker(go_params,detector=self)
                # Create the necessary sky_map dictionary
                go_params,map_struct = SYU.CreateSkyMapStruct(go_params,LisabetaFileName=H5FileLocAndName)

                # Get tiling structs
                moc_structs = gwemopt.moc.create_moc(go_params, map_struct=map_struct)
                tile_structs = gots.moc(go_params,map_struct,moc_structs)
            elif TileStrat=="greedy":
                # Adjust gwemopt params for requested tiling strat -- this will later be integrated into the athena init function...
                go_params["timeallocationType"] = "powerlaw"

                # Run checker to ensure dict is complete according to gwemopt definitions
                go_params = SYU.gou_params_checker(go_params,detector=self)
                # Create the necessary sky_map dictionary
                go_params,map_struct = SYU.CreateSkyMapStruct(go_params,LisabetaFileName=H5FileLocAndName)

                # Get tiling structs
                tile_structs = gots.greedy(go_params,map_struct)
                for telescope in go_params["telescopes"]:
                    go_params["config"][telescope]["tesselation"] = np.empty((0,3))
                    tiles_struct = tile_structs[telescope]
                    for index in tiles_struct.keys():
                        ra, dec = tiles_struct[index]["ra"], tiles_struct[index]["dec"]
                        go_params["config"][telescope]["tesselation"] = np.append(go_params["config"][telescope]["tesselation"],[[index,ra,dec]],axis=0)
            elif TileStrat=="hierarchical":
                # Adjust gwemopt params for requested tiling strat -- this will later be integrated into the athena init function...
                go_params["timeallocationType"] = "powerlaw"

                # Run checker to ensure dict is complete according to gwemopt definitions
                go_params = SYU.gou_params_checker(go_params,detector=self)
                # Create the necessary sky_map dictionary
                go_params,map_struct = SYU.CreateSkyMapStruct(go_params,LisabetaFileName=H5FileLocAndName)

                # Get tiling structs
                tile_structs = gots.hierarchical(go_params,map_struct)
                go_params["Ntiles"] = []
                for telescope in go_params["telescopes"]:
                    go_params["config"][telescope]["tesselation"] = np.empty((0,3))
                    tiles_struct = tile_structs[telescope]
                    for index in tiles_struct.keys():
                        ra, dec = tiles_struct[index]["ra"], tiles_struct[index]["dec"]
                        go_params["config"][telescope]["tesselation"] = np.append(go_params["config"][telescope]["tesselation"],[[index,ra,dec]],axis=0)
                    go_params["Ntiles"].append(len(tiles_struct.keys()))
            elif TileStrat=="ranked":
                # Adjust gwemopt params for requested tiling strat -- this will later be integrated into the athena init function...
                go_params["timeallocationType"] = "powerlaw"

                # Run checker to ensure dict is complete according to gwemopt definitions
                go_params = SYU.gou_params_checker(go_params,detector=self)
                # Create the necessary sky_map dictionary
                go_params,map_struct = SYU.CreateSkyMapStruct(go_params,LisabetaFileName=H5FileLocAndName)

                # Get tiling structs
                go_params["doMinimalTiling"]=True ###### ADDED THIS IN CASE THS HELPS THE CAUSE -- NEED TO VERIFY THAT THIS SAVES THE ROUTINE FROM BREAKING
                moc_structs = gwemopt.rankedTilesGenerator.create_ranked(go_params,map_struct)
                tile_structs = gots.moc(go_params,map_struct,moc_structs)
            elif TileStrat=="galaxy":
                # Adjust gwemopt params for requested tiling strat -- this will later be integrated into the athena init function...
                go_params["timeallocationType"] = "powerlaw"

                # Run checker to ensure dict is complete according to gwemopt definitions
                go_params = SYU.gou_params_checker(go_params,detector=self)
                # Create the necessary sky_map dictionary
                go_params,map_struct = SYU.CreateSkyMapStruct(go_params,LisabetaFileName=H5FileLocAndName)

                # Get tiling structs
                map_struct, catalog_struct = gwemopt.catalog.get_catalog(go_params, map_struct)
                tile_structs = gots.galaxy(go_params,map_struct,catalog_struct)
                for telescope in go_params["telescopes"]:
                    go_params["config"][telescope]["tesselation"] = np.empty((0,3))
                    tiles_struct = tile_structs[telescope]
                    for index in tiles_struct.keys():
                        ra, dec = tiles_struct[index]["ra"], tiles_struct[index]["dec"]
                        go_params["config"][telescope]["tesselation"] = np.append(go_params["config"][telescope]["tesselation"],[[index,ra,dec]],axis=0)

            # # Plot tiles using gwemopt directly...
            # gwemopt.plotting.tiles(go_params, map_struct, tile_structs)

            # Keeping the sub-dict used by gwemopt just in case (for now), insert some value conversions for SYNEX compatibility -- might get rid of this later and stick to gwemopt conventions...
            tile_structs[self.telescope]["Tile Strat"] = TileStrat
            for tileID in range(len(tile_structs[self.telescope].keys())):
                corners=tile_structs[self.telescope][tileID]['corners']
                lambda_range=[np.deg2rad(corners[0][0]),np.deg2rad(corners[1][0])]
                beta_range=[np.deg2rad(corners[0][1]),np.deg2rad(corners[1][1])]
                tile_structs[self.telescope][tileID]["lambda"]=np.deg2rad(tile_structs[self.telescope][tileID]['ra'])-np.pi
                tile_structs[self.telescope][tileID]["beta"]=np.deg2rad(tile_structs[self.telescope][tileID]['dec'])
                tile_structs[self.telescope][tileID]["lambda_range"]=lambda_range
                tile_structs[self.telescope][tileID]["beta_range"]=beta_range
                tile_structs[self.telescope][tileID]["Center post prob"]=tile_structs[self.telescope][tileID]['prob']

            # Output dictionary of tile properties
            TileDict = {"LISA Data File": H5FileLocAndName,
                        "tile_structs": tile_structs,
                        "go_params": go_params,
                        "map_struct": map_struct}

        # Add source properties to dictionary for faster kuiper step if requested
        if SAVE_SOURCE_EM_PROPERTIES_IN_TILE_JSON:
            # grab source parameters from corresponding json file
            source = SYU.GetSourceFromInferenceData(JsonFileLocAndName)

            # Recreate the EM flux and CTR - fime name and gamma hard coded for now... Maybe shoul dbe given at entry line to this function for continuity?
            source.GenerateEMFlux(self)
            ARF_file_loc_name = 'XIFU_CC_BASELINECONF_2018_10_10.arf'
            source.GenerateCTR(ARF_file_loc_name=ARF_file_loc_name,gamma=1.7)

            # Add to tile dict
            TileDict["source_EM_properties"]={
            "xray_flux" : source.xray_flux,
            "xray_time" : source.xray_time,
            "phi" : source.GW_phi,
            "Omega" : source.GW_Omega,
            "r" : source.r,
            "xray_gamma" : source.xray_gamma,
            "xray_phi_0" : source.xray_phi_0,
            "CTR" : source.CTR,
            "GW_freqs" : source.GW_freqs
            }

        # Make output file name if not specified - this will put the file in folder where user is executing their top most script
        if TilePickleName==None:
            TilePickleName = H5FileLocAndName.split("inference_data")[0] + "Tile_files" + go_params["skymap"].split("Skymap_files")[-1]
            TilePickleName = TilePickleName[:-3] + "_" + TileStrat + ".dat"
            print(TilePickleName)

        # Write to file
        with open(TilePickleName, 'wb+') as f:
            pickle.dump(TileDict, f)

        # Write the name of the file to the detector object
        self.TilePickleName = TilePickleName ### Pickle or dat files? Not sure which is better here...

    def GetKuiper(self, TilePickleFile, source=None):
        # Check first that the necessary steps have been done for the source and detector
        # if not hasattr(source, "CTR"):
        #     raise ValueError("Source object does not have a CTR - please call source.GenerateEMFlux followed by source.GenerateCTR...")
        # if not hasattr(detector, "TileJsonName"):
        #     raise ValueError("Detector object does not have a TileJsonName - please call detector.TileSkyArea...")

        # Figure out maximum times we can do in the time available
        # Will include time slew as an extra exit flag in the loop over tiles,
        # cutting out the last m tiles (m<n) once the slew time consumes the equivalent time for m tiles.
        # Define the time to merger
        # Get tiles from json file
        with open(TilePickleFile, 'rb') as f:
            TileDict = pickle.load(f)

        H5FileLocAndName = TileDict["LISA Data File"]
        json_file,H5FileLocAndName=SYU.CompleteLisabetaDataAndJsonFileNames(LISAPosteriorDataFile)
        with open(json_file, 'r') as f:
            input_params = json.load(f)
        f.close()

        if source==None:
            # convert the json params to source class
            source = SYU.ParamsToClasses(input_params)

            # Check if the EM properties were saved in the tile json file, otherwise recalculate. NOTE: recalculation is time intensive...
            if "source_EM_properties" in TileDict:
                source.xray_flux=TileDict["source_EM_properties"]["xray_flux"]
                source.xray_time=TileDict["source_EM_properties"]["xray_time"]
                source.GW_phi=TileDict["source_EM_properties"]["phi"]
                source.GW_Omega=TileDict["source_EM_properties"]["Omega"]
                source.r=TileDict["source_EM_properties"]["r"]
                source.xray_gamma=TileDict["source_EM_properties"]["xray_gamma"]
                source.xray_phi_0=TileDict["source_EM_properties"]["xray_phi_0"]
                source.CTR=TileDict["source_EM_properties"]["CTR"]
                source.GW_freqs=TileDict["source_EM_properties"]["GW_freqs"]
            else:
                # Recreate the EM flux and CTR - fime name and gamma hard coded for now... Maybe shoul dbe given at entry line to this function for continuity?
                source.GenerateEMFlux(self)
                ARF_file_loc_name = 'XIFU_CC_BASELINECONF_2018_10_10.arf'
                source.GenerateCTR(ARF_file_loc_name=ARF_file_loc_name,gamma=1.7, LatTime=self.T_lat)

        T_s = source.xray_time[0]
        if T_s>input_params["waveform_params"]["DeltatL_cut"]:
            TimeToMerger = T_s
        else:
            TimeToMerger = input_params["waveform_params"]["DeltatL_cut"]
        n_times = int(-TimeToMerger//self.T_lat)
        print(n_times, "tiles in remaining time to merger...")

        # Sort the dictionary to contain just the tiles we want
        Extra_keys = ["Tile Strat", "overlap", "LISA Data File", "source_EM_properties", "tile_structs", "go_params", "map_struct"]
        if TileDict["Tile Strat"]=="MaxProb":
            TileDict_reduced = {k: v for k, v in TileDict.items() if k not in Extra_keys and int(k)<n_times}
        else:
            for telescope in TileDict["tile_structs"].keys(): ######### NEED TO ADJUST THIS SO WE CAN HANDLE SEVERAL TELESCOPE INSTANCES... MAYBE WE DON'T NEED TO AND WE CAN FOCUS ON SEVERAL TILING METHODS AT A TIME INSTEAD?
                TileDict_reduced = {k: v for k, v in TileDict["tile_structs"][telescope].items() if int(k)<n_times} ######## I THINK YOU CAN REDUCE COMLEXITY HERE IF WE CAN FEED N_TILES DIRECTLY INTO GWEMOPT- IT HAS A PARAM IN go_params THAT INDICATES HOW MANY TILES TO CALCULATE...
        max_key = [int(k) for k in TileDict_reduced.keys()]
        max_key = max(max_key)
        if max_key<n_times:
            n_times = max_key
        self.tileIDs = [int(k) for k in TileDict_reduced.keys()]
        self.tile_times = [(ID-0.5)*self.T_lat+TimeToMerger for ID in self.tileIDs] # + slew dead time! ### This is the central time

        # Get photon arrival times, cut the majority of the xray time and CTR since we will start tiling at time_to_merger
        from scipy.stats import uniform
        from astropy.stats import kuiper,kuiper_two
        from functools import partial
        import random
        xray_time_merger = [time for time in source.xray_time if time>=TimeToMerger]
        xray_flux_merger = [f for (time,f) in zip(source.xray_time,source.xray_flux) if time>=TimeToMerger]
        phi_merger = [phi/2. for (time,phi) in zip(source.xray_time,source.GW_phi) if time>=TimeToMerger] # Should be orbital phase NOT GW... They say in the paper that they just 'infer' the phi_is from the measured gravitational wave phase evolution
        # phi_merger = [random.random()*2.*np.pi for ii in range(len(phi_merger))]
        # Omega_merger = [Om for (time,Om) in zip(source.xray_time,source.GW_Omega) if time>=TimeToMerger]
        CTR_merger = [CTR for (time,CTR) in zip(source.xray_time,source.CTR) if time>=TimeToMerger]
        # CTR_merger = [np.sin(2.*np.pi*0.0001*time)+1.5 for time in xray_time_merger]
        CTR_sum = sum(CTR_merger)
        probs_merger = [CTR/CTR_sum for CTR in CTR_merger]
        n_photons = int(np.trapz(CTR_merger,xray_time_merger))
        ts = xray_time_merger[0]
        te = ts + self.T_lat

        t_is_merger = np.random.choice(xray_time_merger, n_photons, p=probs_merger).tolist()
        t_is_merger.sort()

        two_pi = 2.*np.pi
        f = scipy.interpolate.interp1d(xray_time_merger,phi_merger) # Omega_merger) #
        phi_is_merger = list(f(t_is_merger)) # [2.*np.pi*0.0001*time for time in t_is_merger]

        # Add backgrounds
        ADD_BACKGROUND = True
        if ADD_BACKGROUND:
            import random
            CTR_bg = [7.4e-5]*len(CTR_merger) # 7.4e-5
            n_bg_photons = int(np.trapz(CTR_bg,xray_time_merger))
            print(n_bg_photons, "background photons added to the whole x-ray timeseries...")
            max_time = -xray_time_merger[0]
            time_lim = [time-xray_time_merger[0] for (time,CTR) in zip(xray_time_merger,CTR_merger) if CTR>0]
            time_lim = time_lim[-1]
            t_is_bg = [random.random()*time_lim-max_time for ii in range(n_bg_photons)] # np.random.choice(xray_time_merger, n_bg_photons).tolist() # default is uniform
            t_is_bg.sort()
            t_is_merger+=t_is_bg
            rand_angles = [random.random()*two_pi for ii in range(n_bg_photons)] # default is uniform
            phi_is_merger+=rand_angles

        bin_pops,bins=np.histogram([p_i%two_pi for p_i in phi_is_merger],bins=100) # int(len(phi_is_merger)/10)
        bin_centres = [bins[ii]+0.5*(bins[1]-bins[0]) for ii in range(len(bins)-1)]
        bin_pops_normed = [bin_pops[ii]/sum(bin_pops) for ii in range(len(bin_pops))]
        bin_pops_normed_cumsum = [sum(bin_pops_normed[:ii]) for ii in range(len(bin_pops_normed))]
        bin_pops_normed_mean = np.mean(bin_pops_normed)
        probs_merger_mean = np.mean(probs_merger)

        ###########    THIS NEEDS CLEANING UP !!    ###########
        kuipers = [0.1]
        self.n_photons = []
        self.n_exposures = 0
        exposure_photons = 0
        exposure_t_is = []
        exposure_phi_is = []
        self.exposure_tiles = []
        self.exposure_xray_time = []
        self.exposure_CTR = []
        self.exposure_tile_probs = []
        self.tile_kuiper_p_val_trace = []
        self.exposure_kuiper_trace = []
        self.exposure_kuiper_p_val_trace = []
        self.detection_kuiper_p_val_trace = []
        self.exposure_tile_xray_time = []
        self.exposure_tile_CTR = []
        self.exposure_tile_probs = []
        self.accum_exposure_photons_trace = []
        self.exposure_photons_trace = []
        self.tile_detection_p_val_trace = []
        exposure_kuiper = 0.01
        exposure_p_val = 1.
        tile_p_val = 1.
        n_photons = 0
        import time
        from operator import add
        t0 = time.time()
        for tile in TileDict_reduced:
            # See if the tile includes the source - update the statistics for the source
            if source.beta<TileDict_reduced[tile]['beta_range'][1] and source.beta>TileDict_reduced[tile]['beta_range'][0] and source.lamda<TileDict_reduced[tile]['lambda_range'][1] and source.lamda>TileDict_reduced[tile]['lambda_range'][0]: # 1: #
                # Calculate tile statistical properties
                tile_CTR = [CTR for (CTR,time) in zip(CTR_merger,xray_time_merger) if time>ts and time<te]
                tile_xray_time = [time for time in xray_time_merger if time>=ts and time<=te]
                t_is = [t_i_merger for t_i_merger in t_is_merger if t_i_merger>ts and t_i_merger<te]
                phi_is = [phi_i%two_pi for (phi_i,time) in zip(phi_is_merger,t_is_merger) if time>ts and time<te]
                n_photons = len(t_is)
                self.n_photons.append(n_photons)

                tile_kuiper, tile_p_val = kuiper(phi_is, partial(uniform.cdf, loc=0.,scale=two_pi)) # loc=min(phi_is),scale=max(phi_is)-min(phi_is))) #
                kuipers.append(tile_kuiper)

                # Update exposure statistics
                self.n_exposures+=1
                self.exposure_tiles.append(int(tile))
                exposure_photons += n_photons
                exposure_t_is += t_is
                exposure_phi_is += phi_is
                bin_pops, bins = np.histogram(exposure_phi_is,bins=np.linspace(0.,two_pi,50))
                bin_centres = [bins[ii]+0.5*(bins[1]-bins[0]) for ii in range(len(bins)-1)]
                s_phi = [bin_pops[ii]/sum(bin_pops) for ii in range(len(bin_pops))]
                S_phi = [sum(s_phi[:ii]) for ii in range(len(s_phi))]
                U_phi = [ii/len(bin_centres) for ii in range(len(bin_centres))]
                exposure_kuiper, exposure_p_val = kuiper(exposure_phi_is, partial(uniform.cdf, loc=0.,scale=two_pi))

                if self.n_exposures>0:
                    detection_p_val = 1.-(1.-exposure_p_val)**(self.n_exposures*518400)
                    tile_detection_p_val = 1.-(1.-tile_p_val)**(self.n_exposures*518400)
                else:
                    detection_p_val = 1.-(1.-exposure_p_val)
                    tile_detection_p_val = 1.-(1.-tile_p_val)
                print("Tile:", tile, self.n_exposures, n_photons, len(exposure_phi_is), tile_kuiper, tile_p_val, exposure_kuiper, exposure_p_val, detection_p_val)

            # Update the traced properties
            self.tile_kuiper_p_val_trace.append(tile_p_val)
            self.tile_detection_p_val_trace.append(tile_detection_p_val)
            self.exposure_kuiper_trace.append(exposure_kuiper)
            self.exposure_kuiper_p_val_trace.append(exposure_p_val)
            self.detection_kuiper_p_val_trace.append(detection_p_val)
            self.accum_exposure_photons_trace.append(exposure_photons)
            self.exposure_photons_trace.append(n_photons)

            # Update the start and end times for the tile latency
            ts = te # + slew time
            te = ts + self.T_lat # times are in seconds to merger

        # exposures
        print(self.n_exposures, "exposures of source by tiles", self.exposure_tiles, "using", TileDict["Tile Strat"], "tiling strategy.")
        # print("Total photons:", sum(self.n_photons), "with", sum(self.n_photons[:-1]), "exposure photons and", self.n_photons[-1],"background photons.")
        print("Accumulated photons during on-source exposures (S+B):", sum(self.n_photons))

        # Plot exposure photons as function of tiling
        PLOT_PHOTON_TRACE = False
        if PLOT_PHOTON_TRACE:
            Times = [t/(24.*60.*60.) for t in self.tile_times]
            plt.plot(Times, self.exposure_photons_trace)
            plt.ylabel(r"Exposure Photons")
            plt.xlabel("Time [days]")
            plt.show()
            plt.grid()

        # Plot exposure photons as function of tiling
        PLOT_EXP_PHOTON_TRACE = False
        if PLOT_EXP_PHOTON_TRACE:
            Times = [t/(24.*60.*60.) for t in self.tile_times]
            plt.plot(Times, self.accum_exposure_photons_trace)
            ax = plt.gca()
            ax.set_yscale('log')
            plt.ylabel(r"Accumulated Exposure Photons")
            plt.xlabel("Time [days]")
            plt.show()
            plt.grid()

        # Plot Kuiper statistic for exposures as function of tiling
        PLOT_KUIPERS = False
        if PLOT_KUIPERS:
            Times = [t/(24.*60.*60.) for t in self.tile_times]
            plt.plot(Times, self.exposure_kuiper_trace)
            plt.ylabel(r"$\mathcal{K}_{exposure}$")
            plt.xlabel("Time [days]")
            plt.show()
            plt.grid()

        # Plot Kuiper p-value for tiles as function of tiling
        PLOT_TILE_KUIPER_PVAL_TRACE = False
        if PLOT_TILE_KUIPER_PVAL_TRACE:
            Times = [t/(24.*60.*60.) for t in self.tile_times]
            plt.plot(Times,self.tile_kuiper_p_val_trace)
            ax = plt.gca()
            ax.set_yscale('log')
            plt.xlabel("Time [days]")
            plt.ylabel("Tile Kuiper p-value")
            plt.show()
            plt.grid()

        # Plot Detection p-value for tiles as function of tiling
        PLOT_TILE_DETECTION_PVAL_TRACE = False
        if PLOT_TILE_DETECTION_PVAL_TRACE:
            Times = [t/(24.*60.*60.) for t in self.tile_times]
            plt.plot(Times,self.tile_detection_p_val_trace)
            ax = plt.gca()
            ax.set_yscale('log')
            plt.xlabel("Time [days]")
            plt.ylabel("Tile Detection Kuiper p-value")
            plt.show()
            plt.grid()

        # Plot Kuiper p-value for exposures as function of tiling
        PLOT_EXP_KUIPER_PVAL_TRACE = False
        if PLOT_EXP_KUIPER_PVAL_TRACE:
            Times = [t/(24.*60.*60.) for t in self.tile_times]
            plt.plot(Times,self.exposure_kuiper_p_val_trace)
            ax = plt.gca()
            ax.set_yscale('log')
            plt.xlabel("Time [days]")
            plt.ylabel("Exposure Kuiper p-value")
            plt.show()
            plt.grid()

        # Plot Kuiper detection p-value for exposures as function of tiling (n_pix and n_exposures)
        PLOT_DET_KUIPER_PVAL_TRACE = False
        if PLOT_DET_KUIPER_PVAL_TRACE:
            Times = [t/(24.*60.*60.) for t in self.tile_times]
            plt.plot(Times, self.detection_kuiper_p_val_trace)
            ax = plt.gca()
            ax.set_yscale('log')
            plt.ylabel("Detection Kuiper p-value")
            plt.xlabel("Time [days]")
            plt.show()
            plt.grid()

    def RunSIXTE(self,TileJsonFile):
        # Figure out maximum times we can do in the time available
        # Will include time slew as an extra exit flag in the loop over tiles,
        # cutting out the last m tiles (m<n) once the slew time consumes the equivalent time for m tiles.
        # Define the time to merger
        json_file = TileJsonFile['LISA Data File']
        json_file = json_file.split("inference_data")[0] + "inference_param_files" + json_file.split("inference_data")[-1]
        with open(json_file, 'r') as f:
            input_params = json.load(f)
        f.close()
        TimeToMerger = input_params["waveform_params"]["DeltatL_cut"]
        n_times = int((-TimeToMerger-self.T_init)//self.T_lat)

        # Get tiles from json file
        with open(TileJsonFile, 'r') as f:
            TileDict = json.load(f)
        f.close()

        # Sort the dictionary to contain just the tiles we want
        Extra_keys = ['Tile Strat', 'overlap', 'LISA Data File']
        TileDict_reduced = {k: v for k, v in TileDict.items() if k not in Extra_keys and int(k)<=n_times}
        max_key = [int(k) for k in TileDict_reduced.keys()]
        max_key = max(max_key)
        if max_key<n_times:
            n_times = max_key

        # Now loop through the SIXTE commands
        for i in range(n_times):
            TileDict[str(i)]['beta']
            TileDict[str(i)]['lambda']
