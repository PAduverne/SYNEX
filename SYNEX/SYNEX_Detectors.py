import astropy.units as u
import astropy.stats as astat
import astropy, astroplan
import numpy as np
import SYNEX.SYNEX_Utils as SYU
from SYNEX.SYNEX_Utils import SYNEX_PATH
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
    detector_go_params : (gwemopt compatible) dict
        See './SYNEX/gwemopt_defaults.py' all possible flags.

    detector_config_struct : (gwemopt compatible) dict
        See './SYNEX/gwemopt_defaults.py' all possible flags.

    ExistentialFileName : Path or string
        Path and file name of file loadable by 'pickle'. Contents are a
        dictionary containing all the kwargs of a saved 'Athena' class object.
        If ExistentialFileName is specified but the file oesn't exist, then it
        will be automatically created so that the class can be loaded again later.
        Automatic save points through SYNEX routines will update this savefile

    NB: GWEMOPT uses serialized dictionaries 'go_params' and 'config_struct'.
        These are constructed with one or more 'detector_go_params' or
        'detector_config_struct' under key values given by the 'telescope' name.
        You might experience errors if you try passing 'detector_go_params' or
        'detector_config_struct' directly to GWEMOPT. You will need to pass classes
        through interface routines like SYNEX_Utils.InitGWEMOPTParamsFromDetectors.py
        to prepare the dictionaries for use in GWEMOPT.
    """

    def __init__(self, **kwargs):
        # Import default gwemopt dict
        from SYNEX.gwemopt_defaults import go_params_default, config_struct_default as detector_go_params, detector_config_struct

        # Default is not to recompute tesselation if the saved '.tess' file exists
        MUTATED=False

        # Check if we are resurrecting a class from a save file
        if "ExistentialFileName" in kwargs.keys():
            # FUSE THIS WITH CONFIG KEY IN GO_PARAMS... Save in two files so we
            # respect gwemopt conventions? Seems overkill but 'gwemop.utils.readParamsFromFile(file)'
            # already exists and can be used for both 'config' and 'params' dictionaries.
            self.ExistentialFileName=kwargs["ExistentialFileName"]

            # Load saved dictionary
            with open(self.ExistentialFileName, 'r') as f:
                SavedDict = pickle.load(f)

            # Use these values as default to modify with remaining keys in '**kwargs'
            detector_go_params = SavedDict["detector_go_params"]
            detector_config_struct = SavedDict["detector_config_struct"]

            # Check if we will modify something later
            if len(kwargs.keys())>1:
                # more than just 'ExistentialFileName' specified
                ValueCheck = [value!=SavedDict[key] for key,value in kwargs.items()]
                if any([ValueCheck]):
                    # Values are changed so recompute tesselation later even if the '.tess' file already exists
                    MUTATED=True

        # Change all things in go_params that are specified in kwargs,
        # and warn the user if any of the things they set are not used...
        print_reminder = False
        for key,value in kwargs.items():
            if key in detector_go_params:
                detector_go_params[key]=value
            else:
                print("'",key,"' not contained in gwemopt 'params' dict...")
                print_reminder = True
        # Now initiate config struct for this telescope
        if "ConfigFileName" in kwargs.keys():
            detector_config_struct = gwemopt.utils.readParamsFromFile(kwargs["ConfigFileName"])
            detector_config_struct["telescope"] = detector_go_params["telescopes"]
        else:
            for key,value in kwargs.items():
                if key in detector_config_struct:
                    detector_config_struct[key]=value
                else:
                    print("'",key,"' not contained in gwemopt 'config_struct' dict...")
                    print_reminder = True

        # Check if tesselation file was set
        if detector_config_struct["tesselationFile"]==None:
            detector_config_struct["tesselationFile"]=SYNEX_PATH+"/gwemopt_tess_files/" + detector_go_params["telescopes"] + ".tess"

        # Check if tesselation file exists and if we want to recompute - otherwise if it exists it will be loaded
        if MUTATED and os.path.isfile(self.detector_config_struct["tesselationFile"]):
            if "NewtesselationFile" in kwargs:
                # Take new filename if given
                self.detector_config_struct["tesselationFile"] = kwargs["NewtesselationFile"]
            else:
                # No new filename given- create it. NOTE: file extension is checked in 'SYU.GWEMOPTPathChecks()' step below...
                try:
                    # Does it already have an extension number? If so, start there...
                    TessFileExt=self.detector_config_struct["tesselationFile"].split("_")[-1] # e.g. '3.tess'
                    TessFileExt = int(TessFileExt.split(".")[0])
                except:
                    # If not,start at 1
                    TessFileExt = 1
                    self.detector_config_struct["tesselationFile"] = ".".join(self.detector_config_struct["tesselationFile"].split(".")[:-1]) + "_1." + self.detector_config_struct["tesselationFile"].split(".")[-1]

                # Find the first version that doesn't exist yet...
                while os.path.isfile(self.detector_config_struct["tesselationFile"]):
                    TessFileExt+=1
                    self.detector_config_struct["tesselationFile"] = "_".join(self.detector_config_struct["tesselationFile"].split("_")[:-1]) + "_" + str(TessFileExt) + "." + self.detector_config_struct["tesselationFile"].split(".")[-1]

        # Set config struct telescope name - for now there is only one 'telescope' in 'telescopes'
        detector_config_struct["telescope"]=detector_go_params["telescopes"]

        # Check that file names are all coherent with SYNEX_PATH and telescope name
        detector_go_params, detector_config_struct = SYU.GWEMOPTPathChecks(detector_go_params,detector_config_struct)

        # Set as class attributes
        self.detector_go_params = detector_go_params
        self.detector_config_struct = detector_config_struct

        # Set save file name if not already there
        if not hasattr(self,"ExistentialFileName"):
            # Default name to include telecope 'name' variable?
            today = date.today()
            d = today.strftime("%d_%m_%Y")
            ExistentialFile = d + "_TelescopeDict.dat"
            ExistentialFile=SYNEX_PATH+"/Saved_Telescope_Dicts/"+ExistentialFile
            self.ExistentialFileName=ExistentialFile

        # If we resurrected with mutation, keep a reference to where this class came from
        if MUTATED:
            self.MutatedFromDetectorFile = self.ExistentialFileName
            if "NewExistentialFileName" in kwargs:
                # Take new filename if given
                self.ExistentialFileName = kwargs["NewExistentialFileName"]
            else:
                # No new filename given- create it. NOTE: file extension left ambiguous
                try:
                    # Does it already have an extension number? If so, start there...
                    ExistentialFileExt=self.ExistentialFileName.split("_")[-1] # e.g. '3.dat' for '4.config'
                    ExistentialFileExt = int(ExistentialFileExt.split(".")[0])
                except:
                    # If not, start at 1
                    ExistentialFileExt = 1
                    self.ExistentialFileName = ".".join(self.ExistentialFileName.split(".")[:-1]) + "_1." + self.ExistentialFileName.split(".")[-1]

                # Find the first version that doesn't exist yet...
                while os.path.isfile(self.ExistentialFileName):
                    ExistentialFileExt+=1
                    self.ExistentialFileName = "_".join(self.ExistentialFileName.split("_")[:-1]) + "_" + str(ExistentialFileExt) + "." + self.ExistentialFileName.split(".")[-1]
            print("Successfully mutated detector:", self.MutatedFromSourceFile)
            print("New savefile for mutation:", self.ExistentialFileName)

        # Check that file paths exist - in case of subdirectory organizational architectures...
        # Tesselation path already checked in 'SYU.GWEMOPTPathChecks()'
        ExistentialPath="/".join(self.ExistentialFileName.split("/")[:-1])
        pathlib.Path(ExistentialPath).mkdir(parents=True, exist_ok=True)

        # Calculate tesselation - detector is saved at the end of this in case we recompute in other codes
        # If 'MUTATED' was False then if the '.tess' file exists the tesselation will be loaded.
        self.ComputeTesselation()

        # Issue reminder of where to find list of gwemopt variables and flags
        if print_reminder:
            print("Some keys given at initiatiation of Athena class are not contained in gwemop params - see 'SYNEX/gwemopt_defaults.py' for full list of possible field names.")

    def ComputeTesselation(self):
        if self.detector_go_params["doSingleExposure"]:
            exposuretime = np.array(self.detector_go_params["exposuretimes"].split(","),dtype=np.float)[0]
            nmag = -2.5*np.log10(np.sqrt(self.detector_config_struct["exposuretime"]/exposuretime))
            self.detector_config_struct["magnitude"] = self.detector_config_struct["magnitude"] + nmag
            self.detector_config_struct["exposuretime"] = exposuretime
        if "tesselationFile" in self.detector_config_struct:
            if not os.path.isfile(self.detector_config_struct["tesselationFile"]):
                if self.detector_config_struct["FOV_type"] == "circle":
                    gwemopt.tiles.tesselation_spiral(self.detector_config_struct)
                elif self.detector_config_struct["FOV_type"] == "square":
                    gwemopt.tiles.tesselation_packing(self.detector_config_struct)
            if self.detector_go_params["tilesType"] == "galaxy":
                self.detector_config_struct["tesselation"] = np.empty((3,))
            else:
                self.detector_config_struct["tesselation"] = np.loadtxt(self.detector_config_struct["tesselationFile"],usecols=(0,1,2),comments='%')

        if "referenceFile" in self.detector_config_struct: ### Not sure what this is but we include it to be complete with GWEMOPT
            from astropy import table
            refs = table.unique(table.Table.read(
                self.detector_config_struct["referenceFile"],
                format='ascii', data_start=2, data_end=-1)['field', 'fid'])
            reference_images =\
                {group[0]['field']: group['fid'].astype(int).tolist()
                for group in refs.group_by('field').groups}
            reference_images_map = {1: 'g', 2: 'r', 3: 'i'}
            for key in reference_images:
                reference_images[key] = [reference_images_map.get(n, n)
                                         for n in reference_images[key]]
            self.detector_config_struct["reference_images"] = reference_images

        location = astropy.coordinates.EarthLocation(self.detector_config_struct["longitude"],self.detector_config_struct["latitude"],self.detector_config_struct["elevation"])
        observer = astroplan.Observer(location=location)
        self.detector_config_struct["observer"] = observer

        # Save it all to file!
        self.ExistentialCrisis()

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

    def ExistentialCrisis(self,NewFileName=None):
        """
        Function to save all class attributes as a dictionary to file,
        making sure to overwrite existing files by the same name. This will
        make source resurrection easier if we do a long analysis run in stages.

        NB: When we get to tiling attributes may not be serializable for json files,
        so here we opt for pickling to '.dat' files instead.
        However, we do not check that the new 'FileName' has the right extension
        or path. Need to do this later.
        """
        if NewFileName!=None:
            # Check new filepath exists...
            NewFilePath="/".join(NewFileName.split("/")[:-1])
            pathlib.Path(NewFilePath).mkdir(parents=True, exist_ok=True)
            # Reset name in class attributes
            self.ExistentialFileName = NewFileName
        # Gather attributes to dict
        MyExistentialDict = self.__dict__
        # Save to file...
        print("Saving detector attributes...")
        with open(self.ExistentialFileName, 'wb') as f:
            pickle.dump(MyExistentialDict, f)
        print("Done.")
