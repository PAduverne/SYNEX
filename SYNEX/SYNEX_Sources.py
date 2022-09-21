import astropy.units as u
import numpy as np
import healpy as hp
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import Planck13, z_at_value # needed only to convert a given distance to redshift at initialization
import astropy.units as u
from astropy.time import Time
import warnings
import os
from datetime import date
import pathlib
import pickle
import json
# import statistics
from scipy import stats
from scipy.stats import norm
from scipy import integrate

import SYNEX.SYNEX_Utils as SYU
from SYNEX.SYNEX_Utils import SYNEX_PATH

import lal
import lalsimulation as lalsim

import lisabeta.pyconstants as pyconstants
import lisabeta.lisa.lisa_fisher as lisa_fisher
import lisabeta.lisa.lisa as lisa
import lisabeta.lisa.lisatools as lisatools

# mpi stuff
try:
    from mpi4py import MPI
except ModuleNotFoundError:
    MPI = None

################### Define base classes and add useful properties/function ###################

# This is where we will add supernova and other source classes when we get around to it.
# Can even inherit from the binary class again for BNS events? Can include
# PM portion too to create the options for BNSM and BNSPM?

class SMBH_Merger:
    """Base Class for frequency domain strains from Binary Black Holes.

    Parameters
    ----------
    M: float
        Total mass of the black hole binary (``m1`` +``m2``)
    q: float
        Mass ratio of the black hole binary (``m1`` /``m2``, ``m1`` <``m2``)
    z: float
        Redshift of the black hole binary
    chi1: float, optional, [-1., 1.]
        The dimensionless spin parameter abs(a/m) for black hole ``m1`` defaults to 0.0.
    chi2: float, optional, [-1., 1.]
        The dimensionless spin parameter abs(a/m) for black hole ``m2`` defaults to 0.0.
    Deltat: float
        Time shift (s)
    dist: float
        Luminosity distance (Mpc)
    inc: float, [0., np.pi]
        Inclination angle (rad)
    phi: float, [-np.pi, np.pi]
        Observer's azimuthal phase (rad)
    lambda: float, [-np.pi, np.pi]
        Source longitude (rad)
    beta: float, [-np.pi/2., np.pi/2.]
        Source latitude (rad)
    psi: float, [0., np.pi]
        Polarization angle (rad)

    Extra Parameters -- these are still from gwent. Need to change to parameters necessary for LISABETA
    to calculate the corresponding waveform.
    ----------------
    f_min: float, optional
        The lowest frequency in natural units (``Mf``, G=c=1) at which the BBH waveform is calculated
    f_max: float, optional
        The highest frequency in natural units (``Mf``, G=c=1) at which the BBH waveform is calculated
    nfreqs: int, optional
        The number of frequencies at which the BBH waveform is calculated
    instrument: object, optional
        If assigned, the optimal frequency (ie. most sensitive frequency) of the detector is used as
        the binary's GW frequency
    f_gw: float, Quantity, optional
        The binary's GW frequency if source is monochromatic
    approximant: str, optional
        the approximant used to calculate the frequency domain waveform of the source.
        Can either be the python implementation of IMRPhenomD (``'pyPhenomD'``, the default) given below,
        or a waveform modelled in LIGO's ``lalsuite`` 's ``lalsimulation`` package.
    lalsuite_kwargs: dict, optional
        More specific user-defined kwargs for the different ``lalsuite`` waveforms

    Notes
    -----
    IMRPhenomD waveforms calibrated for ``q`` = ``m1`` /``m2`` < 18
    lisabeta waveforms assume aligned spins ``chi1``, ``chi2`` = abs(a/m) <= 0.85 or if ``q`` =1 abs(a/m)<0.98
    """

    def __init__(self, **kwargs_in):
        """
        TO DO:
        1. Need to check if read in vals from H5 file are LFrame or SSBFrame,
            and then convert accordingly::
            params_SSBframe=lisatools.convert_Lframe_to_SSBframe(
                params_Lframe, t0=0., frozenLISA=False, LISAconst=pyresponse.LISAconstProposal)

            But both frames should be stored in h5 file... need to make sure we are reading in the right ones.

        2. Need to clean up innit function by importing the defaults in SYNEX_PTMC code and running checks
            for values of keys there. This will be much more efficient and allow better control over changes to
            default values. Seperate these into a seperate file with source and LISA defaults markated clearly.
        """
        # Set verbosity
        self.verbose=kwargs_in.pop("verbose") if "verbose" in kwargs_in else True

        # Make sure to handle case where we are in cluster so we don't write too many files and exceed disk quota
        # self.use_mpi=kwargs["use_mpi"] if "use_mpi" in kwargs else False
        if MPI is not None:
            MPI_size = MPI.COMM_WORLD.Get_size()
            MPI_rank = MPI.COMM_WORLD.Get_rank()
            comm = MPI.COMM_WORLD
            use_mpi=(MPI_size > 1)
        else:
            use_mpi=False
            MPI_rank=0
        self.PermissionToWrite=not use_mpi # MPI_rank==0 # This will not write skymap file since it is memory instensive

        # Check sav file paths incase we are loading on a new system eg cluster
        if "ExistentialFileName" in kwargs_in.keys():
            try:
                ExistentialPath="/".join(kwargs_in["ExistentialFileName"].split("/")[:-1])
                pathlib.Path(ExistentialPath).mkdir(parents=True, exist_ok=True)
            except:
                ExistentialPath=SYNEX_PATH+"/Saved_Source_Dicts/"
                kwargs_in["ExistentialFileName"]=ExistentialPath+kwargs_in["ExistentialFileName"].split("/")[-1] ## In case we are now on a cluster or something and older saved files' paths no longer work
                pathlib.Path(ExistentialPath).mkdir(parents=True, exist_ok=True)
        if "NewExistentialFileName" in kwargs_in.keys():
            try:
                ExistentialPath="/".join(kwargs_in["NewExistentialFileName"].split("/")[:-1])
                pathlib.Path(ExistentialPath).mkdir(parents=True, exist_ok=True)
            except:
                ExistentialPath=SYNEX_PATH+"/Saved_Source_Dicts/"
                kwargs_in["NewExistentialFileName"]=ExistentialPath+kwargs_in["NewExistentialFileName"].split("/")[-1] ## In case we are now on a cluster or something and older saved files' paths no longer work
                pathlib.Path(ExistentialPath).mkdir(parents=True, exist_ok=True)

        # Default assume class is not mutated from another saved class
        # "MUTATED" = key to force new savefile
        MUTATED=False

        if "MUTATED" in kwargs_in:
            # Check if we want to force savefile either way
            MUTATED=kwargs_in["MUTATED"]
            if self.verbose: print("Source mutation set to",MUTATED)
            del kwargs_in["MUTATED"]
            kwargs=kwargs_in
        elif "ExistentialFileName" in kwargs_in.keys() and os.path.isfile(kwargs_in["ExistentialFileName"]):
            # Need to check existential file names before this point...
            self.ExistentialFileName=kwargs_in["ExistentialFileName"]

            # Load saved dictionary
            with open(self.ExistentialFileName, 'rb') as f:
                kwargs = pickle.load(f)

            # Check if we will modify something
            if "NewExistentialFileName" in kwargs_in:
                # requesting a new filename so force mutation to save to this name instead
                MUTATED = True
            elif len(kwargs_in.keys())>1:
                # more than just 'ExistentialFileName' specified
                ValueCheck = [value!=kwargs[key] for key,value in kwargs_in.items() if key in kwargs]
                if any([ValueCheck]):
                    # Values are changed so warn user
                    if self.verbose: print("Loaded source dictionary will be mutated- mutation will be saved to a new file...")
                    MUTATED = True

            # Set loaded dict values to kwargs
            kwargs.update(kwargs_in)
        else:
            kwargs=kwargs_in

        # Set variables - Can condense this to setattr(class,val) later but
        # wanted to see all variables while developing as they might change.
        for key, value in kwargs.items():
            # Base source parameters
            if key == "m1":
                self.m1 = value
            elif key == "m2":
                self.m2 = value
            elif key == "M":
                self.M = value
            elif key == "q":
                self.q = value
            elif key == "z":
                self.z = value
                self.dist = cosmo.luminosity_distance(value).to("Mpc").value
            elif key == "chi1":
                self.chi1 = value
            elif key == "chi2":
                self.chi2 = value
            elif key == "Deltat":
                self.Deltat = value
            elif key == "lambda":
                self.lamda = value # = RA-pi = phi-pi
            elif key == "beta":
                self.beta = value # = DEC = pi/2-theta
            elif key == "inc":
                self.inc = value
            elif key == "psi":
                self.psi = value
            elif key == "dist":
                if hasattr(self,"dist"):
                    if self.verbose: print("Distance already defined by redshift, sticking with the redshift equivalent distance.")
                else:
                    self.dist = value
                    self.z = z_at_value(cosmo.luminosity_distance, value*u.Mpc).value
            elif key == "phi":
                self.phi = value
            # Source parameters for waveform generation
            elif key == "minf":
                self.minf = value
            elif key == "maxf":
                self.maxf = value
            elif key == "t0":
                self.t0 = value
            elif key == "timetomerger_max":
                self.timetomerger_max = value # In YEARS
            elif key=='gps_timetomerger_max':
                # gps time at which we *start* the waveform
                # Equivalent to the time of "timetomerger_max" BEFORE merger happens
                self.gps_timetomerger_max=Time('2033-02-02T00:00:00.00', format='isot', scale='utc').gps
            elif key == "fend":
                self.fend = value
            elif key == "phiref":
                self.phiref = value
            elif key == "fref_for_phiref":
                self.fref_for_phiref = value
            elif key == "force_phiref_fref":
                self.force_phiref_fref = value
            elif key == "toffset":
                self.toffset = value
            elif key == "modes":
                self.modes = value
            elif key == "acc":
                self.acc = value
            elif key == "approximant":
                self.approximant = value
            elif key == 'DeltatL_cut':
                self.DeltatL_cut = value
            elif key == 'Lframe':
                self.Lframe = value
            elif key=='lisabetaFile':
                JsonFileLocAndName,H5FileLocAndName=SYU.CompleteLisabetaDataAndJsonFileNames(value)
                self.H5File=H5FileLocAndName
                self.JsonFile=JsonFileLocAndName
            elif key=='H5File':
                self.H5File=value
            elif key=='JsonFile':
                self.JsonFile=value
            elif key=='sky_map':
                self.sky_map=value
            elif key=='ExistentialFileName':
                self.ExistentialFileName=value
            elif key=='do3D':
                self.do3D=value
            elif key=='EM_Flux_Data':
                self.EM_Flux_Data=value
            elif key=='CTR_Data':
                self.CTR_Data=value

        #########
        ##
        ## check mass params... can we do 'pytools.complete_mass_params(source_params)' intead? Not sure
        ## this has the right default behaviour we want...
        ##
        #########

        # Make sure m1>m2 and q>1
        if hasattr(self,"m1") and hasattr(self,"m2"):
            if self.m1<self.m2:
                if self.verbose: print("Redefining m1>m2")
                tmp = self.m1
                self.m1 = self.m2
                self.m2 = tmp
        if hasattr(self,"q"):
            if self.q<1.:
                if self.verbose: print("Redefining q>=1")
                self.q = 1./self.q

        # Default mass parameters if all four mass params are givien - include a warning
        if all([hasattr(self,"M"), hasattr(self,"q"), hasattr(self,"m1"), hasattr(self,"m2")]):
            if self.m1 + self.m2 != self.M or self.m1/self.m2 != self.q:
                if self.verbose: print("The total mass and/or mass ratio does not match the sum of m1 and m2. Redefining based on m1 and m2")
                self.M = self.m1 + self.m2
                self.q = self.m1/self.m2

        # Default mass parameters if three mass params are givien - include a warning
        if all([hasattr(self,"m1"), hasattr(self,"m2"), hasattr(self,"M"), not hasattr(self,"q")]) or all([hasattr(self,"m1"), hasattr(self,"m2"), not hasattr(self,"M"), hasattr(self,"q")]):
            if hasattr(self,"M") and self.m1 + self.m2 != self.M:
                if self.verbose: print("The total mass does not match the sum of m1 and m2. Redefining based on m1 and m2")
                self.M = self.m1 + self.m2
                self.q = self.m1/self.m2
            if hasattr(self,"q") and self.m1/self.m2 != self.q:
                if self.verbose: print("The mass ratio does not match the sum of m1 and m2. Redefining based on m1 and m2")
                self.M = self.m1 + self.m2
                self.q = self.m1/self.m2
        if all([hasattr(self,"M"), hasattr(self,"q"), hasattr(self,"m1"), not hasattr(self,"m2")]) or all([hasattr(self,"M"), hasattr(self,"q"), not hasattr(self,"m1"), hasattr(self,"m2")]):
            if hasattr(self,"m1") and self.m1 != self.M*(self.q/(1.+self.q)):
                if self.verbose: print("The primary mass m1 does not match the total mass and mass ratio. Redefining based on M and q")
                self.m1 = self.M*(self.q/(1.+self.q))
                self.m2 = self.M/(1.+self.q)
            if hasattr(self,"m2") and self.m2 != self.M/(1.+self.q):
                if self.verbose: print("The secondary mass m2 does not match the total mass and mass ratio. Redefining based on M and q")
                self.m1 = self.M*(self.q/(1.+self.q))
                self.m2 = self.M/(1.+self.q)

        # Default mass parameters if two mass params are givien
        if hasattr(self,"m1") and hasattr(self,"m2"):
            self.M = self.m1 + self.m2
            self.q = self.m1/self.m2
        elif hasattr(self,"M") and hasattr(self,"q"):
            self.m1 = self.M*(self.q/(1.+self.q))
            self.m2 = self.M/(1.+self.q)
        elif hasattr(self,"M") and hasattr(self,"m1"):
            self.m2 = self.M-self.m1
            self.q = self.m1/self.m2
        elif hasattr(self,"M") and hasattr(self,"m2"):
            self.m1 = self.M-self.m2
            self.q = self.m1/self.m2
        elif hasattr(self,"q") and hasattr(self,"m1"):
            self.m2 = self.m1/self.q
            self.M = self.m1 + self.m2
        elif hasattr(self,"q") and hasattr(self,"m2"):
            self.m1 = self.q*self.m2
            self.M = self.m1 + self.m2
        elif hasattr(self,"M"): # Set defaults if only one mass param is givien
            if self.verbose: print("Assuming q=1")
            self.q = 1.
            self.m1 = self.M/2.
            self.m2 = self.m1
        elif hasattr(self,"q"):
            if self.verbose: print("Assuming m1=1e6 M_sol")
            self.m1 = 1.e6
            self.m2 = self.m1/self.q
            self.M = self.m1+self.m2
        elif hasattr(self,"m1"):
            if self.verbose: print("Assuming q=1")
            self.q = 1.
            self.m2 = self.m1
            self.M = self.m1+self.m2
        elif hasattr(self,"m2"):
            if self.verbose: print("Assuming q=1")
            self.q = 1.
            self.m1 = self.m2
            self.M = self.m1+self.m2

        # Default mass parameters if no mass param is givien
        if all([not hasattr(self,"M"), not hasattr(self,"q"), not hasattr(self,"m1"), not hasattr(self,"m2")]):
            if self.verbose: print("Assuming q=1 and m1=1e6 M_sol")
            self.m1 = 1.e6
            self.m2 = self.m1
            self.M = self.m1 + self.m2
            self.q = 1.

        # Default base source parameters
        if not hasattr(self,"z"): # If this is true then distance is also missing. Set default to redshift 3.
            if self.verbose: print("Setting z=3")
            self.z = 3.
            self.dist = cosmo.luminosity_distance(self.z).to("Mpc").value
        if not hasattr(self,"chi1"):
            if self.verbose: print("Setting chi1=0")
            self.chi1 = 0.
        if not hasattr(self,"chi2"):
            if self.verbose: print("Setting chi2=0")
            self.chi2 = 0.
        if not hasattr(self,"lamda"):
            if self.verbose: print("Setting lambda=0")
            self.lamda = 0.
        if not hasattr(self,"beta"):
            if self.verbose: print("Setting beta=0")
            self.beta = 0.
        if not hasattr(self,"inc"):
            if self.verbose: print("Setting inc=0")
            self.inc= 0.
        if not hasattr(self,"psi"):
            if self.verbose: print("Setting psi=0")
            self.psi = 0.
        if not hasattr(self,"Deltat"):
            if self.verbose: print("Setting Deltat=0")
            self.Deltat = 0.
        if not hasattr(self,"phi"):
            if self.verbose: print("Setting phi=0")
            self.phi = 0.

        # Default source parameters for waveform generation
        if not hasattr(self,"minf"):
            self.minf = 1e-5     # Included
        if not hasattr(self,"maxf"):
            self.maxf = 0.5     # Included
        if not hasattr(self,"t0"):
            self.t0 = 0.0     # Included
        if not hasattr(self,"timetomerger_max"):
            self.timetomerger_max = 1.0     # Included
        if not hasattr(self,"fend"):
            self.fend = None     # Included
        if not hasattr(self,"phiref"):
            self.phiref = 0.0     # Included
        if not hasattr(self,"fref_for_phiref"):
            self.fref_for_phiref = 0.0     # Included
        if not hasattr(self,"tref"):
            self.tref = 0.0     # Included
        if not hasattr(self,"fref_for_tref"):
            self.fref_for_tref = 0.0     # Included
        if not hasattr(self,"force_phiref_fref"):
            self.force_phiref_fref = True     # Included
        if not hasattr(self,"toffset"):
            self.toffset = 0.0     # Included
        if not hasattr(self,"modes"):
            self.modes = None     # Included
        if not hasattr(self,"acc"):
            self.acc = 1e-4     # Included
        if not hasattr(self,"approximant"):
            self.approximant = "IMRPhenomHM"     # Included
        if not hasattr(self,"DeltatL_cut"):
            self.DeltatL_cut = None     # Included
        if not hasattr(self,"Lframe"):
            self.Lframe = True
        if not hasattr(self,"H5File"):
            self.H5File=None
        if not hasattr(self,"JsonFile"):
            self.JsonFile=None
        if not hasattr(self,"sky_map"):
            self.sky_map=None
        if not hasattr(self,"ExistentialFileName"):
            # Default name to include source 'name' variable?
            today = date.today()
            d = today.strftime("%d_%m_%Y")
            ExistentialFile = d + "_SourceDict.dat"
            ExistentialFile=SYNEX_PATH+"/Saved_Source_Dicts/"+ExistentialFile
            self.ExistentialFileName=ExistentialFile
        if not hasattr(self,"do3D"):
            self.do3D=False
        if not hasattr(self,"gps_timetomerger_max"):
            self.gps_timetomerger_max=Time('2033-01-01T00:00:00.00', format='isot', scale='utc').gps

        # Check h5 and Json files are congruent...
        if self.JsonFile and not self.H5File:
            JsonFileLocAndName,H5FileLocAndName=SYU.CompleteLisabetaDataAndJsonFileNames(self.JsonFile)
            if os.path.isfile(H5FileLocAndName):
                if self.verbose: print("Using similar H5 file found \n"+H5FileLocAndName)
                self.H5File=H5FileLocAndName
        if not self.JsonFile and self.H5File:
            JsonFileLocAndName,H5FileLocAndName=SYU.CompleteLisabetaDataAndJsonFileNames(self.H5File)
            if os.path.isfile(JsonFileLocAndName):
                if self.verbose: print("Using similar Json file found \n"+JsonFileLocAndName)
                self.JsonFile=JsonFileLocAndName
        if self.JsonFile and self.H5File:
            JsonFileLocAndName,H5FileLocAndName=SYU.CompleteLisabetaDataAndJsonFileNames(self.JsonFile)
            if os.path.isfile(JsonFileLocAndName) and os.path.isfile(H5FileLocAndName):
                if self.verbose: print("Using similar Json file found \n"+JsonFileLocAndName+"and similar H5 file found \n"+H5FileLocAndName)
                self.JsonFile=JsonFileLocAndName
                self.H5File=H5FileLocAndName

        # If we resurrected with mutation, keep a reference to where this class came from
        if MUTATED:
            self.MutatedFromSourceFile = self.ExistentialFileName
            if "NewExistentialFileName" in kwargs:
                # Take new filename if given
                self.ExistentialFileName = kwargs["NewExistentialFileName"]
            else:
                # No new filename given- create it. NOTE: file extension left ambiguous
                try:
                    # Does it already have an extension number? If so, start there...
                    ExistentialFileExt=self.ExistentialFileName.split("_")[-1] # e.g. '3.dat', '4.config', '10.json'
                    ExistentialFileExt=int(ExistentialFileExt.split(".")[0])
                except:
                    # If not, start at 1
                    ExistentialFileExt = 1
                    self.ExistentialFileName = ".".join(self.ExistentialFileName.split(".")[:-1]) + "_1." + self.ExistentialFileName.split(".")[-1]

                # Find the first version that doesn't exist yet...
                while os.path.isfile(self.ExistentialFileName):
                    ExistentialFileExt+=1
                    self.ExistentialFileName = "_".join(self.ExistentialFileName.split("_")[:-1]) + "_" + str(ExistentialFileExt) + "." + self.ExistentialFileName.split(".")[-1]
            if self.verbose: print("Successfully mutated source:", self.MutatedFromSourceFile)
            if self.verbose: print("New savefile for mutation:", self.ExistentialFileName)

        # Check that Skymap file path exists - in case of subdirectory organizational architectures...
        # Json and H5 file paths already checked when names are completed, Existential and NewExistential checks at top of init function
        if self.sky_map!=None:
            try:
                SkyMapPath="/".join(self.sky_map.split("/")[:-1])
                pathlib.Path(SkyMapPath).mkdir(parents=True, exist_ok=True)
            except:
                SkyMapPath=SYNEX_PATH+"/Skymap_files" ## In case we are now on a cluster or something and older saved files' paths no longer work
                self.sky_map=SkyMapPath+"/"+self.sky_map.split("/")[-1]
                pathlib.Path(SkyMapPath).mkdir(parents=True, exist_ok=True)

        # Now check if sky_map needs creating or reading -- adaptation for 3D case needed here...
        if MUTATED and self.H5File!=None and os.path.isfile(self.H5File):
            if self.sky_map==None:
                self.sky_map = self.H5File.split("inference_data")[0] + 'Skymap_files' + self.H5File.split("inference_data")[-1]
                self.sky_map = self.sky_map[:-3] + '.fits'
            self.CreateSkyMapStruct()
        elif self.sky_map!=None and os.path.isfile(self.sky_map):
            self.LoadSkymap()
        elif self.sky_map!=None and not os.path.isfile(self.sky_map) and self.H5File!=None and os.path.isfile(self.H5File):
            self.CreateSkyMapStruct()
        elif self.sky_map==None and self.H5File!=None and self.H5File!=None:
            # Default name
            self.sky_map = self.H5File.split("inference_data")[0] + 'Skymap_files' + self.H5File.split("inference_data")[-1]
            self.sky_map = self.sky_map[:-3] + '.fits'
            if os.path.isfile(self.sky_map):
                self.LoadSkymap()
            else:
                try:
                    self.CreateSkyMapStruct()
                except:
                    if self.verbose: print("No h5 data located- skipping skymap calculation.")
                    self.PostSkyArea = None
                    self.FisherSkyArea = None
        else:
            self.PostSkyArea = self.calculatePostSkyArea()
            self.FisherSkyArea = self.calculateFisherSkyArea()

        # Save it all to file!
        if self.PermissionToWrite: self.ExistentialCrisis()

    def LoadSkymap(self):
        """
        Read a skymap in -- Check that parameters match? Not sure it's needed...

        Should we take "self.do3D" as a flag to only read in probs? Or should we
        fix this based on what we find in the skymap file?
        """
        self.map_struct={}
        try:
            # 3D case
            healpix_data, header = hp.read_map(self.sky_map,field=(0,1,2,3),h=True,verbose=self.verbose) # field=(0,1,2,3,4,5,6,7),h=True,verbose=self.verbose)
            prob_data = healpix_data[0]
            # cumprob_data = healpix_data[1]
            # ipix_keep_data = healpix_data[2]
            # pixarea_data = healpix_data[3]
            # pixarea_deg2_data = healpix_data[4]
            distmu_data = healpix_data[1] # [5]
            distsigma_data = healpix_data[2] # [6]
            norm_data = healpix_data[3] # [7]

            self.map_struct["prob"] = prob_data
            # self.map_struct["cumprob"] = cumprob_data
            # self.map_struct["ipix_keep"] = ipix_keep_data
            # self.map_struct["pixarea"] = pixarea_data
            # self.map_struct["pixarea_deg2"] = pixarea_deg2_data
            self.map_struct["distmu"] = distmu_data # / params["DScale"]
            self.map_struct["distsigma"] = distsigma_data # / params["DScale"]
            self.map_struct["distnorm"] = norm_data
            self.do3D = True
        except:
            # 1D case
            healpix_data, header = hp.read_map(self.sky_map,field=(0),h=True,verbose=self.verbose) # field=(0,1,2,3,4),h=True,verbose=self.verbose)
            prob_data = healpix_data # [0]
            # cumprob_data = healpix_data[1]
            # ipix_keep_data = healpix_data[2]
            # pixarea_data = healpix_data[3]
            # pixarea_deg2_data = healpix_data[4]
            prob_data = prob_data/np.sum(prob_data)
            self.map_struct["prob"] = prob_data
            # self.map_struct["cumprob"] = cumprob_data
            # self.map_struct["ipix_keep"] = ipix_keep_data
            # self.map_struct["pixarea"] = pixarea_data
            # self.map_struct["pixarea_deg2"] = pixarea_deg2_data
            self.do3D = False

        # Pixel locations
        print("loaded map struct probability properties",type(self.map_struct["prob"]), np.shape(self.map_struct["prob"]))
        nside = hp.pixelfunc.get_nside(self.map_struct["prob"])
        npix = hp.nside2npix(nside)
        theta, phi = hp.pix2ang(nside, np.arange(npix))
        ra = np.rad2deg(phi)
        dec = np.rad2deg(0.5*np.pi - theta)
        self.map_struct["ra"] = ra
        self.map_struct["dec"] = dec

        # Get additional data required by GWEMOPT
        sort_idx = np.argsort(self.map_struct["prob"])[::-1]
        csm = np.empty(len(self.map_struct["prob"]))
        csm[sort_idx] = np.cumsum(self.map_struct["prob"][sort_idx])
        self.map_struct["cumprob"] = csm
        self.map_struct["ipix_keep"] = np.where(csm <= 1.)[0] # params["iterativeOverlap"])[0]
        pixarea = hp.nside2pixarea(nside)
        pixarea_deg2 = hp.nside2pixarea(nside, degrees=True)
        self.map_struct["pixarea"] = pixarea
        self.map_struct["pixarea_deg2"] = pixarea_deg2

        # Posterior data sky area :: in sq. deg
        self.PostSkyArea = self.calculatePostSkyArea()

        # Fisher sky area :: in sq. deg
        self.FisherSkyArea = self.calculateFisherSkyArea()

    def CreateSkyMapStruct(self,SkyMapFileName=None):
        # if SkyMapFileName given then preferentially save to this filename.
        # This is checked in SYU.WriteSkymapToFile just before writing to file.
        if SkyMapFileName!=None: self.sky_map = SkyMapFileName

        # Read in data from file
        [infer_params, inj_param_vals, _, _] = SYU.read_h5py_file(self.H5File)
        if np.size(infer_params["lambda"][0])>1:
            nsamples = len(infer_params["lambda"][0])
        else:
            nsamples = len(infer_params["lambda"])

        # True locations for gwewmopt in SSB frame
        with open(self.JsonFile) as f: data = json.load(f)
        if data["source_params"]["Lframe"]:
            # Can't take from "source_params_SSBframe" because this hsn't saved for some reason...
            # So we include explicit conversion here just in case.
            params_SSBframe=lisatools.convert_Lframe_to_SSBframe(data["source_params"], t0=data["waveform_params"]["t0"], frozenLISA=data["waveform_params"]["frozenLISA"])
            l=params_SSBframe["lambda"]
            b=params_SSBframe["beta"]
            d=params_SSBframe["dist"]
        else:
            l=inj_param_vals["source_params_SSBframe"]["lambda"][0]
            b=inj_param_vals["source_params_SSBframe"]["beta"][0]
            d=inj_param_vals["source_params_SSBframe"]["dist"][0]
        self.true_ra = np.rad2deg(l) if l>0. else np.rad2deg(2.*np.pi+l) ## In deg
        self.true_dec = np.rad2deg(b)                                 ## In deg
        self.true_lamdaSSB = l                                        ## In rad
        self.true_betaSSB = b                                         ## In rad
        self.true_distance = d                                        ## in Mpc

        # Convert posterior chains to SSB frame if in L frame (NB: inj_param_vals is an h5 dataset, not a dictionary)
        if data["run_params"]["sample_Lframe"]:
            params_Lframe = {}
            params_Lframe["Deltat"]=infer_params["Deltat"] if "Deltat" in infer_params else inj_param_vals["source_params_Lframe"]["Deltat"][0]
            params_Lframe["lambda"]=infer_params["lambda"] if "lambda" in infer_params else inj_param_vals["source_params_Lframe"]["lambda"][0]
            params_Lframe["beta"]=infer_params["beta"] if "beta" in infer_params else inj_param_vals["source_params_Lframe"]["beta"][0]
            params_Lframe["psi"]=infer_params["psi"] if "psi" in infer_params else inj_param_vals["source_params_Lframe"]["psi"][0]
            params_Lframe["Lframe"]=data["run_params"]["sample_Lframe"]
            params_SSBframe=lisatools.convert_Lframe_to_SSBframe(params_Lframe, t0=data["waveform_params"]["t0"], frozenLISA=data["waveform_params"]["frozenLISA"])
        else:
            params_SSBframe = {}
            params_SSBframe["Deltat"]=infer_params["Deltat"] if "Deltat" in infer_params else inj_param_vals["source_params_SSBframe"]["Deltat"][0]
            params_SSBframe["lambda"]=infer_params["lambda"] if "lambda" in infer_params else inj_param_vals["source_params_SSBframe"]["lambda"][0]
            params_SSBframe["beta"]=infer_params["beta"] if "beta" in infer_params else inj_param_vals["source_params_SSBframe"]["beta"][0]
            params_SSBframe["psi"]=infer_params["psi"] if "psi" in infer_params else inj_param_vals["source_params_SSBframe"]["psi"][0]

        # Convert post angles to theta,phi according to conventions used in hp.projplot... (foudn by trial and error)
        # import matplotlib.pyplot as plt
        # plt.scatter(params_SSBframe["lambda"],params_SSBframe["beta"],label='SSB')
        # plt.scatter(infer_params["lambda"],infer_params["beta"],label='L')
        # plt.scatter(self.true_lamdaSSB, self.true_betaSSB, label='True SSB')
        # plt.scatter(self.lamda, self.beta, label='True L')
        # plt.xlim([-np.pi, np.pi])
        # plt.ylim([-np.pi/2., np.pi/2.])
        # plt.legend()
        # plt.show()
        post_phis = -params_SSBframe["lambda"]
        post_thetas = np.pi/2-params_SSBframe["beta"]

        # Start with basic pixel resolution
        nside=32
        npix = hp.nside2npix(nside)
        post_pix = hp.ang2pix(nside,post_thetas,post_phis)

        # Pixel probability according to location -- np.bincount is faster than np.histogram.
        probs = np.bincount(post_pix,weights=np.array([1./nsamples]*nsamples), minlength=npix)

        # Check what average non-zero bin population is
        non_zero_counts = nsamples*probs[np.where(probs>0)[0]]
        _,minmax,mean,_,_,_=stats.describe(non_zero_counts)

        # Do it all again but with better nside resolution, taking account of tiny posterior islands by capping to a max nside=2048
        max_nside = 2048
        popCheck = 2**int((np.log10(10.*nside**2/mean)/np.log10(4.))//1) ### from npix=12*nside**2 and nside=2**integer
        if popCheck<max_nside and popCheck>8:
            nside=popCheck
            npix=hp.nside2npix(nside)
        elif popCheck>max_nside:
            nside = max_nside
            npix=hp.nside2npix(nside)
        post_pix = hp.ang2pix(nside,post_thetas,post_phis)
        probs = np.bincount(post_pix,weights=np.array([1./nsamples]*nsamples), minlength=npix)

        # Check everything went as planned
        # non_zero_counts = nsamples*probs[np.where(probs>0)[0]]
        # _,minmax,mean,_,_,_=stats.describe(non_zero_counts)
        # if self.verbose: print(stats.describe(non_zero_counts), hp.nside2pixarea(nside,degrees=True), mean//10, nside)

        # Get pixel locations
        pix_thetas, pix_phis = hp.pix2ang(nside, np.arange(npix))
        pix_ras = np.rad2deg(pix_phis) # Ra and Dec for pix are stored in map_struct in degrees NOT radians
        pix_decs = np.rad2deg(0.5*np.pi - pix_thetas)

        # Start the dictionary
        map_struct = {"prob":probs,
                      "ra":pix_ras,
                      "dec":pix_decs}

        # Get distance stats per pixel if we want them AND they were included in inference
        ######## Can we ask for redshift? I think if we infer redshift we also record dist so they are mutually exclusive...
        if self.do3D:
            if not "dist" in infer_params:
                # Check if 'dist' was included in inference
                if self.verbose: print("Distance inference data requested but not found in h5 file... Proceeding with do3D=False.")
                self.do3D=False
            else:
                # Included- bin posteriors according to pixels. Easiest to do with pandas
                import pandas as pd
                ZeroProbPixels=np.where(probs==0.)[0]
                FullPix=np.append(post_pix,ZeroProbPixels)
                FullDist=np.append(infer_params["dist"],np.array([0.]*len(ZeroProbPixels)))
                df = pd.DataFrame({"pixel":FullPix,"dist":FullDist})
                df.sort_values(by=["pixel"],inplace=True)
                distmu = df.groupby(["pixel"]).mean().to_numpy().reshape((npix,)) # reshape because there is a ghost axis from dataframe
                distsigma = df.groupby(["pixel"]).std().to_numpy().reshape((npix,)) # reshape because there is a ghost axis from dataframe
                map_struct["distmu"]=distmu
                map_struct["distsigma"]=distsigma
                # Calculate probability weights due to distance...
                # Gwemopt defines this as N^_i in equation 2 of https://arxiv.org/pdf/1603.07333.pdf -- compare eqn 2 in paper with line 553 of 'gwemopt.utils.py'
                map_struct["distnorm"]=[]
                r = np.linspace(0, 36000)
                for pixID in range(npix):
                    if distmu[pixID]==0. or distsigma[pixID]==0. or np.isnan(distsigma[pixID]) or np.isnan(distmu[pixID]):
                        p_i = 0.
                    else:
                        dp_dr = probs[pixID] * r**2 * norm(distmu[pixID], distsigma[pixID]).pdf(r)
                        p_i=np.trapz(dp_dr,r)
                    map_struct["distnorm"].append(p_i)
                map_struct["distnorm"]=map_struct["distnorm"]/np.sum(map_struct["distnorm"]) # returns ndarray even if we started with a list
                Zeros = np.where(map_struct["distnorm"]==0.)[0]
                map_struct["distnorm"][Zeros]=np.nan
                map_struct["distmu"][Zeros]=np.nan
                map_struct["distsigma"][Zeros]=np.nan

        # Get additional data required by GWEMOPT
        sort_idx = np.argsort(map_struct["prob"])[::-1]
        csm = np.empty(len(map_struct["prob"]))
        csm[sort_idx] = np.cumsum(map_struct["prob"][sort_idx])
        map_struct["cumprob"] = csm
        map_struct["ipix_keep"] = np.where(csm <= 1.)[0] # params["iterativeOverlap"])[0]
        pixarea = hp.nside2pixarea(nside)
        pixarea_deg2 = hp.nside2pixarea(nside, degrees=True)
        map_struct["pixarea"] = pixarea
        map_struct["pixarea_deg2"] = pixarea_deg2

        print("calculated map struct probability properties",type(map_struct["prob"]), np.shape(map_struct["prob"]))

        # create class attribute
        self.map_struct=map_struct

        # Posterior data sky area :: in sq. deg
        self.PostSkyArea = self.calculatePostSkyArea()

        # Fisher sky area :: in sq. deg
        self.FisherSkyArea = self.calculateFisherSkyArea()

        # Save to file
        SYU.WriteSkymapToFile(self.map_struct,self.sky_map,None,self.PermissionToWrite)

    def calculateFisherSkyArea(self,LISA=None,ConfLevel=0.9):
        """
            Explicit function to calculate and return sky area from Fisher matrix.

            NB: Need a LISA object !
        """
        if LISA!=None:
            # Get Fisher covariance matrix
            fishercov = SYU.GetFisher_smbh(self, LISA, **{})

            # Use lisatools in lisabeta for consistent sky area calculation
            FisherSkyArea = lisatools.sky_area_cov(fishercov, sq_deg=True, n_sigma=None, prob=ConfLevel)
        elif os.path.isfile(self.JsonFile):
            # Open json to read in useful dictonaries
            with open(self.JsonFile) as f: data = json.load(f)

            # Read dicts out
            param_dict = data["source_params"]
            waveform_params = data["waveform_params"]
            
            try:
                # Get fishercov from lisabeta functions
                fishercov = lisa_fisher.fisher_covariance_smbh(param_dict, **waveform_params)

                # Use lisatools in lisabeta for consistent sky area calculation
                FisherSkyArea = lisatools.sky_area_cov(fishercov, sq_deg=True, n_sigma=None, prob=ConfLevel)
            except:
                print("Lisabeta Fisher calls failed. If you asked for multimodal and DeltatL_cut, it might be because the cut is too close to the start freqs for each mode and one or more mode cannot initiate. Lisabeta does not contain a check for a non-zero mode signal, so we try again here with (2,2) only... Possibly will get errors when saving skymap. Think about lowering Mtot, for example, as this tends to be the only solution.")
                waveform_params["approximant"]="IMRPhenomD"
                fishercov = lisa_fisher.fisher_covariance_smbh(param_dict, **waveform_params)
                FisherSkyArea = lisatools.sky_area_cov(fishercov, sq_deg=True, n_sigma=None, prob=ConfLevel)
                print("Found area:",FisherSkyArea)
        else:
            FisherSkyArea = None

        return FisherSkyArea

    def calculatePostSkyArea(self,ConfLevel=0.9):
        """
            Explicit function to calculate and return sky area from self.map_struct.

            NB: Need a viable map_struct (eg from full lisabeta inference run)!
        """
        if hasattr(self,"map_struct"):
            # Find non-zero sky_map probabilities
            # NonZeros = np.nonzero(self.map_struct["prob"]>0)[0]
            NonZeros = np.unique(self.map_struct["ipix_keep"])

            # Get out all non-zero prob tiles for sky area calculation
            NonZero_Probs = self.map_struct["prob"][NonZeros]
            # NonZero_ipix_keep = self.map_struct["ipix_keep"][NoneZeros]
            # NonZero_ipix_keep = NonZeros

            # Sort all probabilities
            index = np.argsort(NonZero_Probs)
            index = [index[i] for i in range(len(index)-1,0,-1)] ## Need to be flipped and there will be a better way using np.argsort with a kwarg...
            allTiles_probs_sorted = NonZero_Probs[index]
            # allTiles_ipixs_sorted = NonZero_ipix_keep[index]

            if len(allTiles_probs_sorted)>0:
                cumsum_sorted = np.cumsum(allTiles_probs_sorted)
                CL_eq_len = np.argmax(cumsum_sorted>=ConfLevel)+1
                CL_under_len = CL_eq_len-1

                # Length of cumsum to equivalent area
                if CL_under_len != 0:
                    DecP = 1.-(cumsum_sorted[CL_eq_len-1] - ConfLevel)/(cumsum_sorted[CL_eq_len-1] - cumsum_sorted[CL_under_len-1])
                else:
                    DecP = 1.
                if cumsum_sorted[CL_eq_len-1]==ConfLevel:
                    CL_len = CL_eq_len
                else:
                    CL_len = CL_under_len+DecP

                # Sky area from posteriors
                PostSkyArea = CL_len*self.map_struct["pixarea_deg2"]
            else:
                PostSkyArea = 0.
        else:
            PostSkyArea = None

        return PostSkyArea

    def GenerateEMFlux(self,fstart22=1e-4,TYPE="const",**EM_kwargs):
        """
        Function to calculate EM flux of source.

        https://arxiv.org/pdf/1801.02266.pdf ???

        PARAMS:
        -------

        TYPE : string {"const","lalsim"}
            Type of EM calculation to do.
             -- "const" calculates a rough approx based on a typical
                binary with seperation 10 times sum of Schwartzchild radii
             -- "lalsim" uses lalsimulation to calculate a full time dependent
                waveform and use an approximation similar to idea paper for
                variable flux. NB: This method is HUGELY memory intensive so avoid
                if possible when on a cluster with restrictive memory limits (<several GB per source calculation)
        """

        if TYPE=="const":
            # Calculate flux - xray_flux formula is *DETECTOR* frame using *SOURCE* frame properties
            # NPoints=int(-self.DeltatL_cut//(0.5*(1.+self.z))) if self.DeltatL_cut!=None else 100000
            # xray_time=[(NPoints-ii)*(self.DeltatL_cut)/(NPoints*(1.+self.z)) for ii in range(NPoints)] # seconds to merger
            NPoints=int(-self.DeltatL_cut//0.5) if self.DeltatL_cut!=None else 100000
            xray_time=[(NPoints-ii)*self.DeltatL_cut/NPoints for ii in range(NPoints)] # seconds to merger
            t_start_flux=xray_time[0]-1.
            t_end_flux=xray_time[-1]+1.
            r_sch_1 = 2950.*self.m1 # Schwartzchild rad of primary in meters using m1 in SOLAR MASSES
            r_sch_2 = 2950.*self.m2 # Schwartzchild rad of primary in meters using m2 in SOLAR MASSES
            r=[(r_sch_1+r_sch_2)*10. for _ in range(NPoints)] ### 10 equivalent Sch. radii
            Period=[np.sqrt(4.*np.pi*np.pi*r_i**3/(pyconstants.G_SI*self.M*pyconstants.MSUN_SI)) for r_i in r] # Kepler approx
            Orb_freqs = [1./T for T in Period]
            Orb_Omega=[2.*np.pi/T for T in Period]
            v=[r_i*om_i for r_i,om_i in zip(r,Orb_Omega)] # in m/s
            Orb_phi=(integrate.cumtrapz(Orb_Omega, [-t for t in xray_time], initial=0.))
            Orb_phi_mod=self.psi-(Orb_phi[-1]%(2*np.pi)) # psi is GW polarization at merger in range [0,2*pi]
            Orb_phi=[(ph+Orb_phi_mod)%(2*np.pi) for ph in Orb_phi]
            BolFac = [0.1+(r_sch_1+r_sch_2)/r_i for r_i in r] # This should be normalized as per idea paper but is this a sum of radii or average or equivalent for EoB?
            mu_dopp = [3.*np.sin(self.inc)*np.cos(Orb_phi[ii])*v[ii]/pyconstants.C_SI+1. for ii in range(len(Orb_phi))] # Speed normalized to c=1
            if self.Lframe:
                M_tot = self.M/(1.+self.z)
            else:
                M_tot = self.M
            EddLum = (1.26e38)*M_tot*(1.+self.z) # /(918.*2.) # https://en.wikipedia.org/wiki/Eddington_luminosity -- total mass in solar mass units in source frame
            xray_flux = [mu_dopp[ii]*BolFac[ii]*EddLum/(4.*np.pi*self.dist*1e6*pyconstants.PC_SI*self.dist*1e6*pyconstants.PC_SI*10000.) for ii in range(len(BolFac))] # erg/s/cm^2
            GW_freqs = [2./T for T in Period]
            GW_Omega=[4.*np.pi/T for T in Period]
            GW_phi=[ph*2. for ph in Orb_phi]
        elif TYPE=="lalsim":
            # Empty LALParams dict
            LALParams = lal.CreateDict();
            modearray = lalsim.SimInspiralCreateModeArray()
            lalsim.SimInspiralModeArrayActivateMode(modearray, 2, 2);
            lalsim.SimInspiralModeArrayActivateMode(modearray, 2, 1);
            seobflags = lal.CreateDict();
            lal.DictInsertINT4Value(seobflags, "SEOBNRv4P_SpinAlignedEOBversion", 4);
            # Generate P-frame modes m<0 with the symmetry hP_l-m ~ (-1)^l hP_lm*
            lal.DictInsertINT4Value(seobflags, "SEOBNRv4P_SymmetrizehPlminusm", 1);
            # Use numerical or analytical derivatives of the Hamiltonian
            # Default is numerical with the flag 1
            NumericalOrAnalyticalHamiltonianDerivative = lalsim.SimInspiralWaveformParamsLookupEOBChooseNumOrAnalHamDer(LALParams);
            # Removed test of value of NumericalOrAnalyticalHamiltonianDerivative
            if (NumericalOrAnalyticalHamiltonianDerivative == lalsim.FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_NUMERICAL):
                lal.DictInsertINT4Value(seobflags, "SEOBNRv4P_HamiltonianDerivative", lalsim.FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_NUMERICAL)
            else:
                lal.DictInsertINT4Value(seobflags, "SEOBNRv4P_HamiltonianDerivative", NumericalOrAnalyticalHamiltonianDerivative)
            # Extension of Euler angles post-merger: simple precession around final J at a rate set by QNMs
            lal.DictInsertINT4Value(seobflags, "SEOBNRv4P_euler_extension", lalsim.FLAG_SEOBNRv4P_EULEREXT_QNM_SIMPLE_PRECESSION);
            # Z-axis of the radiation frame L
            lal.DictInsertINT4Value(seobflags, "SEOBNRv4P_Zframe", lalsim.FLAG_SEOBNRv4P_ZFRAME_L);
            # No debug output
            lal.DictInsertINT4Value(seobflags, "SEOBNRv4P_debug", 0);
            chi1_x=0.
            chi1_y=0.
            chi2_x=0.
            chi2_y=0. # assuming aligned spins with orbital momentum...
            # if self.verbose: print(type(phi_c), type(self.Deltat), type(self.m1*pyconstants.MSUN_SI), type(self.m2*pyconstants.MSUN_SI), type(fstart22), type(self.dist*1e6*pyconstants.PC_SI), type(self.inc), type(chi1_x), type(chi1_y), type(self.chi1), type(chi2_x), type(chi2_y), type(self.chi2), type(modearray), type(seobflags))
            # if self.verbose: print((phi_c), (self.Deltat), (self.m1*pyconstants.MSUN_SI), (self.m2*pyconstants.MSUN_SI), (fstart22))
            # if self.verbose: print((self.dist*1e6*pyconstants.PC_SI), (self.inc))
            # if self.verbose: print(chi1_x, chi1_y, self.chi1)
            # if self.verbose: print(chi2_x, chi2_y, self.chi2)
            # if self.verbose: print((modearray), (seobflags))
            # I assume it takes parameters in the source frame...?
            # if self.Lframe: # second input var (delta_t) was set to 1.5... can we increase this? Should 1st parameter be self.psi?
            #     if self.verbose: print("Lframe; rescaling to source fram using z =",self.z)
            #     hplus, hcross, hIlm, hJlm, dyn_Low, dyn_Hi, dyn_all, t_vec_modes, hP22_amp, hP22_phase, hP21_amp, hP21_phase, hP33_amp, hP33_phase, hP44_amp, hP44_phase, hP55_amp, hP55_phase, alphaJ2P, betaJ2P, gammaJ2P, AttachPars = lalsim.SimIMRSpinPrecEOBWaveformAll(self.phi, 1.5, self.m1*pyconstants.MSUN_SI/(1.+self.z), self.m2*pyconstants.MSUN_SI/(1.+self.z), fstart22/(1.+self.z), self.dist*1e6*pyconstants.PC_SI, self.inc, chi1_x, chi1_y, self.chi1, chi2_x, chi2_y , self.chi2, modearray, seobflags)
            # else:
            #     hplus, hcross, hIlm, hJlm, dyn_Low, dyn_Hi, dyn_all, t_vec_modes, hP22_amp, hP22_phase, hP21_amp, hP21_phase, hP33_amp, hP33_phase, hP44_amp, hP44_phase, hP55_amp, hP55_phase, alphaJ2P, betaJ2P, gammaJ2P, AttachPars = lalsim.SimIMRSpinPrecEOBWaveformAll(self.phi, 1.5, self.m1*pyconstants.MSUN_SI, self.m2*pyconstants.MSUN_SI, fstart22, self.dist*1e6*pyconstants.PC_SI, self.inc, chi1_x, chi1_y, self.chi1, chi2_x, chi2_y , self.chi2, modearray, seobflags)
            _, _, _, _, _, _, dyn_all, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = lalsim.SimIMRSpinPrecEOBWaveformAll(self.phi, 1.5, self.m1*pyconstants.MSUN_SI, self.m2*pyconstants.MSUN_SI, fstart22, self.dist*1e6*pyconstants.PC_SI, self.inc, chi1_x, chi1_y, self.chi1, chi2_x, chi2_y , self.chi2, modearray, seobflags)
            # The output here is a single REAL8Vector (I assume in the source frame...?), to be reshaped into the format:
            # t, x, y, z, px, py, pz, s1x, s1y, s1z, s2x, s2y, s2z, phiD, phi, vx, vy, vz, r, phi, pr, pphi, Omega, s1dotZ, s2dotZ, H
            # Spins in units of M^2
            nptdyn = int(len(dyn_all.data) / 26)
            dyn_arr = np.reshape(dyn_all.data, (26, nptdyn))

            # Take out the relevant parameters
            t = dyn_arr[0]
            x = dyn_arr[1]
            y = dyn_arr[2]
            z = dyn_arr[3]
            px = dyn_arr[4]
            py = dyn_arr[5]
            pz = dyn_arr[6]
            s1x = dyn_arr[7]
            s1y = dyn_arr[8]
            s1z = dyn_arr[9]
            s2x = dyn_arr[10]
            s2y = dyn_arr[11]
            s2z = dyn_arr[12]
            r = dyn_arr[18]
            phi = dyn_arr[19]
            pr = dyn_arr[20]
            pphi = dyn_arr[21]
            Omega = dyn_arr[22]
            H = dyn_arr[25]

            # Binary orbit properties
            t = [time for time in t] # *(1.+self.z) for time in t] # Assume time is in source frame and need to convert to detector frame?!
            dt = np.diff(t)
            dr = np.diff(r)
            r = r.tolist()
            r_dot = [(dr[ii]/dt[ii])/2. for ii in range(len(dr))]
            x_dot = px.tolist()
            y_dot = py.tolist()
            z_dot = pz.tolist()

            # GW properties are in *DETECTOR* frame
            GW_Omega = [Om*(1.+self.z) for Om in Omega] # [Om*2.*(1.+self.z) for Om in Omega]
            GW_freqs = [Om/(2.*np.pi) for Om in Omega]
            GW_phi = [phi_i for phi_i in phi] # [phi_i*2. for phi_i in phi]

            # Orbital velocity
            v = [np.sqrt(x_dot[ii]**2 + y_dot[ii]**2 + z_dot[ii]**2) for ii in range(len(x_dot))] # [np.sqrt(x_dot[ii]**2 + y_dot[ii]**2 + z_dot[ii]**2) for ii in range(len(x_dot))]
            v_vec = [[x_dot[ii],y_dot[ii],z_dot[ii]] for ii in range(len(x_dot))] # nx3 list

            # take out the largest component in v_vec
            v_abs = [np.abs(v_comp) for v_comp in v_vec]
            sign_comp = [np.sign(v_comp) for v_comp in v_vec]
            v_max = np.max(v_abs,axis=1)
            v_max = [v_max[ii]*np.sign(v_vec[ii][np.where(v_abs[ii]==v_max[ii])[0][0]]) for ii in range(np.shape(v_vec)[0])]
            import copy
            v_nomax = copy.deepcopy(v_vec)
            for ii in range(len(v_max)):
                v_nomax[ii].remove(v_max[ii])

            # set the sign of total speed depending on what the components of velocity are doing
            tmp = [np.sign(v_comp)==-1. for v_comp in v_vec]
            for ii in range(len(v_max)):
                if all([np.sign(v_comp)==-1. for v_comp in v_vec[ii]]):
                    v[ii] = -v[ii]
                elif all([np.sign(v_comp)==-1. for v_comp in v_nomax[ii]]) and np.sign(v_max[ii])==1. and (np.mean(v_nomax[ii])>np.abs(v_max[ii])):
                    v[ii] = -v[ii]
                elif all([np.sign(v_comp)==1. for v_comp in v_nomax[ii]]) and np.sign(v_max[ii])==-1. and (np.mean(v_nomax[ii])<np.abs(v_max[ii])):
                    v[ii] = -v[ii]

            # Calculate flux - xray_flux formula is *DETECTOR* frame using *SOURCE* frame properties
            BolFac = [0.1+1./r_i for r_i in r]
            # mu_dopp = [3.*np.sin(self.inc)*np.cos(GW_phi[ii]/2.)*v[ii]+1. for ii in range(len(GW_phi))] # phi and v here are *orbital* phase and velocity
            mu_dopp = [3.*np.sin(self.inc)*np.cos(GW_phi[ii]/2.)*v[ii]+1. for ii in range(len(GW_phi))]
            if self.Lframe:
                M_tot = self.M/(1.+self.z)
            else:
                M_tot = self.M
            EddLum = (1.26e38)*M_tot*(1.+self.z) # /(918.*2.) # https://en.wikipedia.org/wiki/Eddington_luminosity -- total mass in solar mass units
            xray_time = [(t[ii]-t[-1])*(M_tot*pyconstants.MSUN_SI*(pyconstants.G_SI/(pyconstants.C_SI**3))) for ii in range(len(t))] # corrected for diff's taken for velocity To be really accurate we should also interpolate phi and r to the right times...
            xray_flux = [mu_dopp[ii]*BolFac[ii]*EddLum/(4.*np.pi*self.dist*1e6*pyconstants.PC_SI*self.dist*1e6*pyconstants.PC_SI*10000.) for ii in range(len(BolFac))] # erg/s/cm^2 - this is orbital properties in *SOURCE* frame, then scaled by lum dist to account for redshifting to *DETECTOR* frame

            # Start time - LISABETA takes the Lframe flag appropriately. Need to convert to *SOURCE* frame
            param_dict,waveform_params,extra_params = SYU.ClassesToParams(self,detector=None,**EM_kwargs) ## Do we need to change the approximant here from the saved one in source to EOB one?
            waveform_params["DeltatL_cut"] = None # Need to set this to 0 to get the remaining signal for photon phases
            # if self.verbose: print("Params from 'ClassesToParams':",param_dict,waveform_params,extra_params)
            if param_dict["m1"]==param_dict["m2"]:
                param_dict["m1"] = param_dict["m1"]*(1.+1e-6) # lisabeta needs q>1, so just adjust lightly so that m1>m2...
            waveform_params["modes"]=[(2,2)] # cumulative SNR functions work for 22 mode only right now
            cumulsnr = lisa.CumulSNRLISATDI_SMBH(param_dict, **waveform_params)
            tf = cumulsnr['tf']
            cumul_SNR = cumulsnr['SNRcumul']
            # t_start_flux = -7.5*24.*60.*60. #
            t_start_flux = [time for (time,SNR) in zip(tf,cumul_SNR) if SNR<10][-1] # Grab the last time before we reach SNR detection threshold - this is negative seconds to merger
            # Return is always Lframe, now decide if this is larger or smaller than the cut time
            t_start_flux/=(1.+self.z) # This is slightly off - can you look up their relation between time and frequency they quote from the paper? Bottom of page 3
            if self.DeltatL_cut != None:
                if self.verbose: print(r"T$_{SNR=10}$",(1.+self.z)*t_start_flux/(24.*60.*60.),r"t$_{cut}$", self.DeltatL_cut/(24.*60.*60.))
                if (1.+self.z)*t_start_flux<self.DeltatL_cut:
                    t_start_flux = self.DeltatL_cut # Because we don't start at a time that is before the sky area is provided
                # Here we can either code the latency time to be the period of the GW wave or we can start when the period is comparable with the latency time requested...

            # End Time - *SOURCE* frame
            ISCO_freq = 4400.*(1.+self.z)/self.M # 4400Hz is *DETECTOR* frame GW ISCO frequency
            if self.Lframe:
                ISCO_freq*=(1.+self.z) # extra factor for total mass in Lframe
            t_end_flux = [time for (time,GW_freq) in zip(xray_time,GW_freqs) if GW_freq<=ISCO_freq]
            if isinstance(t_end_flux,list) and len(t_end_flux)>0:
                t_end_flux = t_end_flux[-1]
            elif isinstance(t_end_flux,list):
                t_end_flux=xray_time[-1]

        # Set the class variables using the start and end times - convert everything *EXCEPT* xray_flux to *DETECTOR* frame
        EM_Flux_Data={}
        EM_Flux_Data["EM_Approx"] = TYPE
        EM_Flux_Data["xray_flux"] = [flux for (flux,time) in zip(xray_flux,xray_time) if time>t_start_flux and time<t_end_flux] # if time>t_start_flux and time<t_end_flux else 0. for (flux,time) in zip(xray_flux,xray_time)]
        EM_Flux_Data["xray_time"] = [time*(1.+self.z) for time in xray_time if time>t_start_flux and time<t_end_flux]
        EM_Flux_Data["GW_phi"] = [phi for (phi,time) in zip(GW_phi,xray_time) if time>t_start_flux and time<t_end_flux] # GW_phi
        EM_Flux_Data["GW_Omega"] = [Om/(1.+self.z) for (Om,time) in zip(GW_Omega,xray_time) if time>t_start_flux and time<t_end_flux] # [Om/(1.+self.z) for Om in GW_Omega]
        EM_Flux_Data["r"] = [rad for (rad,time) in zip(r,xray_time) if time>t_start_flux and time<t_end_flux] # r
        EM_Flux_Data["GW_freqs"] = [f/(1.+self.z) for (f,time) in zip(GW_freqs,xray_time) if time>t_start_flux and time<t_end_flux] # [f/(1.+self.z) for f in GW_freqs]
        self.EM_Flux_Data = EM_Flux_Data

        # Save it
        self.ExistentialCrisis()

    def GenerateCTR(self,ARF_file_loc_name,gamma=1.7):
        CTR_Data={}
        CTR_Data["xray_gamma"] = gamma
        if not hasattr(self,"EM_Flux_Data"):
            raise ValueError("No EM flux generated. Need to call GenerateEMFlux with a detector object before running this function.")

        # CTR_Data["xray_phi_0"] = [(6.242e8)*self.EM_Flux_Data["xray_flux"][ii]*(2.-gamma)/(10.**(2.-gamma)-0.2**(2.-gamma)) for ii in range(len(self.EM_Flux_Data["xray_flux"]))]
        CTR_Data["xray_phi_0"] = [(6.242e8)*flu*(2.-gamma)/(10.**(2.-gamma)-0.2**(2.-gamma)) for flu in self.EM_Flux_Data["xray_flux"]]
        from astropy.io import fits
        hdul = fits.open(ARF_file_loc_name)
        # hdul.info()
        ARF = hdul[1].data[:]
        N = len(ARF)
        E_bins = [ARF[ii][0] for ii in range(N)] # Units are keV
        E_bins.append(ARF[-1][1]) # So that this length is N+1 where N is the length of ARF_func
        ARF_func = [ARF[ii][2] for ii in range(N)]
        dE_bins = np.diff(E_bins)

        # Now compute the CTR - photons s^-1
        integrand = [dE_bins[ii]*ARF_func[ii]*((E_bins[ii]+E_bins[ii+1])/2.)**(-1.*gamma) for ii in range(N)]
        integral = sum(integrand)/(1.+self.z) # We need another 1/(1+z) here otherwise the photon rate is too high. Maybe this is a redshifting of energy bins?
        CTR_Data["CTR"] = [integral*el for el in CTR_Data["xray_phi_0"]] # CTR
        # if self.verbose: print("Integral checks:",np.mean(integrand),integral,1.+self.z)
        # import matplotlib.pyplot as plt
        # EBinsPlot=[((E_bins[ii]+E_bins[ii+1])/2.) for ii in range(N)]
        # plt.loglog(EBinsPlot,integrand,label=r"Convolution")
        # plt.xlabel(r"Energy [keV]",fontsize="xx-small")
        # plt.legend(fontsize="xx-small")
        # plt.show()
        # plt.loglog(EBinsPlot,ARF_func,label=r"ARF")
        # plt.loglog(EBinsPlot,[((E_bins[ii]+E_bins[ii+1])/2.)**(-1.*gamma) for ii in range(N)],label=r"Energy relationship")
        # plt.xlabel(r"Energy [keV]",fontsize="xx-small")
        # plt.legend(fontsize="xx-small")
        # plt.show()
        # plt.semilogy([t/(24.*60.*60.) for t in self.EM_Flux_Data["xray_time"]],CTR_Data["xray_phi_0"],label=r"X-ray $\phi_0$")
        # plt.xlabel(r"Time to Merger [d]",fontsize="xx-small")
        # plt.legend(fontsize="xx-small")
        # plt.show()
        # plt.semilogy([t/(24.*60.*60.) for t in self.EM_Flux_Data["xray_time"]],CTR_Data["CTR"],label=r"CTR")
        # plt.xlabel(r"Time to Merger [d]",fontsize="xx-small")
        # plt.legend(fontsize="xx-small")
        # plt.show()

        # Store data in source object
        self.CTR_Data=CTR_Data

        # Save it
        self.ExistentialCrisis()

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
        # Gather attributes to dict and take care of "lamda" -> "lambda"
        MyExistentialDict = dict(self.__dict__)
        MyExistentialDict["lambda"]=self.lamda
        del MyExistentialDict["lamda"]

        if "map_struct" in MyExistentialDict: del MyExistentialDict["map_struct"]
        if "CTR_Data" in MyExistentialDict: del MyExistentialDict["CTR_Data"]
        if "EM_Flux_Data" in MyExistentialDict: del MyExistentialDict["EM_Flux_Data"]

        # Save to file...
        if self.verbose: print("Saving source attributes to:",self.ExistentialFileName)
        with open(self.ExistentialFileName, 'wb') as f:
            pickle.dump(MyExistentialDict, f)
        if self.verbose: print("Done.")
