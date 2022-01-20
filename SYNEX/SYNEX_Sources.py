import astropy.units as u
import numpy as np
import healpy as hp
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import Planck13, z_at_value # needed only to convert a given distance to redshift at initialization
import astropy.units as u
import warnings
import os

import SYNEX.SYNEX_Utils as SYU

import lal
import lalsimulation as lalsim

import lisabeta.pyconstants as pyconstants
import lisabeta.lisa.lisa as lisa

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

    def __init__(self, **kwargs):
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
                self.dist = cosmo.luminosity_distance(self.z).to("Mpc").value
            elif key == "chi1":
                self.chi1 = value
            elif key == "chi2":
                self.chi2 = value
            elif key == "Deltat":
                self.Deltat = value
            elif key == "lambda":
                self.lamda = value # // RA
            elif key == "beta":
                self.beta = value # // DEC
            elif key == "inc":
                self.inc = value
            elif key == "psi":
                self.psi = value
            elif key == "dist":
                if hasattr(self,"dist"):
                    print("Distance already defined by redshift, sticking with the redshift equivalent distance.")
                else:
                    self.dist = value
                    self.z = z_at_value(cosmo.luminosity_distance, value*u.Mpc)
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
                self.timetomerger_max = value
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
                if os.path.isfile(JsonFileLocAndName) and os.path.isfile(H5FileLocAndName):
                    self.H5File=H5FileLocAndName
                    self.JsonFile=H5FileLocAndName
                elif not os.path.isfile(JsonFileLocAndName) and os.path.isfile(H5FileLocAndName):
                    print("Warning: no Json file found- setting json and h5 filenames to None...")
                    self.H5File=None
                    self.JsonFile=None
                elif os.path.isfile(JsonFileLocAndName) and not os.path.isfile(H5FileLocAndName):
                    print("Warning: no h5 data file found- setting json and h5 filenames to None...")
                    self.H5File=None
                    self.JsonFile=None
                else:
                    print("Warning: no json or h5 data file found- setting json and h5 filenames to None...")
                    self.H5File=None
                    self.JsonFile=None
            elif key=='sky_map':
                self.sky_map=value

        # Make sure m1>m2 and q>1
        if hasattr(self,"m1") and hasattr(self,"m2"):
                if self.m1<self.m2:
                    print("Redefining m1>m2")
                    tmp = self.m1
                    self.m1 = self.m2
                    self.m2 = tmp
        if hasattr(self,"q"):
                if self.q<1.:
                    print("Redefining q>=1")
                    self.q = 1./self.q

        # Default mass parameters if all four mass params are givien - include a warning
        if all([hasattr(self,"M"), hasattr(self,"q"), hasattr(self,"m1"), hasattr(self,"m2")]):
            print("You have specified all mass parameters- this is uneeded but I will continue.")
            if self.m1 + self.m2 != self.M or self.m1/self.m2 != self.q:
                print("The total mass and/or mass ratio does not match the sum of m1 and m2. Redefining based on m1 and m2")
                self.M = self.m1 + self.m2
                self.q = self.m1/self.m2

        # Default mass parameters if three mass params are givien - include a warning
        if all([hasattr(self,"m1"), hasattr(self,"m2"), hasattr(self,"M"), not hasattr(self,"q")]) or all([hasattr(self,"m1"), hasattr(self,"m2"), not hasattr(self,"M"), hasattr(self,"q")]):
            print("You have specified three mass parameters- this is uneeded but I will continue.")
            if hasattr(self,"M") and self.m1 + self.m2 != self.M:
                    print("The total mass does not match the sum of m1 and m2. Redefining based on m1 and m2")
                    self.M = self.m1 + self.m2
                    self.q = self.m1/self.m2
            if hasattr(self,"q") and self.m1/self.m2 != self.q:
                    print("The mass ratio does not match the sum of m1 and m2. Redefining based on m1 and m2")
                    self.M = self.m1 + self.m2
                    self.q = self.m1/self.m2
        if all([hasattr(self,"M"), hasattr(self,"q"), hasattr(self,"m1"), not hasattr(self,"m2")]) or all([hasattr(self,"M"), hasattr(self,"q"), not hasattr(self,"m1"), hasattr(self,"m2")]):
            print("You have specified three mass parameters- this is uneeded but I will continue.")
            if hasattr(self,"m1") and self.m1 != self.M*(self.q/(1.+self.q)):
                    print("The primary mass m1 does not match the total mass and mass ratio. Redefining based on M and q")
                    self.m1 != self.M*(self.q/(1.+self.q))
                    self.m2 != self.M/(1.+self.q)
            if hasattr(self,"m2") and self.m2 != self.M/(1.+self.q):
                print("The secondary mass m2 does not match the total mass and mass ratio. Redefining based on M and q")
                self.m1 != self.M*(self.q/(1.+self.q))
                self.m2 != self.M/(1.+self.q)

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
                print("Assuming q=1")
                self.q = 1.
                self.m1 = self.M/2.
                self.m2 = self.m1
        elif hasattr(self,"q"):
                print("Assuming m1=1e6 M_sol")
                self.m1 = 1.e6
                self.m2 = self.m1/self.q
                self.M = self.m1+self.m2
        elif hasattr(self,"m1"):
                print("Assuming q=1")
                self.q = 1.
                self.m2 = self.m1
                self.M = self.m1+self.m2
        elif hasattr(self,"m2"):
                print("Assuming q=1")
                self.q = 1.
                self.m1 = self.m2
                self.M = self.m1+self.m2

        # Default mass parameters if no mass param is givien
        if all([not hasattr(self,"M"), not hasattr(self,"q"), not hasattr(self,"m1"), not hasattr(self,"m2")]):
                print("Assuming q=1 and m1=1e6 M_sol")
                self.m1 = 1.e6
                self.m2 = self.m1
                self.M = self.m1 + self.m2
                self.q = 1.

        # Default base source parameters
        if not hasattr(self,"z"): # If this is true then distance is also missing. Set default to redshift 3.
            print("Setting z=3")
            self.z = 3.
            self.dist = cosmo.luminosity_distance(self.z).to("Mpc").value
        if not hasattr(self,"chi1"):
            print("Setting chi1=0")
            self.chi1 = 0.
        if not hasattr(self,"chi2"):
                print("Setting chi2=0")
                self.chi2 = 0.
        if not hasattr(self,"lamda"):
                print("Setting lambda=0")
                self.lamda = 0.
        if not hasattr(self,"beta"):
                print("Setting beta=0")
                self.beta = 0.
        if not hasattr(self,"inc"):
                print("Setting inc=0")
                self.inc= 0.
        if not hasattr(self,"psi"):
                print("Setting psi=0")
                self.psi = 0.
        if not hasattr(self,"Deltat"):
                print("Setting Deltat=0")
                self.Deltat = 0.
        if not hasattr(self,"phi"):
                print("Setting phi=0")
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
                self.approximant = "IMRPhenomD"     # Included
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

        # Extra useful params
        self.true_ra = np.rad2deg(self.lamda+np.pi)   ## In deg
        self.true_dec = np.rad2deg(self.beta)         ## In deg
        self.true_distance = self.dist                ## in Mpc

        # Now check if sky_map needs creating or reading -- adaptation for 3D case needed here...
        if self.sky_map!=None and os.path.isfile(self.sky_map):
            self.LoadSkymap()
        elif self.sky_map!=None and not os.path.isfile(self.sky_map):
            self.CreateSkyMapStruct()
        elif self.sky_map==None and self.H5File!=None:
            # Default name
            self.sky_map = self.H5File.split("inference_data")[0] + 'Skymap_files' + self.H5File.split("inference_data")[-1]
            self.sky_map = self.sky_map[:-3] + '.fits'
            if os.path.isfile(self.sky_map):
                self.LoadSkymap()
            else:
                self.CreateSkyMapStruct()

    def LoadSkymap(self):
        """
        Read a skymap in -- Check that parameters match? Not sure it's needed...
        Most basic sky map here. Need to add 3d option...
        """
        self.map_struct={}
        prob_data, header = hp.read_map(self.sky_map, field=0, verbose=False,h=True)
        prob_data = prob_data / np.sum(prob_data)

        self.map_struct["prob"] = prob_data
        map_nside = hp.pixelfunc.get_nside(self.map_struct["prob"])
        npix = hp.nside2npix(nside)
        theta, phi = hp.pix2ang(nside, np.arange(npix))
        ra = np.rad2deg(phi)
        dec = np.rad2deg(0.5*np.pi - theta)

        self.map_struct["ra"] = ra
        self.map_struct["dec"] = dec

    def CreateSkyMapStruct(self,nside=None,SkyMapFileName=None):
        # default pixelation if no nside given
        if nside==None:
            nside=128 # 1024

        # if SkyMapFileName given; then preferentially save to this filename.
        # This is checked in SYU.WriteSkymapToFile just before writing to file.
        if SkyMapFileName!=None:
            self.sky_map = SkyMapFileName

        # Read in data from file
        [infer_params, _, _, _] = SYU.read_h5py_file(self.H5File)
        labels = ["lambda","beta"]
        if np.size(infer_params[labels[0]][0])>1:
            nsamples = len(infer_params["lambda"][0])
        else:
            nsamples = len(infer_params["lambda"])

        # Get pixel locations from healpy params
        npix = hp.nside2npix(nside)
        pix_thetas, pix_phis = hp.pix2ang(nside, np.arange(npix))
        pix_ras = np.rad2deg(pix_phis) # Ra and Dec for pix are stored in map_struct in degrees NOT radians
        pix_decs = np.rad2deg(0.5*np.pi - pix_thetas)
        pixels_numbers=hp.ang2pix(nside, pix_thetas, pix_phis)

        # Convert post angles to theta,phi
        post_phis = infer_params["lambda"]+np.pi
        post_thetas = np.pi/2.-infer_params["beta"]
        post_pix = hp.ang2pix(nside,post_thetas,post_phis)

        # Probabilities of pixels
        probs,bin_edges = np.histogram(post_pix, bins=np.arange(npix+1), density=True)

        # Make the dictionary
        map_struct = {"prob":probs,
                      "ra":pix_ras,
                      "dec":pix_decs}

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

        # create class attribute
        self.map_struct=map_struct

        # Save to file
        SYU.WriteSkymapToFile(self.map_struct,self.sky_map,None)

    def GenerateEMFlux(self,detector,fstart22=1e-4,**EM_kwargs):
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
        # print(type(phi_c), type(self.Deltat), type(self.m1*pyconstants.MSUN_SI), type(self.m2*pyconstants.MSUN_SI), type(fstart22), type(self.dist*1e6*pyconstants.PC_SI), type(self.inc), type(chi1_x), type(chi1_y), type(self.chi1), type(chi2_x), type(chi2_y), type(self.chi2), type(modearray), type(seobflags))
        # print((phi_c), (self.Deltat), (self.m1*pyconstants.MSUN_SI), (self.m2*pyconstants.MSUN_SI), (fstart22))
        # print((self.dist*1e6*pyconstants.PC_SI), (self.inc))
        # print(chi1_x, chi1_y, self.chi1)
        # print(chi2_x, chi2_y, self.chi2)
        # print((modearray), (seobflags))
        # I assume it takes parameters in the source frame...?
        # if self.Lframe: # second input var (delta_t) was set to 1.5... can we increase this? Should 1st parameter be self.psi?
        #     print("Lframe; rescaling to source fram using z =",self.z)
        #     hplus, hcross, hIlm, hJlm, dyn_Low, dyn_Hi, dyn_all, t_vec_modes, hP22_amp, hP22_phase, hP21_amp, hP21_phase, hP33_amp, hP33_phase, hP44_amp, hP44_phase, hP55_amp, hP55_phase, alphaJ2P, betaJ2P, gammaJ2P, AttachPars = lalsim.SimIMRSpinPrecEOBWaveformAll(self.phi, 1.5, self.m1*pyconstants.MSUN_SI/(1.+self.z), self.m2*pyconstants.MSUN_SI/(1.+self.z), fstart22/(1.+self.z), self.dist*1e6*pyconstants.PC_SI, self.inc, chi1_x, chi1_y, self.chi1, chi2_x, chi2_y , self.chi2, modearray, seobflags)
        # else:
        #     hplus, hcross, hIlm, hJlm, dyn_Low, dyn_Hi, dyn_all, t_vec_modes, hP22_amp, hP22_phase, hP21_amp, hP21_phase, hP33_amp, hP33_phase, hP44_amp, hP44_phase, hP55_amp, hP55_phase, alphaJ2P, betaJ2P, gammaJ2P, AttachPars = lalsim.SimIMRSpinPrecEOBWaveformAll(self.phi, 1.5, self.m1*pyconstants.MSUN_SI, self.m2*pyconstants.MSUN_SI, fstart22, self.dist*1e6*pyconstants.PC_SI, self.inc, chi1_x, chi1_y, self.chi1, chi2_x, chi2_y , self.chi2, modearray, seobflags)
        hplus, hcross, hIlm, hJlm, dyn_Low, dyn_Hi, dyn_all, t_vec_modes, hP22_amp, hP22_phase, hP21_amp, hP21_phase, hP33_amp, hP33_phase, hP44_amp, hP44_phase, hP55_amp, hP55_phase, alphaJ2P, betaJ2P, gammaJ2P, AttachPars = lalsim.SimIMRSpinPrecEOBWaveformAll(self.phi, 1.5, self.m1*pyconstants.MSUN_SI, self.m2*pyconstants.MSUN_SI, fstart22, self.dist*1e6*pyconstants.PC_SI, self.inc, chi1_x, chi1_y, self.chi1, chi2_x, chi2_y , self.chi2, modearray, seobflags)
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
        param_dict,waveform_params,extra_params = SYU.ClassesToParams(self,detector,**EM_kwargs) ## Do we need to change the approximant here from the saved one in source to EOB one?
        waveform_params["DeltatL_cut"] = None # Need to set this to 0 to get the remaining signal for photon phases
        # print("Params from 'ClassesToParams':",param_dict,waveform_params,extra_params)
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
            print(r"T$_{SNR=10}$",(1.+self.z)*t_start_flux/(24.*60.*60.),r"t$_{cut}$", self.DeltatL_cut/(24.*60.*60.))
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

        # # Plot some tests of frequency versus time
        # # Calculate their frequency - time relation for spinless systems -- ALL IN *SOURCE* frame then convert to *DETECTOR* frame
        # M_reduced = (self.m1*self.m2/(self.M*(1.+self.z)))*pyconstants.MSUN_SI*(pyconstants.G_SI/(pyconstants.C_SI**3))
        # eta = M_reduced/(M_tot*pyconstants.MSUN_SI*(pyconstants.G_SI/(pyconstants.C_SI**3)))
        # M_chirp = (M_reduced**(3/5))*((M_tot*pyconstants.MSUN_SI*(pyconstants.G_SI/(pyconstants.C_SI**3)))**(2/5))
        # t_of_f_paper = [-((5/256)*M_chirp*(np.pi*M_chirp*f/2)**(-8/3))*(1.+1.33333*(743/336+11*eta/4)*(np.pi*M_chirp*f/2)**(2/3)-(np.pi*M_chirp*f/2)*32*np.pi/5 + 2*(3058673/1016064 + 5429*eta/1008 + 617*eta*eta/144)*(np.pi*M_chirp*f/2)**(4/3)) for f in GW_freqs]
        # import matplotlib.pyplot as plt
        # import matplotlib.pylab as pylab
        # params = {'legend.fontsize': 8, # 'x-large',
        #          'axes.labelsize': 8, # 'x-large',
        #          'xtick.labelsize': 4, # 'x-large',
        #          'ytick.labelsize': 4, # 'x-large'}
        #          'lines.markersize': 2}
        # pylab.rcParams.update(params)
        # times = [time*(1.+self.z)/(24.*60.*60.) for time in xray_time]
        # freqs_plot = [f/(1.+self.z) for f in GW_freqs]
        # times_paper = [(1.+self.z)*time*M_tot*pyconstants.MSUN_SI*(pyconstants.G_SI/(pyconstants.C_SI**3))/(24.*60.*60.) for time in t_of_f_paper]
        # freqs_plot_paper = [f/(1.+self.z) for f in GW_freqs]
        # t_end_flux_paper = [time for (time,GW_freq) in zip(times_paper,GW_freqs) if GW_freq<ISCO_freq]
        # if isinstance(t_end_flux_paper,list) and len(t_end_flux_paper)>0:
        #     t_end_flux_paper = t_end_flux_paper[-1]
        # elif isinstance(t_end_flux_paper,list):
        #     t_end_flux_paper=t_end_flux_paper[-1]
        # t_start_flux_paper = [time_p for (time_p,time) in zip(times_paper,xray_time) if time<t_start_flux][-1]
        # print("Detector frame test t$_{start}$:",t_start_flux_paper,"Detector frame test t$_{end}$:",t_end_flux_paper)
        # plt.plot(times, freqs_plot)
        # plt.plot(times_paper, freqs_plot_paper)
        # plt.xlabel('Time from Merger [d]')
        # plt.ylabel('GW Frequency [Hz]')
        # ax = plt.gca()
        # ax.set_yscale('log')
        # plt.grid()
        # plt.show()

        # Set the class variables using the start and end times - convert everything *EXCEPT* xray_flux to *DETECTOR* frame
        self.xray_flux = [flux for (flux,time) in zip(xray_flux,xray_time) if time>t_start_flux and time<t_end_flux] # if time>t_start_flux and time<t_end_flux else 0. for (flux,time) in zip(xray_flux,xray_time)]
        self.xray_time = [time*(1.+self.z) for time in xray_time if time>t_start_flux and time<t_end_flux]
        self.GW_phi = [phi for (phi,time) in zip(GW_phi,xray_time) if time>t_start_flux and time<t_end_flux] # GW_phi
        self.GW_Omega = [Om/(1.+self.z) for (Om,time) in zip(GW_Omega,xray_time) if time>t_start_flux and time<t_end_flux] # [Om/(1.+self.z) for Om in GW_Omega]
        self.r = [rad for (rad,time) in zip(r,xray_time) if time>t_start_flux and time<t_end_flux] # r
        self.GW_freqs = [f/(1.+self.z) for (f,time) in zip(GW_freqs,xray_time) if time>t_start_flux and time<t_end_flux] # [f/(1.+self.z) for f in GW_freqs]

        # Plot cumulative SNR
        PLOT_CUMUL_SNR = False
        if PLOT_CUMUL_SNR:
            import matplotlib.pyplot as plt
            import matplotlib.pylab as pylab
            params = {'legend.fontsize': 8, # 'x-large',
                     'axes.labelsize': 8, # 'x-large',
                     'xtick.labelsize': 4, # 'x-large',
                     'ytick.labelsize': 4, # 'x-large'}
                     'lines.markersize': 2}
            pylab.rcParams.update(params)
            times = [t/(24.*60.*60.) for t in tf]
            plt.plot(times, cumul_SNR)
            ax = plt.gca()
            ax.set_yscale('log')
            plt.xlabel('Time [d]')
            plt.ylabel('Cumul SNR')
            plt.grid()
            plt.show()

        # Plot radial distance
        DO_RADIUS_PLOT = False
        if DO_RADIUS_PLOT:
            print("Plotted system parameters:", self.M, self.q, self.z)
            import matplotlib.pyplot as plt
            import matplotlib.pylab as pylab
            params = {'legend.fontsize': 8, # 'x-large',
                     'axes.labelsize': 8, # 'x-large',
                     'xtick.labelsize': 4, # 'x-large',
                     'ytick.labelsize': 4, # 'x-large'}
                     'lines.markersize': 2}
            pylab.rcParams.update(params)
            times = [time/(24.*60.*60.) for time in self.xray_time]
            plt.plot(times, self.r)
            ax = plt.gca()
            ax.set_yscale('log')
            plt.xlabel('Time [d]')
            plt.ylabel('EOB Rad [R$_S$]')
            plt.grid()
            plt.show()

        # Plot GW frequency as function of time
        DO_GW_FREQ_PLOT = False
        if DO_GW_FREQ_PLOT:
            import matplotlib.pyplot as plt
            import matplotlib.pylab as pylab
            params = {'legend.fontsize': 8, # 'x-large',
                     'axes.labelsize': 8, # 'x-large',
                     'xtick.labelsize': 4, # 'x-large',
                     'ytick.labelsize': 4, # 'x-large'}
                     'lines.markersize': 2}
            pylab.rcParams.update(params)
            times = [time/(24.*60.*60.) for time in self.xray_time]
            plt.plot(times, self.GW_freqs)
            plt.xlabel('Time from Merger [d]')
            plt.ylabel('GW Frequency [Hz]')
            ax = plt.gca()
            ax.set_yscale('log')
            plt.grid()
            plt.show()

        # Plot x-ray flux and velocity
        DO_XRAY_PLOT = True
        if DO_XRAY_PLOT:
            print("Plotted system parameters:", self.M, self.q, self.z)
            import matplotlib.pyplot as plt
            import matplotlib.pylab as pylab
            params = {'legend.fontsize': 8, # 'x-large',
                     'axes.labelsize': 8, # 'x-large',
                     'xtick.labelsize': 4, # 'x-large',
                     'ytick.labelsize': 4, # 'x-large'}
                     'lines.markersize': 2}
            pylab.rcParams.update(params)
            times = [time/(24.*60.*60.) for time in self.xray_time]
            plt.plot(times, self.xray_flux)
            plt.xlabel('Time from Merger [d]')
            plt.ylabel('X-ray Flux [erg s$^{-1}$ cm$^-2$]')
            ax = plt.gca()
            ax.set_yscale('log')
            plt.grid()
            plt.show()

        # Plot EOB radial velocity for fun
        DO_RAD_VEL_PLOT = False
        if DO_RAD_VEL_PLOT:
            print("Plotted system parameters:", self.M, self.q, self.z)
            import matplotlib.pyplot as plt
            import matplotlib.pylab as pylab
            params = {'legend.fontsize': 8, # 'x-large',
                     'axes.labelsize': 8, # 'x-large',
                     'xtick.labelsize': 4, # 'x-large',
                     'ytick.labelsize': 4, # 'x-large'}
                     'lines.markersize': 2}
            pylab.rcParams.update(params)
            plt.plot(xray_time, r_dot)
            plt.xlabel('Time [s]')
            plt.ylabel('EOB Rad Vel [c]')
            plt.grid()
            plt.show()

    def GenerateCTR(self,ARF_file_loc_name,gamma=1.7):
        self.xray_gamma = gamma
        if not hasattr(self,"xray_flux"):
            raise ValueError("No x-ray flux generated. Need to call GenerateEMFlux with a detector object before running this function.")

        self.xray_phi_0 = [(6.242e8)*self.xray_flux[ii]*(2.-self.xray_gamma)/(10.**(2.-self.xray_gamma)-0.2**(2.-self.xray_gamma)) for ii in range(len(self.xray_flux))]
        from astropy.io import fits
        hdul = fits.open(ARF_file_loc_name)
        hdul.info()
        ARF = hdul[1].data[:]
        N = len(ARF)
        E_bins = [ARF[ii][0] for ii in range(N)] # Units are keV
        E_bins.append(ARF[-1][1]) # So that this length is N+1 where N is the length of ARF_func
        ARF_func = [ARF[ii][2] for ii in range(N)]
        dE_bins = np.diff(E_bins)

        # Now compute the CTR - photons s^-1
        integrand = [dE_bins[ii]*ARF_func[ii]*((E_bins[ii]+E_bins[ii+1])/2.)**(-1.*self.xray_gamma) for ii in range(N)]
        integral = sum(integrand)/(1.+self.z) # We need another 1/(1+z) here otherwise the photon rate is too high. Maybe this is a redshifting of energy bins?
        # CTR = [integral*self.xray_phi_0[ii] for ii in range(len(self.xray_phi_0))]
        self.CTR = [integral*self.xray_phi_0[ii] for ii in range(len(self.xray_phi_0))] # CTR
