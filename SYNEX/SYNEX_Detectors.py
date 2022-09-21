import astropy.units as u
import astropy.stats as astat
from astropy.time import Time
import astropy, astroplan
import numpy as np
import os
import SYNEX.SYNEX_Utils as SYU
from SYNEX.SYNEX_Utils import SYNEX_PATH
from datetime import date
import pathlib

# mpi stuff
try:
    from mpi4py import MPI
except ModuleNotFoundError:
    MPI = None


##############
#
# Possible structural change is to bring all lisabeta inference master functions
# into the LISA class since they are exclusively for LISA inference. Would make this
# easier to read and emphasise to future developers that lisabeta is for LISA only.
#
##############


class LISA:
    """
    Class to make a Space Based Instrument using the Interferometer base class

    Parameters
    ----------



    """

    def __init__(self, **kwargs):
        # Set verbosity
        self.verbose=kwargs.pop("verbose") if "verbose" in kwargs else True

        # Now unpack other kwargs -- Need to flatten this block of code cause its ugly af
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
            elif key == "gps_science_start":
                self.gps_science_start = value
            elif key == "mission_duration":
                self.mission_duration = value

        # Defaults if not anything not specified in kwargs
        if not hasattr(self,"tmin"):
            self.tmin = None
        if not hasattr(self,"tmax"):
            self.tmax = None
        if not hasattr(self,"TDI"):
            self.TDI = "TDIAET"
        if not hasattr(self,"order_fresnel_stencil"):
            self.order_fresnel_stencil = 0
        if not hasattr(self,"LISAconst"):
            self.LISAconst = "Proposal"
        if not hasattr(self,"responseapprox"):
            self.responseapprox = "full"
        if not hasattr(self,"frozenLISA"):
            self.frozenLISA = False
        if not hasattr(self,"TDIrescaled"):
            self.TDIrescaled = True
        if not hasattr(self,"gps_science_start"):
            self.gps_science_start = Time('2033-01-01T00:00:00.00', format='isot', scale='utc').gps
        if not hasattr(self,"mission_duration"):
            self.mission_duration = 3.

        # Now replace LISA noise params if mentioned in kwargs (replacing any LISAnoise existing dict)
        if not hasattr(self,"LISAnoise"): self.LISAnoise = {}
        self.LISAnoise["InstrumentalNoise"]=kwargs["InstrumentalNoise"] if "InstrumentalNoise" in kwargs else "SciRDv1"
        self.LISAnoise["InstrumentalNoise"]=kwargs["InstrumentalNoise"] if "WDbackground" in kwargs else True
        self.LISAnoise["InstrumentalNoise"]=kwargs["InstrumentalNoise"] if "WDduration" in kwargs else self.mission_duration
        self.LISAnoise["InstrumentalNoise"]=kwargs["InstrumentalNoise"] if "lowf_add_pm_noise_f0" in kwargs else 0.
        self.LISAnoise["InstrumentalNoise"]=kwargs["InstrumentalNoise"] if "lowf_add_pm_noise_alpha" in kwargs else 2.
