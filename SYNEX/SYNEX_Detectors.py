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
