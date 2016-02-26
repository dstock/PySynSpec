######################################################################################
#
# Python Reimplementation of Spectra Factory (by Jan Cami)
# 
# This file is a 
# ( ) Translation of an original .pro file (called ...)
# ( ) Reimplementation of an original file (called ...)
# (X) New file
#
# COMMENTS:
# This file defines the various classes needed to implement spectrafactory.
#
#
#
# VERSION HISTORY:
# Created: 04/02/2016 (DJS)
# V0.001: Can create Tau Profiles for HITRAN04 linelists 18/02/2016
# V0.01: Can make output templates whoch match IDL versions for arbitrary
#        resolutions, v_turbs, and linelists from HITRAN04 database.
#        The line chunking routine in the tau calculation is very basic
#        and needs updating and steamlining. In particular the leading 
#        chunk of each wavelength grid is an arbitrarily huge array of 
#        zeros which we're still handling.However it seems to be correct.
#        Extensively tested using CO 626 (molno=2, isono=1).
#        
#        Split various classes into separate files (26/02/16, DJS)
#
######################################################################################

import astropy.constants as apc


class Line:
    """Container for all the info about one absorption line."""
    # Assumes wave given in microns
    def __init__(self, species, iso, wave, eA, epp=0, strength=0): #will add many more, also type checking
        self.spec = species
        self.iso = iso
        self.wave = wave # currently in wno, add um
        self.waveum = 10000.0/wave
        self.eA = eA
        self.freq = apc.c.cgs.value * 1e4 / self.waveum
        self.epp = epp
        self.strength = strength