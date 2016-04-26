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
# This file defines the chunk object, which is intended to serve as a container
# for chunks of linelists. One day merge this back into a linelist object I guess
# but for now trying just to get it to work. 20/4/2016
#
#
#
# VERSION HISTORY:
# Created: 20/04/2016 (DJS)
#
# This whole thing would work better as a linelist object if we have an appropriate constructor.
######################################################################################

from globals import HT04_globals
from getQ import getQ
import astropy.constants as apc
import numpy as np


class chunk():
    def __init__(self, wave, freq, strength, lineinds, gridinds, nchunks, ll_name, spec, iso, outwithgrid, epp):
        self.nchunks = nchunks
        self.waveum = wave
        self.freq = freq
        self.strength = strength
        self.lineinds = lineinds
        self.gridinds = gridinds
        self.specs_calced = 0
        self.ll_name = ll_name
        self.spec = spec
        self.iso = iso
        self.outwithgrid = outwithgrid
        self.epp = epp
        
    def calc_specifics(self, Temp):
        """A separate method to calculate the specific line list properties based on an input T."""
        if self.specs_calced == 0:
            #make sure we don't inadvertently try and do this twice
            if self.ll_name == 'HITRAN04':
                self.Temp = Temp
                self.specs_calced = 1
                #lets make sure the relevant temperature is now carried around with the linelist.                
                
                props = HT04_globals(self.spec, self.iso)
                
                if Temp == 296.0 and self.ll_name == 'HITRAN04':
                    Q=props.Q296
                else:
                    Q=getQ(self.spec, self.iso, self.ll_name, Temp)    
                
                c2 = (apc.h.cgs.value*apc.c.cgs.value)/apc.k_B.cgs.value
     
                E_temp = -1.0 * self.epp * c2 / Temp
                #print E_temp
                w_temp = -1.0 * (10000/self.waveum) * c2 / Temp
                #print w_temp
                self.strength = self.strength * (props.abund/ Q) * (np.exp(E_temp) * (1.0-np.exp(w_temp))) * apc.c.cgs.value
        
        else:
            print 'Linelist already transformed to appropriate strengths.'
