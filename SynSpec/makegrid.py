######################################################################################
#
# Python Reimplementation of Spectra Factory (by Jan Cami)
# 
# This file is a 
# ( ) Translation of an original .pro file (called ...)
# (X) Reimplementation of an original file (called make_grid.pro)
# ( ) New file
#
# COMMENTS:
# This routine makes arbitrary wavelength grids.
#
#
#
# VERSION HISTORY:
# Created: 12/02/2016 (DJS)
# Made the wave array explicitly have the type float64: 19/02/2016
# Corrected such that wave grids with arbitrary resolutions actually work
#
######################################################################################

import numpy as np
import numpy.ma as ma

import SynspecSettings as ss
import astropy.constants as apc



class grid():
    def __init__(self, wave_start=ss.wavestart, wave_end=ss.waveend, v_turb=ss.vturb, 
                 resolution=ss.resolution, oversample=ss.oversample, xtitle="", 
                 ytitle="", edges=False): 
        
        """ Need to be careful with this, basically v_turb is always the dominant factor, however
        the IDL has no traps for putting v_turb = 0, which would crash the whole thing. In time
        we will need to update this such that there is always an underlying element of thermal
        broadening such that the line profiles cannot go to zero width. For now, v_turb is going
        to be the key parameter.  I have had to remove the option to return just a grid. Always an
        object now."""
        
        if resolution == ss.resolution:
            resolution = apc.c.cgs.value/(v_turb*1.0e5)
        # Behaviour as follows:
        # If resolution is set to be other than default value, it is used (i.e. making output grids)
        # If resolution is left at default value, we recalculate it based on the default v_turb (i.e. making default grid)
                    
        
        R = np.float(resolution) * np.float(oversample)
        f = -(1.0 + 2.0 * R) / (1.0 - 2.0 * R)
        n_points = (np.floor( (np.log(wave_end/wave_start)) / (np.log( f )))) + 1
        
        #print 'R: ', R, ' f: ', f, ' top:', (np.log(wave_end/wave_start)), ' bottom: ', np.log(f)
        
        #print 'makegrid n: ',n_points, ' resolution: ', resolution
        
        wave = np.arange(n_points, dtype=np.float64)
        wave.fill(wave_start)
        factor = f**np.arange(n_points, dtype=np.float64)
        wave = wave*factor
        
        # Checked against IDL version

        flux = np.arange(n_points)
        freq = (apc.c.cgs.value * 1.0e4) / wave
        
        dl = wave/resolution

        xstart = wave - dl/2.
        xend = wave + dl/2.

        if edges == True:
            """ Copy code to take care of edges here"""
            idx = ma.where(xstart < wave_start)
            if idx[0] != []:
                print "in here"
                xstart[idx] = wave_start
                dl[idx] = xend[idx]-xstart[idx]
            idx = ma.where(xstart > wave_end)
            if idx[0] != []:
                print "in here"
                xend[idx] = wave_end
                dl[idx] = xend[idx]-xstart[idx]
    
        dnu = (apc.c.cgs.value * 1.0e4 * dl) / wave**2

        self.wave = wave # this wave in um
        self.flux = flux
        self.xtitle = xtitle
        self.ytitle = ytitle
        self.n_points = n_points
        self.resolution = resolution
        self.oversample = oversample
        self.freq = freq
        self.dl = dl
        self.dnu = dnu
        self.xstart = xstart
        self.xend = xend


#test = grid(edges=True)
#print test.wave