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
#
######################################################################################

import csv
import astropy.constants as apc
import numpy as np
import sys
import os.path
from getQ import getQ
import pickle
from math import floor, pi
from scipy.io import readsav

from globals import HT04_globals
import SynspecSettings as ss
from create_filename import create_filename
from makegrid import grid
import matplotlib.pyplot as plt

c2 = (apc.h.cgs.value*apc.c.cgs.value)/apc.k_B.cgs.value
 
 
 
 
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

class LineList:
    """A basic line list constructor"""
    def __init__(self, lines):  # need to add options
        self.lines = lines
        spec = np.array([x.spec for x in lines])
        iso  = np.array([x.iso for x in lines])
        
        self.freq = np.array([x.freq for x in lines])
        self.wave = np.array([x.wave for x in lines])
        self.waveum = np.array([x.waveum for x in lines])
        self.eA =  np.array([x.eA for x in lines])
        self.epp = np.array([x.epp for x in lines])
        self.strength = np.array([x.strength for x in lines]) 
    
        if np.var(spec) == 0.0:
            self.spec = np.median(spec)
        else:
            sys.exit('Trying to construct linelist from multiple species') # Kill it if wrong

        if np.var(iso) == 0.0:
            self.iso = np.median(iso)
        else:
                sys.exit('Trying to construct linelist from multiple isotopologues') # Kill it if wrong

    def addline(self):
        print 'To be implemented later if necessary'
    
    def removeline(self):
        print 'To be implemented later if necessary'

    def ll_reverse(self):
        'Linelists normally presented as sorted in wavenumber. We want sorted by um.'
        self.freq = self.freq[::-1]
        self.wave = self.wave[::-1]
        self.waveum = self.waveum[::-1]
        self.eA = self.eA[::-1]
        self.epp = self.epp[::-1]
        self.strength = self.strength[::-1] 
    
class SpecificLineList(LineList):
    """Container for methods fetching hitran04 line lists"""
    
    specs_calced = 0
    
    def __init__(self, molno, isono, ll_name, regen=False):
        
        print molno, isono
        self.ll_name = ll_name
        props = HT04_globals(molno, isono)
        
        if ll_name == 'HITRAN04':
            filename = create_filename(molno, isono, ll_name, 'linelist')
            #filename = 'filereadtest/01_01_hit04.par'
            
            if os.path.isfile(filename) == False:
                sys.exit('File not found:'+filename)
            
            if os.path.isfile(filename+'.pickle') == False or regen==True:  # Here I am checking whether we already have an object for this linelist.
                #read in the file - first find number of rows
                with open(filename) as csvfile:
                    reader = csv.reader(csvfile)
                    i=0
                    for row in reader:
                        #print i
                        i+=1
        
                        lines = np.arange(i,dtype=object)
        
                #now loop through and fill the array with lines
                with open(filename) as csvfile:
                    reader = csv.reader(csvfile)
                    i=0
                    for row in reader:
                        #print i
                        #print('{:2d}{:1d}{:12.6f}{:10.3e}{:10.3e}{:5.4f}{:5.4f}{:10.4f}{:4.2f}{:8.6f}{:15c}{:15c}{:15c}{:15c}{:6d}{:12d}{:1c}{:7.1f}{:7.1f}'.format(row[0]))
                        data = row[0]
                        thisLine = Line( np.uint(data[0:2]), np.uint(data[2:3]), np.float64(data[3:15]), np.float64(data[25:35]), strength=np.float64(data[15:25]), epp=np.float64(data[45:55]))
                        lines[i] = thisLine
                        i+=1

                outfile = open(filename+'.pickle','w')
                pickle.dump(lines, outfile)
                outfile.close
            else:
                infile = open(filename+'.pickle','r') #Can only get to this block if the .pickle file exists. Restore it, its faster.
                lines = pickle.load(infile)

            self.wave = np.array([x.wave for x in lines]) # currently this is actually wno
            self.waveum = np.array([x.waveum for x in lines])
            self.epp = np.array([x.epp for x in lines])
            self.eA = np.array([x.eA for x in lines])
            self.freq = np.array([x.freq for x in lines])
            E_temp = -1.0*self.epp * c2 / 296.0
            w_temp = -1.0*self.wave * c2 / 296.0
            self.strength = np.array([x.strength for x in lines]) * (props.Q296/props.abund) / (np.exp(E_temp) * (1.0-np.exp(w_temp)) )
            print np.array([x.strength for x in lines]) 
            spec = np.array([x.spec for x in lines])
            iso = np.array([x.iso for x in lines])
        
        
            if np.var(spec) == 0.0:
                self.spec = int(np.median(spec))
            else:
                sys.exit('Trying to construct linelist from multiple species') # Kill it if wrong

            if np.var(iso) == 0.0:
                self.iso = int(np.median(iso))
            else:
                sys.exit('Trying to construct linelist from multiple isotopologues') # Kill it if wrong
      
           
        
        else:
            print 'other Linelists Not Implemented'             
    
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
                
     
                E_temp = -1.0 * self.epp * c2 / Temp
                print E_temp
                w_temp = -1.0 * self.wave * c2 / Temp
                print w_temp
                self.strength = self.strength * (props.abund/ Q) * (np.exp(E_temp) * (1.0-np.exp(w_temp))) * apc.c.cgs.value
                #I have no idea why Jan multiplies by C here, but he does, so lets copy it.
                
                strengths_jan = readsav('/home/dstock/sf/idl/code/ff.xdr')
                
                print "My Calcs:", self.strength
                print 'My epp', self.epp
                print "Jan's Calcs:", strengths_jan.ff

                print self.strength[0:12]/strengths_jan.ff

                print strengths_jan.ff[0]/self.strength[0]

                #sys.exit('get outta town')
                

                          
class Spectrum:
#    """Container for actual spectra"""
#    #Convert Jan's default spectra into something python can read and make it default
    def __init__(self, wave_start=ss.wavestart, wave_end=ss.waveend, v_turb=ss.vturb, 
                 resolution=ss.resolution, oversample=ss.oversample, xtitle="", 
                 ytitle="", edges=False):
        """ This function will set up a default wavelength array (in wavenumbers) as well as a container for actual spectra """
        self.grid = grid(wave_start, wave_end, v_turb, resolution, oversample, xtitle, ytitle, edges)
        # This makes sure that when we make a spectrum object with no args, it defaults to containing the basic wave grid.

    
    def get_tau(self, molno, isono, ll_name, Temp, regen=False):
        """General form of Jan's solution to this problem:
            make a 2d array of strengths, delta_nu's and line centers
            make a 2d tau array with form (npoints * nlines) where the npoints represent the array of frequencies
            subtract the line center array (npoints * nlines) where every row is just frequency of the line center from the tau array 
            tau is now equal to nu - nu0
            divide the whole thing by your 2d delta nu array 
            tau now = (nu-nu0)/dnu
            square the whole thing
            mask out the crap
            do e(-tau)
            do a bit of renomalization
            all done"""
    
        thislinelist = SpecificLineList(molno, isono, ll_name, regen)
        thislinelist.calc_specifics(Temp)
        #thislinelist.ll_reverse()
        
        print thislinelist.strength.size, type(thislinelist.strength.size)
        print self.grid.n_points
        print  type(thislinelist.strength)

        mastertau = self.grid.flux
        mastertau[:] = 0.0
        
        delta_nud = thislinelist.wave * ss.vturb * 1e5 
    
        # assume that every incoming linelist is sorted.
        #split spectrum into chunks of say 100 lines
        #lets just make it work - no overlaps yet.
        
        
        nlines=100
        cycle = 0
        iters = floor( thislinelist.strength.size/nlines )
        
        print thislinelist.strength.size
        print floor(thislinelist.strength.size/nlines)
        #sys.exit()
        a_size=nlines

        thisgrid = self.grid.freq[::-1]
        
        for cycle in range(int(iters)+1):
            print 'cycle=', cycle

            print 'line number range=', (cycle)*nlines, (cycle+1)*nlines, thislinelist.strength.size

            if cycle == int(iters):
                print 'in here'
                thislooplines_s = thislinelist.strength[(cycle)*nlines:]
                thislooplines_w = thislinelist.freq[(cycle)*nlines:]
                thisdeltanu = delta_nud[(cycle)*nlines:]
                
                a_size = thislooplines_s.size
            else:
                thislooplines_s = thislinelist.strength[(cycle)*nlines:(cycle+1)*nlines]
                thislooplines_w = thislinelist.freq[(cycle)*nlines:(cycle+1)*nlines]
                thisdeltanu = delta_nud[(cycle)*nlines:(cycle+1)*nlines]
            
            
            print 'wavelength range=', thislooplines_w[0], thislooplines_w[-1]
            
            if cycle==0:
                linerange = [thislooplines_w[0],thislooplines_w[-1] ]
            else:
                linerange = [thislinelist.freq[(cycle*nlines)-1], thislooplines_w[-1] ]
            
            print 'linerange=', linerange
            #wavecov = linerange[1] - linerange[0]
            
            print self.grid.freq[0], self.grid.freq[-1]
            
            if cycle==0:
                thiswavegrid = thisgrid[(thisgrid < linerange[1])]
            elif cycle==int(iters):
                thiswavegrid = thisgrid[(thisgrid > linerange[0])]
            else:
                thiswavegrid = thisgrid[(thisgrid > linerange[0]) & (thisgrid < linerange[1])]
           
                      
            print 'Nlines, Npoints', thislooplines_s.size, thiswavegrid.size
        
            
            ff2d = np.vstack([thislooplines_s]*thiswavegrid.size)
            lines2d = np.vstack([thislooplines_w]*thiswavegrid.size)
            dnu_2d = np.vstack([thisdeltanu]*thiswavegrid.size)
            wavegrid2d = np.transpose(np.vstack([thiswavegrid]*a_size))
            
            temp = wavegrid2d - lines2d
            temp = temp/dnu_2d
            temp = temp**2.0
            temp = -temp
            temp = np.exp( temp )
            
            #TODO: Mask crap here
            
            temp = temp * ff2d
            temp = temp / dnu_2d
            temp = temp / np.sqrt(pi)
            
            thistau = np.sum(temp, axis=1)
            
            if cycle==0:
                mastertau[(thisgrid < linerange[1])] = thistau
            elif cycle==int(iters):
                mastertau[(thisgrid > linerange[0]) ] = thistau
            else:    
                mastertau[(thisgrid > linerange[0]) & (thisgrid < linerange[1])] = thistau
            
    
        
        jan = readsav('/home/dstock/sf/idl/code/templates/CO2/626/HITRAN04/v3.000/gauss/tau/CO2_626_HITRAN04_v3.000_NTH_gauss_T100.000.xdr')
        plt.plot(self.grid.wave, mastertau[::-1]*1.1, 'b', jan.wave, jan.tau, 'r')
        plt.show()
    
    
    
    



data = Spectrum()
data.get_tau(2,1,'HITRAN04',100, regen=False)
    
    
    
    
    
    
    
    
    
    
    