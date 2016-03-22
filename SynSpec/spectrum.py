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
# This file defines Spectrum object, which contains wavelength and flux grids.
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
import numpy as np
import sys
import os.path
import pickle
from math import floor, pi
from scipy.io import readsav
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy import signal
from scipy import interpolate
import time
import numdifftools as nd


from linelist import SpecificLineList
import SynspecSettings as ss
from create_filename import create_filename
from makegrid import grid

from b_nu import B_nu
from _dbus_bindings import String

c2 = (apc.h.cgs.value*apc.c.cgs.value)/apc.k_B.cgs.value
                           
class Spectrum:
#    """Container for actual spectra"""
#    #Convert Jan's default spectra into something python can read and make it default
    def __init__(self, wave_start=ss.wavestart, wave_end=ss.waveend, v_turb=ss.vturb, 
                 resolution=ss.resolution, oversample=ss.oversample, xtitle="", 
                 ytitle="", edges=False):
        """ This function will set up a default wavelength array (in wavenumbers) as well as a container for actual spectra """
        self.grid = grid(wave_start, wave_end, v_turb, resolution, oversample, xtitle, ytitle, edges)
        # This makes sure that when we make a spectrum object with no args, it defaults to containing the basic wave grid.
        self.spectrumtype='init'
        
    
    def get_tau(self, molno, isono, ll_name, Temp, regen=False, vturb = ss.vturb):
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
        
        self.Temp = Temp
        
        if self.spectrumtype != 'init':
            sys.exit('This spectrum object is not blank.')
        
        #lets make sure that we haven't done this before..
        filename = create_filename(molno, isono, ll_name, 'tau', vturb, Temp=self.Temp)
        
        if os.path.isfile(filename) == False or regen == True:
            
            thislinelist = SpecificLineList(molno, isono, ll_name, regen)
            thislinelist.calc_specifics(Temp)
            #thislinelist.ll_reverse()
            
            print thislinelist.strength.size, type(thislinelist.strength.size)
            print self.grid.n_points
            print  type(thislinelist.strength)
    
            mastertau = self.grid.flux
            mastertau[:] = 0.0
            
            delta_nud = thislinelist.wave * vturb * 1e5 
        
            # assume that every incoming linelist is sorted.
            #split spectrum into chunks of say 100 lines
            #lets just make it work - no overlaps yet.
            
            
            nlines=100
            cycle = 0
            iters = floor( thislinelist.strength.size/nlines )
            
            #print thislinelist.strength.size
            #print floor(thislinelist.strength.size/nlines)
            #sys.exit()
            a_size=nlines
    
            thisgrid = self.grid.freq[::-1]
            
            for cycle in range(int(iters)+1):
                print 'cycle=', cycle
    
                print 'line number range=', (cycle)*nlines, (cycle+1)*nlines, thislinelist.strength.size
    
                if cycle == int(iters):
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
                
                #print self.grid.freq[0], self.grid.freq[-1]
                
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
                
                mask = ma.where(temp < 500)
                                
                temp[mask] = temp[mask] * ff2d[mask]
                temp[mask] = temp[mask] / dnu_2d[mask]
                temp[mask] = temp[mask] / np.sqrt(pi)
                
                thistau = np.sum(temp, axis=1)
                
                if cycle==0:
                    mastertau[(thisgrid < linerange[1])] = thistau
                elif cycle==int(iters):
                    mastertau[(thisgrid > linerange[0]) ] = thistau
                else:    
                    mastertau[(thisgrid > linerange[0]) & (thisgrid < linerange[1])] = thistau
                
        
            self.tau = mastertau[::-1] #The -1 reverses the array to align it with the um wavelength grid
            outfile = open(filename,'w')
            pickle.dump(mastertau, outfile)
            outfile.close()
        else:
            infile = open(filename,'r') #Can only get to this block if the .pickle file exists. Restore it, its faster.
            mastertau = pickle.load(infile)
            infile.close()
            self.tau = mastertau[::-1] #The -1 reverses the array to align it with the um wavelength grid
            print 'restored tau pickle'           
                
        #jan = readsav('/home/dstock/sf/idl/code/templates/CO2/626/HITRAN04/v3.000/gauss/tau/CO2_626_HITRAN04_v3.000_NTH_gauss_T100.000.xdr')
        #plt.plot(self.grid.wave, mastertau[::-1]*1.1, 'b')#, jan.wave, jan.tau, 'r')
        #plt.show()
        
        self.spectrumtype='tau'
    
    
    
    
    def do_rt(self, N):
        'Sort out the radiative transfer mess.'
        
        if self.spectrumtype != 'tau':
            sys.exit('Wrong type of spectrum input to do_rt, expected tau, found'+self.spectrumtype)
        
        I_back = Spectrum()
        S_func = Spectrum()
        
        I_back.grid.flux = B_nu(I_back.grid.wave, 3e3)
        S_func.grid.flux = B_nu(I_back.grid.wave, self.Temp)
        
        thistau = self.tau*N
        
        I_0 = S_func.grid.flux + np.exp(-thistau)*(I_back.grid.flux-S_func.grid.flux)
        
        masternorm = I_0/I_back.grid.flux
        masteremi = S_func.grid.flux * (1e0 - np.exp(thistau))
        # This turned out not to be wrong, it was just friday afternoon and I wasn't thinking clearly.
               
        self.norm = masternorm
        self.emi = masteremi
        self.spectrumtype = 'RT'
        
    def regrid(self, resolution=ss.resolution, oversample=ss.oversample):
        if self.spectrumtype != 'RT':
            sys.exit('Wrong type of spectrum input to do_rt, expected RT, found'+self.spectrumtype)
        
        #create a wavegrid with the correct resolution:
        outspec = Spectrum(resolution=resolution, oversample=oversample)
        #then we map the calculated spectrum (i.e. that contained in the 'self' on which this method operates
        # onto this new spectral grid, and add the new spectral grid wave and flux to the original spectrum object.
        # this way our output objects will always have a copy of their wavegrid.
        
        #this takes two steps: 1) smooth the output to the right R, 2) resample the smoothed grid onto the output wavegrid.                        
        
        #Steps: make kernel, convolve data with kernel, regrid to smaller grid.
        
        #Jans kernel code: (assuming a constant resolution grid..)
        
        # At this stage, we should have at least a grid of constant resolution
        x = self.grid.wave
        n= np.size(self.grid.wave)

        # Make the basic convolution kernel. This is the same for all but one
        # method. 
        
        w_r = x[n/2]
        sigma = (w_r/float(resolution))/(2e0 * np.sqrt(2e0*np.log(2)))
        mm = np.finfo(float).tiny
        sigma_cut = np.sqrt(-2.0 * np.log(mm))
        idx_gauss = x[(x > w_r - sigma_cut*sigma) & (x <= w_r + sigma_cut*sigma)]
        n_gauss = np.size(idx_gauss)
        idx_gauss = np.arange(n_gauss) - n_gauss/2 + n/2
        new = (x[idx_gauss]-w_r)/sigma
        
        ww = np.exp(-5e-1 * new**2)/(np.sqrt(2.0*np.pi)*sigma)

        # Confirmed kernel matches Jan's. - Step one
        #kernel = readsav('/home/dstock/sf/idl/code/kernel.xdr')
        #plt.plot(ww, 'b', kernel.ww, 'ro' )
        #plt.show()
        
        filtered = signal.fftconvolve(self.norm, ww, 'same')/sum(ww)
        print 'new grid size:', np.size(filtered), np.size(self.grid.wave)
        
        self.smoothednorm = filtered
           
        #Signal in 'filtered' matches the line shape produced by IDL code exactly. Also very fast. - Step two
        
        out_flux = interpolate.griddata(self.grid.wave, filtered, outspec.grid.wave, method='linear')      
        print 'New grid size:', np.size(out_flux), np.size(outspec.grid.wave) 
        
        self.templatewave = outspec.grid.wave
        self.templatenorm= out_flux
        
    def get_tau2(self, molno, isono, ll_name, Temp, regen=False, vturb = ss.vturb):
        #Fixing the problems of the original get_tau method.
        # 1: Chunks with weird sizes and lots of zeros. --- Have removed the leading and trailing zeros.
        # 2: Save/ restore the chunks.
        #  -- How does Jan decide on chunk sizes?
        # 3: work out the necessary overlap region.
        
        self.Temp = Temp
        
        if self.spectrumtype != 'init':
            sys.exit('This spectrum object is not blank.')
        
        #lets make sure that we haven't done this before..
        filename = create_filename(molno, isono, ll_name, 'tau', vturb, Temp=self.Temp)
        
        if os.path.isfile(filename) == False or regen == True:
            
            thislinelist = SpecificLineList(molno, isono, ll_name, regen)
            thislinelist.calc_specifics(Temp)
    
            mastertau = self.grid.flux
            mastertau[:] = 0.0
            
            delta_nud = thislinelist.wave * vturb * 1e5 
        
            # assume that every incoming linelist is sorted.
            #split spectrum into chunks of say 100 lines            
            
            nlines=100
            cycle = 0
            iters = floor( thislinelist.strength.size/nlines )
            
            over = 2e9
            
            #print thislinelist.strength.size
            #print floor(thislinelist.strength.size/nlines)
            #sys.exit()
            a_size=nlines
    
            thisgrid = self.grid.freq[::-1]
            thisfreq = thislinelist.freq
            thiswave = thislinelist.waveum                      
            #first = np.vstack((thisgrid, mastertau))
            #second = np.vstack((thislinelist.freq,thislinelist.strength))
            
            #print first.shape, second.shape
            
            #t1 = time.clock()
            #test = np.concatenate([first, second], axis=1)
            #test2 = np.sort(test, kind="quicksort")
            #t2=time.clock()
            #tic = t2-t1
            
            freqs_cum = np.cumsum(thisfreq)/np.sum(thisfreq)
                        
            # working on finding the breakpoints between the different sets of lines such that we can make a recursive algorithm to split linelists into chunks.
            dx = np.diff(thisfreq)
            length = dx.size

            test = n_max(dx, 10)
           
            print test
            
                       
            freqs = [thisfreq[x[1]] for x in test]
            inds = [x[1] for x in test]
            vals = [x[0] for x in test]
            
            print length
            
            #print freqs
            ends =  sorted(inds)
            
            #The inds are basically the end of each chunk, so create a new array for the beginnings:
         
            starts = [x+1 for x in ends]
            ends.append(length) # the last line has to be the end of the last chunk
            starts.insert(0, 0) # the first chunk starts with the first line..
                        
            sizes = [ x-y for x,y in zip(ends,starts) ]#(ends-starts)
            n_points = [np.size(thisgrid[ (thisgrid < thisfreq[x]) & (thisgrid > thisfreq[y]) ]) for x,y in zip(ends,starts)] 
            #find n points in overall grid between line start[i] and end[i]
            print sizes
            print n_points
            
            arr_sizes = [x*y for x,y in zip(sizes,n_points)]
            print arr_sizes
            
            #print freqs_cum.size, thisfreq[0:-1].size
            
            plt.plot(thisfreq, freqs_cum, '-ob', thisfreq[0:-1], dx/np.max(dx), '-or', freqs, vals/np.max(vals), 'og')
            for x in range(0, len(starts)):
                plt.plot( [thisfreq[starts[x]], thisfreq[ends[x]]], [0.5, 0.5], 'r')
                        
            plt.show()
            
            sys.exit()
            
            for cycle in range(int(iters)+1):
                print 'cycle=', cycle
    
                print 'line number range=', (cycle)*nlines, (cycle+1)*nlines, thislinelist.strength.size
    
                thislooplines_s = thislinelist.strength[(cycle)*nlines: (cycle+1)*nlines]
                thislooplines_w = thislinelist.freq[(cycle)*nlines: (cycle+1)*nlines]
                thisdeltanu = delta_nud[(cycle)*nlines: (cycle+1)*nlines]
                
                a_size =  thislooplines_s.size
                               
                print 'wavelength range=', thislooplines_w[0], thislooplines_w[-1]

                linerange = [thislooplines_w[0], thislooplines_w[-1] ]
                
                print 'linerange=', linerange
                #wavecov = linerange[1] - linerange[0]
                
                #print self.grid.freq[0], self.grid.freq[-1]

                #if cycle==0: sys.exit()
            
                thiswavegrid = thisgrid[(thisgrid > linerange[0]-over) & (thisgrid < linerange[1]+over)]
               
                test = np.concatenate([thiswavegrid, thislooplines_w])
                test2 = np.sort(test, kind='quicksort')
               
                print thiswavegrid - thislooplines_w[50]
                print thislooplines_w

                plt.plot(test2, np.exp(-((test2 - thislooplines_w[50])/thisdeltanu[50])**2), "-bo", thislooplines_w[50], 1e0, "ro" )
                plt.show()
                
                sys.exit()
               
                          
                print 'Nlines, Npoints', thislooplines_s.size, thiswavegrid.size
            
                if thiswavegrid.size == 0:
                    continue
                else:        
                    ff2d = np.vstack([thislooplines_s]*thiswavegrid.size)
                    lines2d = np.vstack([thislooplines_w]*thiswavegrid.size)
                    dnu_2d = np.vstack([thisdeltanu]*thiswavegrid.size)
                    wavegrid2d = np.transpose(np.vstack([thiswavegrid]*a_size))
                    
                    temp = wavegrid2d - lines2d
                    temp = temp/dnu_2d
                    temp = temp**2.0
                    temp = -temp
                    temp = np.exp( temp )
                    
                    mask = ma.where(temp < 500)
                                    
                    temp[mask] = temp[mask] * ff2d[mask]
                    temp[mask] = temp[mask] / dnu_2d[mask]
                    temp[mask] = temp[mask] / np.sqrt(pi)
                    
                    thistau = np.sum(temp, axis=1)
                    
                    
                    mastertau[(thisgrid > linerange[0]-over) & (thisgrid < linerange[1]+over)] += thistau
                
        
            self.tau = mastertau[::-1] #The -1 reverses the array to align it with the um wavelength grid
            outfile = open(filename,'w')
            pickle.dump(mastertau, outfile)
            outfile.close()
        else:
            infile = open(filename,'r') #Can only get to this block if the .pickle file exists. Restore it, its faster.
            mastertau = pickle.load(infile)
            infile.close()
            self.tau = mastertau[::-1] #The -1 reverses the array to align it with the um wavelength grid
            print 'restored tau pickle'           
                
        #jan = readsav('/home/dstock/sf/idl/code/templates/CO2/626/HITRAN04/v3.000/gauss/tau/CO2_626_HITRAN04_v3.000_NTH_gauss_T100.000.xdr')
        #plt.plot(self.grid.wave, mastertau[::-1]*1.1, 'b')#, jan.wave, jan.tau, 'r')
        #plt.show()
        
        self.spectrumtype='tau'


def n_max(arr, n):
    indices = arr.ravel().argsort()[-n:]
    indices = (np.unravel_index(i, arr.shape) for i in indices)
    return [[arr[i], i[0]] for i in indices]
