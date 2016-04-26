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
# This file defines the final template object containing all of the calculated quantities.
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

import os.path
import pickle
from scipy.io import readsav

import SynspecSettings as ss
from create_filename import create_filename
import matplotlib.pyplot as plt
from linelist import SpecificLineList

from spectrum import Spectrum

  
class Template:
    def __init__(self, molno, isono, Temp, N, regen=False, resolution=ss.resolution, oversample=ss.oversample, ll_name='HITRAN04', v_turb=ss.vturb):
        #First of all lets see if we have made the requested template before
        filename = create_filename(molno, isono, ll_name, 'template', vturb=v_turb, Temp=Temp, N=N, resolution=resolution)
        
        if os.path.isfile(filename) == False or regen == True:
            #Stuff to (re)generate template
            data = Spectrum()
            
            #check if this spectrum has been chunked
            
            #chunkinfo filename:
            fname = create_filename(molno, isono, ll_name, "chunkinfo")
            
            if not os.path.isfile(fname) or regen == 1:
                thislinelist = SpecificLineList()
                thislinelist.readlines(molno, isono, ll_name, regen)
                           
            data.get_tau_chunks(molno, isono, ll_name, Temp, regen)
            
            data.do_rt(N)
            data.regrid(resolution=resolution, oversample=oversample)    
            self.data = data
            
            outfile = open(filename,'w')
            pickle.dump(data, outfile)
            outfile.close()
            print 'saved template pickle'
            
        else:
            infile = open(filename,'r') #Can only get to this block if the .pickle file exists. Restore it, its faster.
            self.data = pickle.load(infile)
            infile.close()
            print 'restored template pickle'
     



template = Template(2,1,1010.0, 10**19.0, ll_name='CDSD', regen=True, resolution=120.0)


    
#jan = readsav('/home/dstock/sf/idl/code/templates/CO2/626/HITRAN04/v3.000/gauss/templates/R120.00/CO2_626_HITRAN04_v3.000_NTH_gauss_T1000.000_N18.500_R120.00_O2.xdr')
#wave = readsav('/home/dstock/sf/idl/templates/masters/wave_R120.00_O2_umstart_0.350000_umend_3000.000000.xdr')

#jan = readsav('/home/dstock/sf/idl/code/templates/CO2/626/HITRAN04/v3.000/gauss/tau/CO2_626_HITRAN04_v3.000_NTH_gauss_T1000.000.xdr')

#plt.plot(template.data.grid.wave, template.data.tau, 'r')#, jan.wave, jan.tau, 'b')#, data.grid.wave, data.smoothednorm, 'g')

plt.plot(template.data.templatewave, template.data.templatenorm)

#plt.plot(template.data.tau/jan.tau)
plt.xlim(0,500)

plt.show()    
    