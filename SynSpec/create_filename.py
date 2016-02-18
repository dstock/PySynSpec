######################################################################################
#
# Python Reimplementation of Spectra Factory (by Jan Cami)
# 
# This file is a 
# ( ) Translation of an original .pro file (called ...)
# (X) Reimplementation of an original file (called synspec_create_filename.pro)
# ( ) New file
#
# COMMENTS:
# 
# General fucntion to return the proper PATH + filename for a given
# calculation/template to be found/saved. 
#
# VERSION HISTORY:
# Created: 05/02/2016 (DJS)
#
######################################################################################


import SynspecSettings as ss
import sys




def create_filename(molno, isono, ll_name, switch, vturb=ss.vturb, want_thermal=0, 
                    wavestart=ss.wavestart, waveend=ss.waveend, profile=ss.profile, 
                    resolution=ss.resolution, oversample=ss.oversample ):
    # Also: chunkID, Temp, N
    # Jan's version used all of these - I am choosing only to include the ones I need so far
    
    
    
    if switch == 'linelist':
        if ll_name == 'HITRAN04':
            molstring = '{:0=2d}'.format(molno)
            isostring = '{:0=2d}'.format(isono)
            fileroot = molstring + "_" + isostring + "_hit04.par"
            filename = ss.linedir + ll_name + '/' + fileroot
        
        if ll_name == 'HITRAN08':
            print 'not implemented yet'

        if ll_name == 'HITRAN12':
            print 'not implemented yet'


        
    if switch == 'masterwno': 
        print 'not invented here'
 
    if switch == 'masterwno': 
            print 'not invented here'

    if switch ==  'chunkinfo': 
        print 'not invented here'
    
    if switch ==  'chunks': 
        print 'not invented here'
    
    if switch ==  'lines': 
        print 'not invented here'
    
    if switch == 'tau':
        print 'not invented here'
    
    if switch ==  'template': 
        print 'not invented here'
    
    if switch ==  'template_masterwave':
        print 'not invented here'

    print filename

    return filename
