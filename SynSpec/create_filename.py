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
# Added tau file paths: 19/02/2016
#
######################################################################################

import os

import SynspecSettings as ss
from molDB import molProps



def create_filename(molno, isono, ll_name, switch, vturb=ss.vturb, want_thermal=0, 
                    wavestart=ss.wavestart, waveend=ss.waveend, profile=ss.profile, 
                    resolution=ss.resolution, oversample=ss.oversample, Temp=1e3):
    # Also: chunkID,  N
    # Jan's version used all of these - I am choosing only to include the ones I need so far
    # Not sure about the initialization of Temp here. I am going to make it such that it never gets 
    # referenced though.
    
    
    
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

        
    if switch == 'tau' or switch == 'template':
        moldata = molProps(molno, isono)
        
        #Preserving Jan's directory structure
        temp_dir = ss.templatedir
        
        moldir = temp_dir +  moldata.molname # no slash here because tempdir already has a trailing slash.
        isodir = moldir + '/' + moldata.isocode
        linedir = isodir + '/' + ll_name 
        vdir = linedir + '/' + 'v'+str(vturb)
        profdir = vdir + '/' + profile

        filename = moldata.molname + '_' + moldata.isocode + '_' + ll_name + '_' + 'v'+str(vturb)
        if want_thermal == 1:
            filename = filename + '_TH_'
        else:
            filename = filename + '_NTH_'
        
        filename  = filename + profile + '_' + 'T'+str(Temp)

        if switch == 'tau':
            totaldir = profdir + '/' + 'tau/'
            if os.path.isdir(totaldir) == False: #I think you could use an try/catch here, but this is neater.
                os.makedirs(totaldir)            
            filename  = totaldir+filename+'.pickle'
        
        if switch == 'template':
            totaldir = profdir + '/' + 'templates/'
            if os.path.isdir(totaldir) == False: #I think you could use an try/catch here, but this is neater.
                os.makedirs(totaldir)            
            filename  = totaldir+filename+'.pickle'
    
        
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
    
    if switch ==  'template_masterwave':
        print 'not invented here'

    print filename

    return filename


#print create_filename(2, 1, 'HITRAN04', 'tau', ss.vturb)

