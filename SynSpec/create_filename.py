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
import numpy as np
import sys

import SynspecSettings as ss
from molDB import molProps
from fileinput import filename
import glob


def create_filename(molno, isono, ll_name, switch, vturb=ss.vturb, want_thermal=0, 
                    wavestart=ss.wavestart, waveend=ss.waveend, profile=ss.profile, 
                    resolution=ss.resolution, oversample=ss.oversample, Temp=1e3, N=1e18, chunkID=-1):
    # Also: chunkID,  N
    # Jan's version used all of these - I am choosing only to include the ones I need so far
    # Not sure about the initialization of Temp here. I am going to make it such that it never gets 
    # referenced though.
    
    #default chunkid=-1 for reasons explained in the 'chunks' case
    
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
        
        if ll_name == 'CDSD':
            #There's more than one file per isotope here !!
            # But its all CO
            # So check
            if molno != 2 or ((isono > 4) or (isono < 1)):
                sys.exit("CDSD Database contains only CO (molno=2, isono=[1,4]), you tried molno: " + str(molno) +', isono: '+str(isono))
            
            moldata = molProps(molno, isono)
            isostring = moldata.isocode

            dir = ss.linedir + ll_name + '/'
            fnames = glob.glob(dir+isostring+'*'+'.cdsd')
            fnames = sorted(fnames)
            
            filename = fnames
            # WARNING: because of the way CDSD works this is actually an array (really a list) of filenames with length between 2 and 6.
        
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
            filename  = totaldir+filename+'_N'+str(np.log10(N))+'_R'+str(resolution)+'.pickle'
    
        
    if switch == 'masterwno': 
        print 'not invented here'
 
    if switch == 'masterwno': 
            print 'not invented here'


    if switch ==  'chunks' or 'chunkinfo' or 'lines': 
        # Create the filename for the chunk[info] files. 
        moldata = molProps(molno, isono)
        temp_dir = ss.templatedir

        moldir = temp_dir +  moldata.molname
        isodir = moldir + '/' + moldata.isocode
        linedir = isodir + '/' + ll_name 

        chunkdir = linedir +'/' + 'linechunks/'

        if (switch == 'chunkinfo' and chunkID >= 0) or (switch == 'chunks' and chunkID < 0):
            #         if format ne 'chunkinfo' XOR keyword_set(chunkID) then stop
            #The IDL code caught this error with an XOR, but what we want to know is: are the inputs consistent?
            # if we want chunkinfo, then chunkid shouldn't have been initialized, likewise
            # if we want chunks, chunkid should be greater than zero.
            sys.exit('Something went wrong with the chunk filenames.' + str(chunkID) + switch)
                
        if vturb == ss.vturb:
            vstring='' 
            vstring2 = ''
        else:
            vstring = '_' + str( vturb)
            vstring2 = vstring+'_'

        if os.path.isdir(chunkdir) == False: #I think you could use an try/catch here, but this is neater.
            os.makedirs(chunkdir)  

        if switch == 'chunkinfo':
            filename = chunkdir + 'chunkinfo' + vstring + '.dat'
        if switch == 'chunks':
            filename = chunkdir+'chunk'+vstring2+str(chunkID)+".lines"

    
    
    if switch ==  'lines': 
        print 'not invented here'
    
    if switch ==  'template_masterwave':
        print 'not invented here'

    print filename

    return filename


#print create_filename(2, 1, 'HITRAN04', 'tau', ss.vturb)

