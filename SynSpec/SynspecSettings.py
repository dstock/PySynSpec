######################################################################################
#
# Python Reimplementation of Spectra Factory (by Jan Cami)
# 
# This file is a 
# (X) Translation of an original .pro file (called synspec_startup.dat)
# ( ) Reimplementation of an original file (called ...)
# ( ) New file
#
# COMMENTS:
# This is a container for all the program defaults. 
#
#
#
# VERSION HISTORY:
#
# Created: 04/02/2016 - DJS
#
######################################################################################


    ## make this read a file to get a list of settings and defaults
profile      = 'gauss'       # Default profile for optical depth
wavestart    = 0.35          # Min wave for full range (in microns)
waveend      = 3000.0        # Max wave for full range (in microns)
maxexp       = 5e2           # Max value for exponential calculations.
overwrite    = 0	         # 0=no, 1=template o/w, 2=tau overwrite 
vturb        = 3e0           # v_turb (for splitting linelist in chunks)
maxmem       = 1e8           # Maxmem in bytes (so 1d9 is 1 Gbyte)
extraline    = 1e2	         # extra space for significant lines (see synspec_get_lines.pro)
resolution   = 3e2	         # Default instrumental resolution
oversample   = 2	         # Default instrumental oversampling rate. 
exp2         = 17            # Default power of two for fft operations
exp2_con     = 12	         # Default power of two for convol operations
smoothmethod = 'blk_con'     # Default smoothing method. 
linedir      = '/home/dstock/sf/data/'  # Location of linelists
moldir       = ''
isodir       = ''
vdir         = ''
profdir      = ''
taudir       = ''
templatedir  = ''
resdir       = ''
masterdir    = ''
chunkdir     = ''

