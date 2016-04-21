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
######################################################################################




class chunk():
    def __init__(self, wave, freq, strength, lineinds, gridinds, nchunks):
        self.nchunks = nchunks
        self.wave = wave
        self.freq = freq
        self.strength = strength
        self.lineinds = lineinds
        self.gridinds = gridinds
