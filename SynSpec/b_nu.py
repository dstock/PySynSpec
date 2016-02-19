#########################################################################################
#
# PySpecFac - A reimplementation of IDL SpectraFactor (by Jan Cami)
# 
# This file is a 
# (X) Translation of an original .pro file (called b_nu.pro)
# ( ) Reimplementation of an original file (called ...)
# ( ) New file
#
# COMMENTS:
# 
# 28/04/1999 : adapted to implement correction as in Engelke, AJ 104, 1248 (1992) (JC)
# 09/05/1999 : removed bug that returned the T-array when using Engelke option (JC)
# 02/02/2016 : Converted to Python (DS)
#
# lambda assumed in micron !
#
#########################################################################################


def B_nu( wave, T, engelke=0 ):
    import astropy.constants as apc 
    import numpy as np
    import sys
    #Importing functions like this seems to be uncommon? Might have speed penalty?
    
        #Check that wave and T are both floats:
    stat1 = wave.dtype == 'float64'   
    print wave.dtype, stat1
    #Weird, because I want wave to be an numpy array its type is in a string
    stat2 = type(T) == float
    print type(T), stat2
    if stat2 == False or stat2 == False:
        sys.exit('Types not correct in function B_nu') # Kill it if wrong

    if engelke == 1:
            T_eff = T
            T = 0.738 * T_eff * (1.0 + (79450.0/(wave * T_eff)))**0.182
    
    freq =  apc.c.cgs.value * 1e4 / wave

    e_factor = 1.0/( np.exp( (apc.h.cgs.value * freq )/( apc.k_B.cgs.value * T ) ) - 1.0)
    
    b_nu = 2.0 * apc.h.cgs.value * freq**3 * e_factor / (apc.c.cgs.value)**2

    return b_nu
