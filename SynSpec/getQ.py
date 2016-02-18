######################################################################################
#
# Python Reimplementation of Spectra Factory (by Jan Cami)
# 
# This file is a 
# ( ) Translation of an original .pro file (called ...)
# (X) Reimplementation of an original file (called synspec_get_Q.pro)
# ( ) New file
#
# COMMENTS:
# This routine finds the appropriate partition function (Q) for the given input molecule. 
#
#
#
# VERSION HISTORY:
# Created: 09/02/2016 (DJS)
#
######################################################################################''

from molDB import molProps
import numpy as np
import csv
from pylab import *
from scipy.interpolate import interp1d
from scipy.signal.bsplines import cubic
from scipy.interpolate import InterpolatedUnivariateSpline
import time


def getQ(molno, isono, ll_name, Temp):
    moldata = molProps(molno, isono)

    if ll_name == "HITRAN04":
        idx_string = moldata.molname+"_"+moldata.isocode
        #print idx_string
        
        filename = "/home/dstock/sf/data/HITRAN04/GlobalData/parsum.dat"
   
        #now read the file:
        filedata = np.array(list(csv.reader(open(filename), delimiter=' ', skipinitialspace=True)))

        qdata = np.empty(shape=[2, filedata.size-1])
        idx_no = filedata[0].index(idx_string)
           
        for i in range(1,filedata.size):
            cleanedline = [var for var in filedata[i] if var] # if an element is true, whack it in a new array
            qdata[0,i-1] = cleanedline[0] 
            qdata[1,i-1] = cleanedline[idx_no]
        
        if Temp < qdata[0,0]:
            order=1
            s=InterpolatedUnivariateSpline(qdata[0,:], qdata[1,:],k=order)
            if s(Temp) < 0.0:
                sys.exit('Partition function extrapolation to low temperature failed. (Q < 0)')  
            return s(Temp)
        elif Temp > qdata[0,-1]:
            sys.exit('Extrapolation to high temperatures not implemented.')
        else:
            f = interp1d(qdata[0,:], qdata[1,:])#, kind='cubic') #Cubic makes this slooooooow.
            return f(Temp)
  
        
        
        
#test = getQ(1,1,'HITRAN04',25.0)
#print test