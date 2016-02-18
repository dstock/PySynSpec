'''
Created on 2016-02-08

@author: dstock
'''

from sfClasses import Spectrum
import time
import platform
print platform.architecture()

data = Spectrum()
tick = time.clock()
data.get_tau(2,1,'HITRAN04',1000)
tock = time.clock()

print tock-tick, '<='


