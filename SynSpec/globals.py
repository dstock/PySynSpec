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
# This file contains the `globals' classes, which basically fetches information
# relevant to entire linelists. Mainly this means partition functions.
#
#
# VERSION HISTORY:
# Created: 08/02/2016 (DJS)
#
######################################################################################

import numpy as np
import csv
import sys

from molDB import molProps

class HT04_globals():
    def __init__(self, molno, isono):
        
        # Replace this filename with a call to create filename
        filename = '/home/dstock/sf/data/HITRAN04/GlobalData/molparam.txt'
        
        moldata = molProps(molno, isono)
        
        # First read the file containing all the global Values
        with open(filename) as csvfile:
            reader = csv.reader(csvfile)
            right_molno = False
            molmissing = True
            isomissing = True
            i=0
            for row in reader:
                if i != 0:
                    molstring = row[0][0:8]
                    if molstring != "        ":
                        temp_molno = row[0][8:10]
                                          
                        if int(temp_molno.strip(")")) == molno:
                            right_molno=True
                            molmissing=False
                        else:
                            right_molno=False #Important to reset this
                   
                    else:
                        isostring = row[0][8:12]
                        if int(isostring) == int(moldata.isocode) and right_molno:
                            self.molno = molno
                            self.isono = isono
                            self.abund = float(row[0][14:25]) # figure out range
                            self.Q296 = float(row[0][29:39]) # figure out range
                            self.gj = int(row[0][43:44])
                            self.mm = float(row[0][49:59])
                            isomissing=False
                            #need to copy Jans manipulation of these variables. Presumably T is involved.
                            #load up with Selfs, we found the one!
                            break
                    
                                     
                i+=1
            if molmissing or isomissing: sys.exit("""Molecule or Isotopolog not found in molparam.txt. If you ever see this,
                                                     something is very wrong as the combination specified DOES exist in
                                                     the molecular database list moltranslation.txt""") # In practice you should never, ever see this message.

#test = HT04_globals(2,9,1000) # some test cases
#         181  1.99983E-03    1.7511E+02    1     20.014811         
#print test.abund, test.Q296, test.gj, test.mm
        