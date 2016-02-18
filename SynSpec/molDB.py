######################################################################################
#
# Python Reimplementation of Spectra Factory (by Jan Cami)
# 
# This file is a 
# ( ) Translation of an original .pro file (called ...)
# (X) Reimplementation of an original file (called synspec_get_dbinfo.pro)
# ( ) New file
#
# COMMENTS:
# This routine performs the database lookup between our internal numbering system and 
# the overall properties of each molecule found in the moltran file.
#
#
# VERSION HISTORY:
# Created: 09/02/2016 (DJS)
#
######################################################################################

import csv
import numpy as np
import sys

class molProps():
    """ This class explictly pulls out the properties of ONE molecule. I have kept this separate from
    the following such that you have the option of pulling out the full moltranslation database using
    molDB or just finding the properties of the molecule we're currently working on."""
    def __init__(self, molno, isono):
        data = molDB()
        molmissing = True
        isomissing = True
        for i in range(1,len(data.dbdata['spec'])):
            #print  i, data.dbdata['spec'][i], molno
            if np.uint(data.dbdata['spec'][i]) == molno:
                molmissing = False
                #print i, data.dbdata['isono'][i], isono
                if np.uint(data.dbdata['isono'][i]) == isono:
                    isomissing = False
                    self.spec = molno
                    self.molname = data.dbdata['molname'][i]
                    self.isono = isono
                    self.isocode = data.dbdata['isocode'][i]
                    self.isoname = data.dbdata['isoname'][i]
                    self.LL = data.dbdata['LL'][i]
                    break
        
        if molmissing or isomissing: sys.exit("Combination of molecule number and isotopologue number not found in MolTran.txt. ")       
                # The way I have done this very explicitly assumes that there will be precisely ONE match.

class molDB():
    """ This class loads the entire molecule translation database."""
    dbdata = {'spec':[], 'molname':[], 'isono':[], 'isocode':[], 'isoname':[], 'LL':[]}
    
    filename = "/home/dstock/sf/data/general/moltranslation.txt"
    
    filedata = np.array(list(csv.reader(open(filename), dialect='excel-tab'))) # when the time comes try skipinitialspace=True
    #print  type(filedata), filedata.shape, filedata[2]
    
    for i in range(filedata.size):
        cleanedline = [var for var in filedata[i] if var] # if an element is true, whack it in a new array
        # this gets rid of the empty tab fields in the array.
        #print cleanedline        
        dbdata['spec'].append(cleanedline[0]) 
        dbdata['molname'].append(cleanedline[1])
        dbdata['isono'].append(cleanedline[2])
        dbdata['isocode'].append(cleanedline[3])
        dbdata['isoname'].append(cleanedline[4])
        dbdata['LL'].append(cleanedline[5])
 
    #print type(dbdata)
    
#test = molDB()
#print test.dbdata['molname'][1]

#test2= molProps(3,6)
#print test2.spec, test2.molname, test2.isono, test2.isocode, test2.isoname,  test2.LL
#test2.
