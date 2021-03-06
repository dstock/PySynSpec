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
# This file defines the linelist object, along with readers for HITRAN04 and eventually others.
#
#
#
# VERSION HISTORY:
# Created: 04/02/2016 (DJS)
# V0.001: Can create Tau Profiles for HITRAN04 linelists 18/02/2016
# V0.01: Can make output templates whoch match IDL versions for arbitrary
#        resolutions, v_turbs, and linelists from HITRAN04 database.
#        The line chunking routine in the tau calculation is very basic
#        and needs updating and steamlining. In particular the leading 
#        chunk of each wavelength grid is an arbitrarily huge array of 
#        zeros which we're still handling.However it seems to be correct.
#        Extensively tested using CO 626 (molno=2, isono=1).
#        
#        Split various classes into separate files (26/02/16, DJS)
#
######################################################################################

import csv
import astropy.constants as apc
import numpy as np
import sys
import os.path
import pickle
from scipy.io import readsav

from makegrid import grid
from line import Line
from globals import HT04_globals
import SynspecSettings as ss
from create_filename import create_filename
from getQ import getQ
from find_chunks import find_chunks
from chardet.latin1prober import FREQ_CAT_NUM
from chunk import chunk
from SynspecSettings import chunkdir


c2 = (apc.h.cgs.value*apc.c.cgs.value)/apc.k_B.cgs.value
 

class LineList:
    """A basic line list constructor"""
    def __init__(self, lines):  # need to add options
        #This implicitly expects `lines' to be an array of line objects.
        self.lines = lines
        spec = np.array([x.spec for x in lines])
        iso  = np.array([x.iso for x in lines])
        
        self.freq = np.array([x.freq for x in lines])
        self.wno = np.array([x.wave for x in lines])
        self.waveum = np.array([x.waveum for x in lines])
        self.eA =  np.array([x.eA for x in lines])
        self.epp = np.array([x.epp for x in lines])
        self.strength = np.array([x.strength for x in lines]) 
    
        if np.var(spec) == 0.0:
            self.spec = np.median(spec)
        else:
            sys.exit('Trying to construct linelist from multiple species') # Kill it if wrong

        if np.var(iso) == 0.0:
            self.iso = np.median(iso)
        else:
                sys.exit('Trying to construct linelist from multiple isotopologues') # Kill it if wrong

    def addline(self):
        print 'To be implemented later if necessary'
    
    def removeline(self):
        print 'To be implemented later if necessary'

    def ll_reverse(self):
        'Linelists normally presented as sorted in wavenumber. We want sorted by um.'
        self.freq = self.freq[::-1]
        self.wno = self.wno[::-1]
        self.waveum = self.waveum[::-1]
        self.eA = self.eA[::-1]
        self.epp = self.epp[::-1]
        self.strength = self.strength[::-1] 
    
class SpecificLineList(LineList):
    """Container for methods fetching hitran04 line lists"""
    
    #Comments:
    #This needs to be refactored such that __init__ initializes all the variables
    # and then we have class methods which load up those initial variables using class methods.
    
    #We can write the new init method and leave the old one in place as the second definition
    # of init rebinds the function.
    
    def __init__(self):
        self.ll_name = ""
        self.specs_calced = False
        self.wno = []
        self.waveum = []
        self.freq = []
        self.epp = []
        self.eA = []
        self.iso = []
        self.spec = 0
        self.iso = 0
    
    def readlines(self, molno, isono, ll_name, regen=False):
        
        print molno, isono
        self.ll_name = ll_name
        props = HT04_globals(molno, isono)
        self.specs_calced = False

        
        if ll_name == 'HITRAN04':
            filename = create_filename(molno, isono, ll_name, 'linelist')
            #filename = 'filereadtest/01_01_hit04.par'
            
            if os.path.isfile(filename) == False:
                sys.exit('File not found:'+filename)
            
            if os.path.isfile(filename+'.pickle') == False or regen==True:  # Here I am checking whether we already have an object for this linelist.
                #read in the file - first find number of rows
                with open(filename) as csvfile:
                    reader = csv.reader(csvfile)
                    i=0
                    for row in reader:
                        #print i
                        i+=1
        
                lines = np.arange(i,dtype=object)
        
                #now loop through and fill the array with lines
                with open(filename) as csvfile:
                    reader = csv.reader(csvfile)
                    i=0
                    for row in reader:
                        #print i
                        #print('{:2d}{:1d}{:12.6f}{:10.3e}{:10.3e}{:5.4f}{:5.4f}{:10.4f}{:4.2f}{:8.6f}{:15c}{:15c}{:15c}{:15c}{:6d}{:12d}{:1c}{:7.1f}{:7.1f}'.format(row[0]))
                        data = row[0]
                        thisLine = Line( np.uint(data[0:2]), np.uint(data[2:3]), np.float64(data[3:15]), np.float64(data[25:35]), strength=np.float64(data[15:25]), epp=np.float64(data[45:55]))
                        lines[i] = thisLine
                        i+=1

                outfile = open(filename+'.pickle','w')
                pickle.dump(lines, outfile)
                outfile.close()
            else:
                infile = open(filename+'.pickle','r') #Can only get to this block if the .pickle file exists. Restore it, its faster.
                lines = pickle.load(infile)
                infile.close()


            #
            #
            #
            self.lines = lines[::-1] # Trying to call attention to this line, where I explicitly reverse the linelist to make it go with increasing lambda.
            # 
            #
            #
            self.wno = np.array([x.wno for x in self.lines]) # currently this is actually wno
            self.waveum = np.array([x.waveum for x in self.lines])
            self.epp = np.array([x.epp for x in self.lines])
            self.eA = np.array([x.eA for x in self.lines])
            self.freq = np.array([x.freq for x in self.lines])
            E_temp = -1.0*self.epp * c2 / 296.0
            w_temp = -1.0*self.wno * c2 / 296.0
            self.strength = np.array([x.strength for x in self.lines]) * (props.Q296/props.abund) / (np.exp(E_temp) * (1.0-np.exp(w_temp)) )
            #print np.array([x.strength for x in self.lines]) 
            spec = np.array([x.spec for x in self.lines])
            iso = np.array([x.iso for x in self.lines])
        
        
            if np.var(spec) == 0.0:
                self.spec = int(np.median(spec))
            else:
                sys.exit('Trying to construct linelist from multiple species') # Kill it if wrong

            if np.var(iso) == 0.0:
                self.iso = int(np.median(iso))
            else:
                sys.exit('Trying to construct linelist from multiple isotopologues') # Kill it if wrong
      
            self.create_chunks()
      
        elif ll_name == 'CDSD':
            filename = create_filename(molno, isono, ll_name, 'linelist')
            # this filename is really a list of n filenames where 2 <= n <= 5. So we have to loop over them and construct the linelist
            
            # Filename is sorted coming from create_filename
            # In order to get things in increasing wavelength space we run this backwards
            filename = filename[::-1] 
            
            #
            #
            # We have to make some choices here about how to do this. We can either load in each file and add them to one overall linelist
            # or we could read in each file, make that into a separate linelist, chunk that linelist, then read in the next file, chunk that
            # and keep track of the overall number of chunks (i.e. make sure the chunk indices make sense.
            # The second method is probably the most extensible based on the idea that we will eventually use line lists bigger than even CDSD.
            # For the sake of providing a template method for larger linelists, we will do the second way even though we might get away with the
            # first in this instance.
            #
            init_no = 0
            
            for thisfile in filename: # This is a pythonic loop over each element of the filename list
                # !!!!!!!! In this loop 'thisfile' takes the place of the usual 'filename' !!!!!!!!
                
                # Read in the file
                # Process it into a sub-linelist
                # Chunk the sublinelist
                # Pass the number of chunks in this sub-linelist to the next iteration so the chunk numbers are sequential
                # Update all chunks with final value for nchunks.
                print molno, isono
                Qref = getQ(molno, isono, ll_name, 1000.0)
                print Qref
                
                #First of all, check THIS file exists
                if os.path.isfile(thisfile) == False:
                    sys.exit('File not found:'+file)
                
                print thisfile
                
                # This block is the same as in the HITRAN04 bit. In principle it could be abstracted into its own method.
                # However at the time of writing I'm not sure how much it needs to be change so I have done it this way
                # Possible todo: abstract this into 'filereader' function.
                
                #For each linelist we check if we already read and pickled it.
                if os.path.isfile(thisfile+'.pickle') == False or regen==True:  # Here I am checking whether we already have an object for this linelist.
                    #read in the file - first find number of rows
                    with open(thisfile) as csvfile:
                        reader = csv.reader(csvfile)
                        i=0
                        for row in reader:
                            #print i
                            i+=1
            
                    lines = np.arange(i,dtype=object)
                    print i
            
                    #now loop through and fill the array with lines
                    with open(thisfile) as csvfile:
                        reader = csv.reader(csvfile)
                        i=0
                        for row in reader:
                            data = row[0]
                            thisLine = Line( molno, isono, np.float64(data[4:15]), np.float64('NaN'), strength=np.float64(data[16:25]), epp=np.float64(data[45:55]))
                            lines[i] = thisLine
                            i+=1
                            #in future for even bigger linelists, I think you could break here, process what you have so far and continue. In practice
                            # my 1.8ghz/12gb ram laptop can process 700k+ linelists without exploding or swapping, so I suspect it would only be relevant for the
                            # really huge linelists.
                            
                            
                    outfile = open(thisfile+'.pickle','w')
                    pickle.dump(lines, outfile)
                    outfile.close()
                else:
                    infile = open(thisfile+'.pickle','r') #Can only get to this block if the .pickle file exists. Restore it, its faster.
                    lines = pickle.load(infile)
                    infile.close()                
                

                #
                #
                self.lines = lines[::-1] # Trying to call attention to this line, where I explicitly reverse the linelist to make it go with increasing lambda.
                # 
                #
                #
                self.wno = np.array([x.wno for x in self.lines]) # currently this is actually wno
                self.waveum = np.array([x.waveum for x in self.lines])
                self.epp = np.array([x.epp for x in self.lines])
                self.eA = np.array([x.eA for x in self.lines])
                self.freq = np.array([x.freq for x in self.lines])
                E_temp = -1.0*self.epp * c2 / 1000.0
                w_temp = -1.0*self.wno * c2 / 1000.0
                self.strength = np.array([x.strength for x in self.lines]) * (Qref/props.abund) / (np.exp(E_temp) * (1.0-np.exp(w_temp)) )
                #print np.array([x.strength for x in self.lines]) 
                spec = np.array([x.spec for x in self.lines])
                iso = np.array([x.iso for x in self.lines])                
                
                if np.var(spec) == 0.0:
                    self.spec = int(np.median(spec))
                else:
                    sys.exit('Trying to construct linelist from multiple species') # Kill it if wrong

                if np.var(iso) == 0.0:
                    self.iso = int(np.median(iso))
                else:
                    sys.exit('Trying to construct linelist from multiple isotopologues') # Kill it if wrong
                
                
                
                self.create_chunks(init_no=init_no)
                
                init_no = init_no + self.nchunks
                print init_no

            
    
            #needs to be a loop in here over all chunks which updates the nchunks attribute of each chunk.
            # handily this should be contained within init_no after the last file has been chunked.
            
            chunkstoupdate = True
            index = 0
            
            while chunkstoupdate:
                thischunkfilename = create_filename(molno, isono, ll_name, "chunks", chunkID=index)
                infile = open(thischunkfilename+'.pickle','r')
                thischunk = pickle.load(infile)
                infile.close()
                
                thischunk.nchunks = init_no
                
                outfile = open(thischunkfilename+'.pickle','w')
                pickle.dump(thischunk, outfile)
                outfile.close()
                              
                index = index+1
                if index == init_no -1:
                    chunkstoupdate = False
            
            
            print "done reading CDSD Files and Creating Chunks"
            
        else:
            print 'other Linelists Not Implemented'             
    
    def calc_specifics(self, Temp):
        """A separate method to calculate the specific line list properties based on an input T."""
        #Keep this method for when the chunk object gets merged back into linelist.
        
        if self.specs_calced == False:
            #make sure we don't inadvertently try and do this twice
            if self.ll_name == 'HITRAN04':
                self.Temp = Temp
                self.specs_calced = True
                #lets make sure the relevant temperature is now carried around with the linelist.                
                
                props = HT04_globals(self.spec, self.iso)
                
                if Temp == 296.0 and self.ll_name == 'HITRAN04':
                    Q=props.Q296
                else:
                    Q=getQ(self.spec, self.iso, self.ll_name, Temp)    
                
     
                E_temp = -1.0 * self.epp * c2 / Temp
                #print E_temp
                w_temp = -1.0 * self.wno * c2 / Temp
                #print w_temp
                self.strength = self.strength * (props.abund/ Q) * (np.exp(E_temp) * (1.0-np.exp(w_temp))) * apc.c.cgs.value
                #I have no idea why Jan multiplies by C here, but he does, so lets copy it.
                # Note from 27/02/16: the C corrects the units of strength to energy density.
                
                #strengths_jan = readsav('/home/dstock/sf/idl/code/ff.xdr')
                
                #print "My Calcs:", self.strength
                #print 'My epp', self.epp
                #print "Jan's Calcs:", strengths_jan.ff

                #print self.strength[0:12]/strengths_jan.ff

                #print strengths_jan.ff[0]/self.strength[0]

    def create_chunks(self, init_no=0):
        # init_no is intended to function as the index of the first chunk
        # This allows us to independently chunk linelists that are stored in
        # multiple files.
        
        thisfreq = self.freq
               
        chunks = find_chunks(thisfreq, plot=0)
                
        starts = chunks[0]
        ends = chunks[1]
        
        nchunks = len(starts)
        
        self.nchunks = nchunks
        
        print "*****************************"
        print starts
        print ends
        print len(starts)
        print len(ends)
        print "*****************************"
         
        
        newgrid = grid()
        thisgrid = newgrid.waveum#[::-1]
        
        sumfilename = create_filename(self.spec, self.iso, self.ll_name, "chunkinfo")

        summary = open(sumfilename, 'w')
        
        #headerstring = str(format='(A4,4A15,A25,2A15,A25)', "ID", "idx1", "idx2", "wnogridstart", "wnogridend", "ngridpoints", "wnolinestart", "wnolineend", "nlines")
        
        h1 = "ID"
        h2 = "idx1"
        h3 = "idx2"
        h4 = "wnogridstart"
        h5 = "wnogridend"
        h6 = "ngridpoints"
        h7 = "wnolinestart"
        h8 = "wnolineend"
        h9 = "nlines"
        
        headerstring =  "{:>4} {:>15} {:>15} {:>15} {:>15} {:>25} {:>15} {:>15} {:>25}".format(h1, h2, h3, h4, h5, h6, h7, h8, h9)

                
        summary.write(headerstring+"\n")
        
        for i in range(0, len(starts)):
            print i
            fname = create_filename(self.spec, self.iso, self.ll_name, "chunks", chunkID=i+init_no)
            chunk = self.extract_chunks(starts[i], ends[i], thisgrid, nchunks)
            outfile = open(fname+'.pickle','w')
            pickle.dump(chunk, outfile)
            outfile.close()
        
            if chunk.outwithgrid == 0:
                summarystring = "{:>4} {:>15} {:>15} {:>15} {:>15} {:>25} {:>15} {:>15} {:>25}".format(i+init_no, chunk.gridinds[0], chunk.gridinds[1], 
                            newgrid.wno[chunk.gridinds[0]], newgrid.wno[chunk.gridinds[1]], (chunk.gridinds[1]- chunk.gridinds[0]), 
                            self.waveum[starts[i]], self.waveum[ends[i]], ends[i]-starts[i] )+"\n"
            else:
                summarystring = "{:>4} {:>15} {:>15} {:>15} {:>15} {:>25} {:>15} {:>15} {:>25}".format(i+init_no, "NA", "NA", 
                            "NA", "NA", "NA", 
                            self.waveum[starts[i]], self.waveum[ends[i]], ends[i]-starts[i] )+"\n"
            summary.write(summarystring)
        
        
        summary.close()
        
        #This chunking routine works differently from the original, in that it now stores the chunks in increasing wavelength space. I.e.
        # chunk one contains the lowest wavelength lines.



    def extract_chunks(self, start, end, grid, nchunks):
        
        #Need to make something here which makes `grid' into just a section of the wavegrid which 
        # corresponds to the given chunk of lines.
        mask = np.ma.where( (grid < self.waveum[end]) & (grid > self.waveum[start]) )
        
        if len(mask[0]) == 0:
            outwithgrid = True
            gridinds = [float('NaN'),float('NaN')]
        else:    
            gridinds = [mask[0][0], mask[0][-1]]
            outwithgrid = False
        
        grid  = grid[mask] 
        #print grid, grid.size
        return chunk(self.waveum[start:end], self.freq[start:end], self.strength[start:end], [start, end],
                      gridinds, nchunks, self.ll_name, self.spec, self.iso, outwithgrid, self.epp[start:end])
         
    



#testing
#list1 = SpecificLineList(2, 1, 'HITRAN04', True)
#list1.calc_specifics(1000)
#list1.create_chunks()        

