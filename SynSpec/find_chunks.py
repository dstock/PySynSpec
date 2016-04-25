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
# The idea here is to work out the location of the best chunk sizes recursively.
#
#
# VERSION HISTORY:
# Created: 23/03/2016 (DJS)
#
######################################################################################

import numpy as np
import matplotlib.pyplot as plt
import sys

from makegrid import grid


""" These functions take an arbitrary linelist and find the start and end indices of the most appropriate chunks to take
for later processing. The first function, find_chunks, looks at the whole linelist and used the diff function to make n
chunks which contain minimal blank wavelength space, i.e. it finds the biggest line clusters in wavelength space.
The second function, takes the chunks from the first method which exceed a certain size, and splits them in half until their 
n_lines * n_wavepoints value is below a certain threshold. This is done recursively, i.e. split_chunks calls itself!
"""


def find_chunks(linefreqs, plot=0):
        
        
    thisgrid = grid() # With no arguments this returns the default wavegrid. Need to preserve that behaviour.
    thisgrid = thisgrid.freq#[::-1]
    
    freqs_cum = np.cumsum(linefreqs)/np.sum(linefreqs)
    

    # working on finding the breakpoints between the different sets of lines such that we can make a recursive algorithm to split linelists into chunks.
    dx = -np.diff(linefreqs) # This has to be negative because our cumulative distribution is actually the wrong way round.
    length = dx.size

    test = n_max(dx, 20)
                   
    freqs = [linefreqs[x[1]] for x in test]
    inds = [x[1] for x in test]
    vals = [x[0] for x in test]

    if plot==1:    
        plt.plot(sorted(vals)[::-1])
        plt.show()
        print sorted(vals)[::-1]/np.max(vals)


    # If you plot the distribution above, you see that the size of the gaps goes down to less than 10% of the largest gap within 20.
    # Therefore we only include the gaps bigger than 10% of the size of the biggest gap in our initial cuts.

    inds = [x for x,y in zip(inds, vals) if y/np.max(vals) > 0.1]
    freqs = [x for x,y in zip(freqs, vals) if y/np.max(vals) > 0.1]


    #print length
    
    #print freqs
    ends =  sorted(inds)
    
    #The inds are basically the end of each chunk, so create a new array for the beginnings:
 
    starts = [x+1 for x in ends]

    
    starts.insert(0, 0) # the first chunk starts with the first line..
    ends.append(length) # the last line has to be the end of the last chunk
                
    sizes = [ x-y for x,y in zip(ends,starts) ]#(ends-starts)
    n_points = [np.size(thisgrid[ (thisgrid > linefreqs[x]) & (thisgrid < linefreqs[y]) ]) for x,y in zip(ends,starts)] 
    #find n points in overall grid between line start[i] and end[i]
    # > and < appear back to front because the grid runs the wrong way.

    
    chunk_sizes = [x*y for x,y in zip(sizes,n_points)]
    #print arr_sizes
    
    print chunk_sizes
    
    thresh = 1250000
    
    need_subdiv = [x > thresh for x in chunk_sizes]
    print need_subdiv
    
            
    # based on some kind of size criteria here we want to split up the chunks that are too big.
    
    print starts
    print ends
    
    oldstarts = starts
    oldends = ends
    
    for i in range(0,len(need_subdiv)):
        if need_subdiv[i] == True:
            subchunks = split_chunks(linefreqs, oldstarts[i], oldends[i], thresh) 
            #Then something here to reintegrate the chunks into overall chunk array
            print subchunks[0], starts
            starts = starts + [i for i in subchunks[0] if i not in starts]
            
            print subchunks[1], ends
            ends = ends + [i for i in subchunks[1] if i not in ends]
        
    
    starts = sorted(starts)
    ends = sorted(ends)
            
    if plot == 1:    
        plt.plot(linefreqs, freqs_cum, '-ob', linefreqs[0:-1], dx/np.max(dx), '-or', freqs, vals/np.max(vals), 'og')
        for x in range(0, len(starts)):
            plt.plot( [linefreqs[starts[x]], linefreqs[ends[x]]], [0.5, 0.5], '-or')
        plt.show()
    
    
    return [starts, ends]
                


    
def split_chunks(linefreqs, start, end, threshold, level=1):
 
    print '-------------'
    print 'Level '+str(level)
 
    thisgrid = grid() # With no arguments this returns the default wavegrid. Need to preserve that behaviour.
    thisgrid = thisgrid.freq

    #we presume that the chunk is too big, so we split it before checking again
    
    starts = [start, start + (end-start)/2]
    ends = [start + (end-start)/2 -1, end]
    
    #print 'from split chunks, level'+str(level)
    print '====='
    print starts, ends
    
    sizes = [ x-y for x,y in zip(ends,starts) ]#(ends-starts)
    n_points = [np.size(thisgrid[ (thisgrid > linefreqs[x]) & (thisgrid < linefreqs[y]) ]) for x,y in zip(ends,starts)] 

    #find n points in overall grid between line start[i] and end[i]
    print sizes
    print n_points
    
    chunk_sizes = [x*y for x,y in zip(sizes,n_points)]
    #print arr_sizes
        
    need_subdiv = [x > threshold for x in chunk_sizes]
    
    initstarts=starts
    initends=ends
    
    print 'subdiv:'
    print need_subdiv

    if need_subdiv[0] == 1: # & level < 3:
        print 'Branch 1'
        print initstarts[0], initends[0]
        subchunks0 = split_chunks(linefreqs, initstarts[0], initends[0], threshold, level=level+1) 
        print subchunks0[0], starts
        starts = sorted(starts + [i for i in subchunks0[0] if i not in starts])
        
        print subchunks0[1], ends
        ends = sorted(ends + [i for i in subchunks0[1] if i not in ends])
                
    if need_subdiv[1] == 1: # & level < 3:
        print 'Branch 2'
        print initstarts[1], initends[1]
        subchunks0 = split_chunks(linefreqs, initstarts[1], initends[1], threshold, level=level+1) 
        print subchunks0[0], starts
        starts = sorted(starts + [i for i in subchunks0[0] if i not in starts])
        
        print subchunks0[1], ends
        ends = sorted(ends + [i for i in subchunks0[1] if i not in ends])
                    

    
    
    return [starts, ends]




def n_max(arr, n):
    indices = arr.ravel().argsort()[-n:]
    indices = (np.unravel_index(i, arr.shape) for i in indices)
    return [[arr[i], i[0]] for i in indices]