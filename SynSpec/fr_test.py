# test file reading

import numpy as np
from numpy import *
import csv
from sfClasses import *




with open('/home/dstock/sf/py/code/filereadtest/01_01_hit04.par') as csvfile:
    reader = csv.reader(csvfile)
    i=0
    for row in reader:
	print i
        #print('{:2d}{:1d}{:12.6f}{:10.3e}{:10.3e}{:5.4f}{:5.4f}{:10.4f}{:4.2f}{:8.6f}{:15c}{:15c}{:15c}{:15c}{:6d}{:12d}{:1c}{:7.1f}{:7.1f}'.format(row[0]))
        data = row[0]
	thisLine = Line(uint(data[0:2]), uint(data[2:3]), float64(data[3:15]), float64(data[25:35]))
	if i == 0:
		lines = np.array([thisLine])
		print lines.size
	else:	
		lines = np.append(lines, thisLine)
		print lines.size, 'balls'
	i+=1


print type(lines), lines.size


