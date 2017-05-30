#this function allows you to make and store MP as function of Nb

import subprocess
import csv
import itertools
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc

#spesify settings to run
Nb_min = 3 	#starts with this
Nb_max = 41	#ends with this
Nb = range(Nb_min,Nb_max+1)
Nb_repeats = [7,15,23,28,31,39] #last of pairs of Nb that give equal number of states (avoids repeated calculations)

FileName = 'MP_CCD_Nh14_Nb{}_{}.csv'.format(Nb_min, Nb_max)

#run command
cmd = ['mpirun', '-n', '1', './CCD', 'MP', '14', "", '0.2', '1e-15', '0', '4']

if os.path.exists(FileName):
	print "Happy day! That data set is already present in this folder! Exiting program."
	sys.exit(0)

#generate data
for nb in Nb:
	if nb not in Nb_repeats:
		temp = []
		cmd[6] = '{}'.format(nb)
		print cmd
		proc = subprocess.check_output(cmd).split("\n")
		temp.append(float(proc[-5].split()[2]))
		temp.append(float(proc[-4].split()[1]))
		temp.append(float(proc[-3].split()[1]))
		print temp
		
		#saves data to output.csv, "w" will write over previous file, "a" will append to file
		with open(FileName, 'a') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerow(temp)
