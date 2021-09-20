#this function allows you to make, store, and plot MP as function of Nh and rho

import subprocess
import csv
import itertools
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc

#spesify settings to run
Nh = [14,54,186,358]
Rho = np.linspace(0.025, 0.35, 14)
Vals = []

#run command
cmd = ['mpirun', '-n', '1', './CCD', 'MP', "", '20', "", '1e-15', '0', '4']

#generate data
for nh in Nh:
	temp = []
	for rho in Rho:
		#print cmd.format(nh,rho)
		cmd[5] = '{}'.format(nh)
		cmd[7] = '{}'.format(rho)
		print cmd
		proc = subprocess.check_output(cmd).split("\n")
		temp.append(float(proc[-3].split()[1]))
	Vals.append(temp)

#saves data to output.csv
with open('output.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(Vals)

#different marker for each graph
marker = itertools.cycle(('s', 'v', 'o', 'h'))
linestyle = itertools.cycle(('-', '--', '-.', ':'))

#plot commands
fig, ax = plt.subplots()

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax.grid(True)

#iterate over all data sets
counter = 0
for ECC in Vals:
	plt.plot(Rho, ECC, label=r"$N_h = {}$".format(Nh[counter]), marker=marker.next(), linestyle = linestyle.next())
	plt.hold('on')
	counter += 1

plt.legend(loc=4)

plt.xlabel(r"$\rho$ [fm$^{-3}$]", fontsize=16)
plt.ylabel(r"$E_{CCD}/A$ [MeV]", fontsize=16)

plt.xlim(0, 0.375)

plt.savefig('density_for_MP.pdf', format='pdf')

plt.show()
plt.tight_layout()
