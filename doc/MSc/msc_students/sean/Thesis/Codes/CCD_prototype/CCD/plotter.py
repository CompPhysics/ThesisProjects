#this is a simple plotter where you manually insert data sets (in case you are interested in spesific cases)

import itertools
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
import matplotlib.ticker as ticker

#xData = range(3,40)
#yData = [14.474197721372551,
#		 14.218743793305565,
#		 14.207942941796345,
#		 13.936112367402592,
#		 13.611507926048816,
#		 13.611507926048816,
#		 13.57451376342269,
#		 13.534148280969257,
#		 13.505725917938069,
#		 13.494319125497038,
#		 13.493846612339542,
#		 13.492441495335218,
#		 13.49155729440426,
#		 13.49155729440426,
#		 13.49146272954267,
#		 13.491331913824345,
#		 13.491295036250108,
#		 13.491289243495189,
#		 13.491285320314288,
#		 13.491277400922622,
#		 13.49127102905164,
#		 13.49127102905164,
#		 13.491267741663021,
#		 13.491265278557359,
#		 13.49126066855264,
#		 13.491259278823192,
#		 13.491259278823192,
#		 13.491257416080982,
#		 13.491256539828935,
#		 13.491256539828935,
#		 13.491256450040053,
#		 13.491256192889837,
#		 13.491256035877987,
#		 13.491255921457979,
#		 13.491255880376345,
#		 13.491255847030002,
#		 13.491255798227725,
#		 13.491255798227725,
#		 13.491255789698593,
#		 13.491255774417763]

def divider(inVar):
	return 4./(14**3*(inVar-14)**3)

xData = [54,
		66,
		114,
		162,
		186,
		246,
		294,
		342,
		358,
		406,
		502,
		514,
		610,
		682,
		730,
		778,
		874,
		922,
		970,
		1030,
		1174,
		1238,
		1382,
		1478,
		1502,
		1598,
		1698,
		1790,
		1850,
		1898,
		2042,
		2090]

yData = [11896*divider(54),
		22228*divider(66),
		106956*divider(114),
		277428*divider(162),
		382044*divider(186),
		734312*divider(246),
		1120520*divider(294),
		1590168*divider(342),
		1765168*divider(358),
		2331528*divider(406),
		3747768*divider(502),
		3944292*divider(514),
		5670148*divider(610),
		7221652*divider(682),
		8374140*divider(730),
		9595476*divider(778),
		12321812*divider(874),
		13833692*divider(922),
		15351788*divider(970),
		17367064*divider(1030),
		22834085*divider(1174),
		25559580*divider(1238),
		32157684*divider(1382),
		37077532*divider(1478),
		38328278*divider(1502),
		43501377*divider(1598),
		49047262*divider(1698),
		54968995*divider(1790),
		58847747*divider(1850),
		62067319*divider(1898),
		72170466*divider(2042),
		75635946*divider(2090)]
		
print len(yData)

#for i in range(len(xData)):
#	xData[i] *= 100

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
scale_y = 1e-2
ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_y))
ax.yaxis.set_major_formatter(ticks_y)
#ax.set_yscale('log', basey=10)

#iterate over all data sets
#counter = 0
#for ECC in Vals:
#	plt.plot(Rho, ECC, label=r"$N_h = {}$".format(Nh[counter]), marker=marker.next(), linestyle = linestyle.next())
#	plt.hold('on')
#	counter += 1

plt.plot(xData, yData, marker='s', linestyle = '-')


matplotlib.rcParams.update({'font.size': 18})


plt.legend(loc=4)

plt.xlabel("Single particle states", fontsize=18)
plt.ylabel(r"Density (\%)", fontsize=18)

#plt.xlim(0, 42)

plt.savefig('T3_density.pdf', format='pdf')

plt.tight_layout()
plt.show()
