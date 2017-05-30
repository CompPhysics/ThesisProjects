import numpy as np
import matplotlib.pyplot as plt

x1,y1 = np.loadtxt("build-CCD_prototype-Desktop-Release/example1.txt",unpack=True)
x2,y2 = np.loadtxt("build-CCD_prototype-Desktop-Release/example2.txt",unpack=True)

fig,ax1 = plt.subplots()
ax1.set_title("Plot title...")    
ax1.set_xlabel('your x label..')
ax1.set_ylabel('your y label...')

ax1.plot(x1,y1, c='r', label='numerical')
#ax1.plot(x2,y2, c='b', label='analytical')

leg = ax1.legend()

plt.show()
