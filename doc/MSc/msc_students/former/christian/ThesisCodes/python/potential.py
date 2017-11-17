import numpy as np
import matplotlib.pylab as plt
import sys
'''
N = 200 

h = 10./N

L = [2.5,0,0]
o = 1.

rho = -5 + np.linspace(0, N, N+1)*h;
V = 0.5*o*o*(np.square(rho)-2*abs(rho)*L[0]+L[0]*L[0]);

for i in range(int(N/2)):
    V[i] += 2;

plt.plot(rho,V)
plt.show()
'''
positionvectors = np.loadtxt('../diagonalization/PlotAndData/Positionvectors.dat')
#potential = np.loadtxt('../diagonalization/PlotAndData/Potential.dat')
constants = np.loadtxt('../diagonalization/PlotAndData/Constants.dat')

potentialSquare = np.loadtxt('../diagonalization/PlotAndData/PotentialSquare.dat')
potentialHarmOmega1 = np.loadtxt('../diagonalization/PlotAndData/PotentialHarmOmega1.dat')
potentialHarmOmega010 = np.loadtxt('../diagonalization/PlotAndData/PotentialHarmOmega010.dat')

omega, nDim, Lx, Ly, Lz, N, numEigFunctions, h = constants
N = int(N)-1

# Position vector:
r = np.zeros((int(nDim), N))
for d in range(int(nDim)):
    r[d] = positionvectors[:,d]

#plt.plot(r[0], potential[1:-1])
plt.plot(r[0], potentialSquare[1:-1,0])
plt.plot(r[0], potentialHarmOmega1[1:-1,0])
plt.plot(r[0], potentialHarmOmega010[1:-1,0])
#plt.title("2D Double Harmonic Oscillator Well")

"""
x1 = [-0.5, 0.5]
x2_1 = [-1.5, -0.5]
x2_2 = [0.5, 1.5]
x3_1 = [-1.5, -2.5]
x3_2 = [-0.5, 0.5]
x3_3 = [1.5, 2.5]

E1 = [0.331966, 0.331966]
E2 = [0.789572, 0.789572]
E3 = [1.212873, 1.212873]

plt.plot(x1, E1)
plt.plot(x2_1, E2, x2_2, E2)
plt.plot(x3_1, E3, x3_2, E3, x3_3, E3)

xP1 = [-0.25]
xP2 = [0.25]
xP3 = [-1.25]
xP4 = [-0.75]
xP5 = [0.75]
xP6 = [1.25]
xP7 = [-2.25]
xP8 = [-1.75]
xP9 = [-0.25]
xP10 = [0.25]
xP11 = [1.75]
xP12 = [2.25]

E1P = [0.331966]
E2P = [0.789572]
E3P = [1.212873]


plt.plot(xP1, E1P, ".k", xP2, E1P, ".k")
plt.plot(xP3, E2P, ".k", xP4, E2P, ".k", xP5, E2P, ".k", xP6, E2P, ".k")
plt.plot(xP7, E3P, ".k", xP8, E3P, ".k", xP9, E3P, ".k", xP10, E3P, ".k", xP11, E3P, ".k", xP12, E3P, ".k")
"""

plt.legend(["Square Well", "HO, Omega = 1.00", "HO, Omega = 0.10"])
#plt.legend(["Square Well", "HO, Omega = 1.00", 
#	    "State 1", "State 2", "State 3", "State 4", "State 5", "State 6",
#	    "Particles"])
plt.xlabel("x")
plt.ylabel("V(x)")
plt.xlim([-10.0, 10.0])
plt.ylim([-0.1, 1.3])
plt.plot(r[0], np.zeros(len(r[0])), "--")



plt.show()

for i in range(N):
    print potential[i,0]




