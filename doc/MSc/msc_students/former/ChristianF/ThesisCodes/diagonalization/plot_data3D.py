from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from math import factorial

# Position in first column(s), eigenfunction in the remainder:
eigenvectors = np.loadtxt('PlotAndData/Eigenvectors.dat')

psiX = np.loadtxt('PlotAndData/SeparateEigenvectorsX.dat')
psiY = np.loadtxt('PlotAndData/SeparateEigenvectorsY.dat')

constants = np.loadtxt('PlotAndData/Constants.dat')
positionvectors = np.loadtxt('PlotAndData/Positionvectors.dat')
supX = np.loadtxt('PlotAndData/SeparateSuperpositionsX.dat')
supY = np.loadtxt('PlotAndData/SeparateSuperpositionsY.dat')

supX[:,0] = supX[:,0]/np.sqrt(np.dot(supX[:,0],supX[:,0]))
supY[:,0] = supY[:,0]/np.sqrt(np.dot(supY[:,0],supY[:,0]))

psiX[:,0] = psiX[:,0]/np.sqrt(np.dot(psiX[:,0],psiX[:,0]))
psiY[:,0] = psiY[:,0]/np.sqrt(np.dot(psiY[:,0],psiY[:,0]))


print(np.dot(supX[:,0],supX[:,0]))
print(np.dot(supY[:,0],supY[:,0]))
print(np.dot(supX[:,0]*supY[:,0],supX[:,0]*supY[:,0]))

omega, nDim, Lx, Ly, Lz, N, numEigFunctions, h = constants

numEigFunctions = int(numEigFunctions)
N = int(N)-1
nDim = int(nDim)

# Position vector:
r = np.zeros((nDim, N))
for d in range(nDim):
    r[d] = positionvectors[:,d]

# Numerically calculated eigenvector(2) (for different quantum numbers n):
psi = np.zeros((numEigFunctions,N))
for i in range(numEigFunctions):
    psi[i] = eigenvectors[:,i]

fig = plt.figure()
# ax = fig.gca(projection='3d')
X, Y = np.meshgrid(r[0], r[1])


def plotfunc(n):
    Z = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            Z[i,j] = (psiX[i,n]*psiY[j,n])**2
            
    return Z

zPsi = plotfunc(0)

def plotsup(n):
    Z = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            Z[i,j] = (supX[i,n]*supY[j,n])**2

    return Z

quant = np.zeros((10,2))
quant[0, 0] = 0; quant[0, 1] = 0
quant[1, 0] = 1; quant[1, 1] = 0
quant[2, 0] = 0; quant[2, 1] = 1
quant[3, 0] = 2; quant[3, 1] = 0
quant[4, 0] = 1; quant[4, 1] = 1
quant[5, 0] = 3; quant[5, 1] = 0
quant[6, 0] = 0; quant[6, 1] = 2


def plotASD(n):
	nx = int(quant[n, 0])
	ny = int(quant[n, 1])
	Z = np.zeros((N, N))
	for i in range(N):
		for j in range(N):
			Z[i,j] = (supX[i,nx]*supY[j,ny])**2

	return Z

zSup = plotASD(0)


# surf = ax.plot_surface(X, Y, zSup-zPsi, rstride=1, cstride=1, cmap=cm.jet,
#                        linewidth=0, antialiased=False)
# fig.colorbar(surf, shrink=0.5, aspect=5)

surf = plt.contourf(X,Y, zSup)
#surf = plt.contourf(X,Y, zPsi-zSup)
plt.colorbar(surf)

plt.xlabel("x")
plt.ylabel("y")

plt.show()
