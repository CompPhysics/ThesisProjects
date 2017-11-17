import numpy as np
from mayavi import mlab
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import sys

def readData(filename):

    infile = open("%s" %filename, 'r')
			
    #positions = [[],[],[]]
    positions = []
    for i in range(int(sys.argv[1])):
        positions.append([])

    for line in infile:
        data = line.split()
        for i in range(len(data)):
            positions[i].append(float(data[i]))
        #positions[1].append(float(data[1]))
        #positions[2].append(float(data[2]))
    
    infile.close()

    return np.asarray(positions)

def plotData(dataFile):
    positions = readData(dataFile)#"positionsN10e6Int.dat")
    #x = positions[0]
    #y = positions[1]
    #z = positions[2]
    p = positions
    #if len(p) == 2:
    x = p[0]
    y = p[1]
    
    #nx, binsx, patchesx = plt.hist(x, bins=100, normed=1)
    #ny, binsyY = np.histogram(y, bins=100, normed=1)
    #bincenters = 0.5*(binsyY[1:]+binsyY[:-1])
    #ny /= 2
    hist, binsx, binsy = np.histogram2d(x, y, bins=100, normed=1)
    print len(hist)

    #xZeros = np.zeros(len(binsy[:-1]))
    #xThrees = xZeros-3
    X, Y = np.meshgrid(binsx[:-1], binsy[:-1])
    #X2, Y2 = np.meshgrid(xZeros, binsy[:-1])
    #X3, Y3 = np.meshgrid(xThrees, binsy[:-1])

    #plt.plot(bincenters, ny)
    #plt.axvline(rMean, color="r", label="Mean r = %f" %rMean)
    fig = plt.figure()
    
    #plt.hold("on")
    ax = fig.gca(projection='3d')
    #ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.07))
    #ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.07))
    #ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.07))
    
    #ax.plot_wireframe(X2, Y2, hist)
    #ax.plot(bincenters, ny, zs=-3, zdir='x')
    
    ax.plot_surface(X, Y, hist, linewidth=0, antialiased=False, rcount=100, ccount=100, cmap=cm.Blues)
    #ax.plot_surface(X, Y, hist, linewidth=0, antialiased=False, rcount=100, ccount=100)   

    xlim = ax.get_xlim()[0]
    levels = [0.0]
    ax.contour(X, Y, hist, levels, zdir='x', offset=xlim)
    #plt.legend()
    #plt.title("N=10, MC cycles=1e6, Alpha=0.5, D=3")
    #plt.xlim([-1.5, 1.5])
    #plt.ylim([-1.5, 1.5])
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    #ax.set_zlim([0, 0.6])    
    #ax.set_zlabel("Number of Samples (Normalized)")

    #s = mlab.mesh(X, Y, hist)
    #mlab.show()
    
    #return n, bins, rMean

if __name__ == "__main__":
    n = []
    bins = []
    rMean = []
    for i in range(2, len(sys.argv)):
        plotData(sys.argv[i])    
    #    n_i, bins_i, rMean_i = plotData(sys.argv[i])
    #    n.append(n_i)
    #    bins.append(bins_i)
    #    rMean.append(rMean_i)
    #    plt.close()    

    #plt.figure()
    #plt.hold("on")
    #lab = ["Harmonic Oscillator", "No Jastrow", "All interactions"]
    #lab2 = ["HO", "NJ", "Ai"]
    #style = ["", "--", ""]
    #clr = ["b", "g", "r"]
    #for i in range(len(n)):
    #    plt.plot(bins[i][:-1], n[i], style[i], label=lab[i])
    #    lab2[i] += ", Mean r = %f" %rMean[i]
    #    plt.axvline(rMean[i], color= clr[i], label=lab2[i])
    #plt.xlabel("r")
    #plt.ylabel("Number of Samples (Normalized)")
    #plt.title("N=20, MC cycles=1e6")
    #plt.legend()
    plt.show()
