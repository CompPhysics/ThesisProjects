import numpy as np
import matplotlib.pyplot as plt
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
    x = positions[0]
    y = positions[1]
    #z = positions[2]
    p = positions
    r = 0.
    
    for i in range(len(p)):
        r += p[i]**2
    
    r = np.sqrt(r)
    #r = x
    #r = y
    #r = np.sqrt(x**2 + y**2 + z**2)
    #r = np.sqrt(r**2)
    rMean = sum(r)/len(r)
    
    #plt.figure()
    n, bins, patches = plt.hist(r, bins=100, normed=1)
    #plt.axvline(rMean, color="r", label="Mean r = %f" %rMean)
    #plt.legend()
    #plt.title("N=10, MC cycles=1e6, Alpha=0.5, D=3")
    #plt.xlabel(r"$r/a_{ho}$")
    #plt.ylabel("Number of Samples (Normalized)")
    
    return n, bins, rMean
    
if __name__ == "__main__":
    n = []
    bins = []
    rMean = []
    for i in range(2, len(sys.argv)):
        n_i, bins_i, rMean_i = plotData(sys.argv[i])
        n.append(n_i)
        bins.append(bins_i)
        rMean.append(rMean_i)
        plt.close()    

    plt.figure()
    #plt.hold("on")
    lab = ["Double Harmonic Oscillator"]
    lab2 = ["DHO"]
    #lab = ["Harmonic Oscillator", "No Jastrow", "All interactions"]
    #lab2 = ["HO", "NJ", "Ai"]
    style = ["", "--", ""]
    clr = ["b", "g", "r"]
    for i in range(len(n)):
        plt.plot(bins[i][:-1], n[i], style[i], label=lab[i])
        lab2[i] += ", Mean r = %f" %rMean[i]
        plt.axvline(rMean[i], color= clr[i], label=lab2[i])
    plt.xlabel("r")
    plt.ylabel("Number of Samples (Normalized)")
    #plt.title("N=20, MC cycles=1e6")
    plt.legend()
    plt.show()





        
        
