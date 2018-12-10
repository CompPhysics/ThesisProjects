from numpy import *

def f(x,y,z):
    return x**2 + y**3 + z**3

N = 4
energies = []
for x in range(N):
    for y in range(N):
        for z in range(N):
            energies.append([x,y,z,f(x,y,z)])
b = {}
for i in energies:
    #print i[3]
    b[i[3]] = 0
    
for i in energies:
    #print i[3]
    b[i[3]] += 1

print b
for i in b:
    print i