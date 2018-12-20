from numpy import *
from numpy.random import *
from time import *

N = 4000
Z = uniform(0,1,(N,N))

r200 = range(N)
val = zeros((N,N))
t0 = clock()
for i in r200:
    for e in r200:
        val[i,e] += Z[e,i]*Z[i,e]
print clock()-t0     
#print val

t0 = clock()
ZZ = Z.T*Z
print clock()-t0
#print Z
#print " "
#print ZZ
#print " "
#print val


#N = 100
#Z = zeros((N,N,N,N))
#t0 = clock()
##ZZ = tensordot(Z,Z.T, axes = ([1,0],[0,1]))
#print ZZ
#print clock()-t0
