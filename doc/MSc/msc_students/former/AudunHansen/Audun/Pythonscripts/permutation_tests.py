#Matrix permutation tests
from time import *
from numpy import *

Np = 2
Nh = 2




Z = random.randint(0,10,(Np**2, Nh**2))
A = zeros((Np**2, Nh**2), dtype = int)
Zai = zeros((Np*Nh, Np*Nh), dtype = int)
#permutation
t0 = clock()
for i in range(Nh):
    for j in range(Nh):
        for a in range(Np):
            for b in range(Np):
                A[a + b*Np, i + j * Nh] = Z[b + a*Np, j + i*Nh]
                Zai[a + i*Np, j + b*Nh] = Z[b + a*Np, j + i*Nh]
print "====Zai===="
print Zai
print "==========="








t1 = clock()
print "Time spent on iterative permutation:", t1-t0
B = zeros((Np**2, Nh**2), dtype = int)
#print linspace(0,Np**2, Np)
Tab = zeros((Np**2), dtype = int)
for b in range(Np):
    Tab[range(b*Np,(b+1)*Np)] = range(b,Np**2, Np)

Tai = zeros((Np*Nh), dtype = int)
for a in range(Nh):
    Tai[range(a*Nh,(a+1)*Nh)] = range(b,Np*Nh, Nh)

Tij = zeros((Nh**2), dtype = int)
for i in range(Nh):
    Tij[range(i*Nh,(i+1)*Nh)] = range(i,Nh**2,Nh)
    
#print "Tb:", Tb
b = range(0,Np**2, Np)
d = range(1,Np**2,Np)
c = range(0,Nh**2, Nh)
#print b
#print d
#print Z[0:1:Np**2, :]
#print Z
#print "Compare:"
#print "ZTb:",Z[Tb]
t0 = clock()
Ztb = Z[Tab][:, Tij]  #Permute ab
t1 = clock()
print "Time spent on permutation:", t1-t0
#Zta = Ztb[:, Ti]
print "Zta:", Ztb
print "A:", A
print "==="
Vphhp = random.randint(0,10,(Np*Nh, Np*Nh))
Vhpph = zeros((Np*Nh, Np*Nh), dtype = int)

for i in range(Nh):
    for j in range(Nh):
        for a in range(Np):
            for b in range(Np):
                Vhpph[a + i*Np, j + b*Nh] = Vphhp[0,0]

#permute 