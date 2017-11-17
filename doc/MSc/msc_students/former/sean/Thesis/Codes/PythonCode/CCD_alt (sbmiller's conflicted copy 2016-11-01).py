from sympy import *
from pylab import *
import numpy as np
import math
import matplotlib.pyplot as plt

#CLASSES

class HEG:	#Homogenous Electron Gas
	def __init__(self):
		pass
	
	#3d box potential
	def makeStateSpace(self): #N = num of particles, NB = size of basis (number of energy levels)
		states = []
		for n2 in range(0,NB+1):	#sum over shells (energy levels)
			for nx in range(-n2,n2+1):
				for ny in range(-n2,n2+1):
					for nz in range(-n2,n2+1):
						if (nx**2 + ny**2 + nz**2 == n2):
							states.append(np.asarray([n2,nx,ny,nz,1]))
							states.append(np.asarray([n2,nx,ny,nz,-1]))
		NS = len(states)
		below_fermi = range(0, N)
		above_fermi = range(N, NS)
		return np.asarray(states), below_fermi, above_fermi, NS
	
	def f(self,p):
		returnVal = self.h0(p)
		for i in range(0,N):
			returnVal += self.assym(p,i,p,i)
		return returnVal
	
	def h0(self,p):
		energy = states[p][0]
		return energy*2*pi**2/(m*L2)
	
	def assym(self,p,q,r,s):
		kp = states[p][1:4]; sp = states[p][4]
		kq = states[q][1:4]; sq = states[q][4]
		kr = states[r][1:4]; sr = states[r][4]
		ks = states[s][1:4]; ss = states[s][4]
		
		if np.prod(kp+kq==kr+ks):
			returnVal = 0
			if not np.prod(kp==kr):
				returnVal += (sp==sr)*(sq==ss)/(np.linalg.norm(kr-kp))**2
			if not np.prod(kp==ks):
				returnVal -= (sp==ss)*(sq==sr)/(np.linalg.norm(ks-kp))**2
			return returnVal/(L1*pi)
		else:
			return 0
	
	def kUnique(self,i):#i is index i in states, not momentum itself
		k = states[i][1:-1]
		dk = 2*max(k) + 1
		kuni = k[0] + k[1]*dk + k[2]*dk**2
		return kuni
		
	def kUnique2(self,k,p):
		kuni = self.kUnique(k) + self.kUnique(p)
		return kuni
		
	def restrictions(self):
		pass
						

#END OF CLASSES

#we use natural units
N 	= 2							#number of particles
NB 	= 4							#number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 can at max have N=14)
rs 	= 1.0						#Wigner Seitz radius
rb 	= 1.#2.68#268*1e-6/27.211			#Bohr radius [MeV^-1]
m 	= 1.#0.511#*1e6/21.211			#electron mass [MeV]
L3 	= 4*pi*N*rs**3/3.				#box length
L2 	= L3**(2./3)
L1 	= L3**(1./3)
g 	= 1.							#pairing strength

system = HEG()

states, below_fermi, above_fermi, NS = system.makeStateSpace()

#print states

def makeFockVector():	#due to HF basis, this is diagonal, therefore we only keep a 1D array
	fockVec = np.zeros(NS)
	for p in range(0,NS):
		fockVec[p] = system.f(p)
	return fockVec
    

    
def makeIntMatrix():	#just fyi, this function took me half a day, I'm not fluent in numpy functions
	blockArrays = np.zeros(((NS-N)**2,3))						#i really need to start using Np = NS-N
	
	#only need two indices due to matrix diagonal symmetry
	index = 0
	for a in range(N,NS):
		for b in range(N,NS):
			blockArrays[index] = [system.kUnique2(a,b),a,b] 	#[k_u, a, b], i can probably use a+b*Np instead of "index", but whatevs
			index += 1
	
	#blockArrays is now sorted in ascending orders of k_u
	blockArrays = blockArrays[blockArrays[:,0].argsort()]		#we got this from the internets
	
	
	#Okay, everythign below this point is just speculative and to remember my trail of thought
	blockArraysNew = np.zeros((np.shape(blockArrays[:][0]),2))
	#BEWARE: for-loop notation might be wrong here
	for val in np.unique(blockArrays[1:][0]):					#val now takes all different values of k_u
		index = np.where(blockArrays[:][0] == val)				#note, you gt a tuple from np.where
		
		
	return np.asarray(intMat_pppp)
	
print NS
Vpppp = makeIntMatrix()
print NS
F = makeFockVector()

#Diagrams
def D3():	#L_a
	factor 		= 0.5
	returnVal 	= np.einsum('abcd,ijcd->ijab', Vpppp, T)
	return factor*returnVal
	
def diagrams():	#sum over all L_i and Q_i diagrams
	returnVal = D3()
	return returnVal

#########################################################################

eps = 1e-10
conFac = 1

#set T with diagrams giving zero contribution
T = np.zeros((N, N, NS-N, NS-N))	#hole-indices first, then particle-indices
for i in range(0,N):
	for j in range(0,N):
		for a in range(N,NS):
			for b in range(N,NS):
				if abs(F[j] + F[i] - F[a] - F[b]) > eps:
					T[i,j,a-N,b-N] = Vhhpp[i,j,a-N,b-N]/(F[j] + F[i] - F[a] - F[b])

#iterate for T
ECCD = 0.25*np.einsum('klcd,klcd',Vhhpp,T)

denMat = np.zeros((N,N,NS-N,NS-N))

for i in range(0,N):
	for j in range(0,N):
		for a in range(N,NS):
			for b in range(N,NS):
				if abs(F[j] + F[i] - F[a] - F[b]) > eps:
					denMat[i,j,a-N,b-N] = 1./(F[i] + F[j] - F[a] - F[b])

while conFac > eps:
	T_new = np.zeros((N, N, NS-N, NS-N))
	D = diagrams()
	T_new = (Vhhpp+D)*denMat
	T = T_new
	ECCD_new = 0.25*np.einsum('klcd,klcd',Vhhpp,T)
	
	#print T
	#print ECCD_new
	
	conFac = abs(ECCD_new - ECCD)
	ECCD = ECCD_new
	print ECCD

if PM_on == true:	
	E_REF = -g*N/4
	for p in below_fermi:
		E_REF += states[p,0]-1
	
	print "g value:", g
	print "E_REF:", E_REF
	print "FCI correlation:", system.FCI() - E_REF
	print ECCD - system.FCI() + E_REF
print "CCD correlation:", ECCD/N, "Hatrees"
