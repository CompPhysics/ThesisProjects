from sympy import *
from pylab import *
import numpy as np
import matplotlib.pyplot as plt

#CLASSES

class MP:	#Minnesota Potential
	def __init__(self):
		pass
		
	def makeStateSpace(self):
		states = []
		for n2 in range(0,NB):	#sum over shells
			for nx in range(-n2-1,n2+1):
				for ny in range(-n2-1,n2+1):
					for nz in range(-n2-1,n2+1):
						if (nx**2 + ny**2 + nz**2 == n2):
							states.append(np.asarray([n2,nx,ny,nz,1,1]))	#sz =  1, tz =  1
							states.append(np.asarray([n2,nx,ny,nz,1,-1]))	#sz =  1, tz = -1
							states.append(np.asarray([n2,nx,ny,nz,-1,1]))	#sz = -1, tz =  1
							states.append(np.asarray([n2,nx,ny,nz,-1,-1]))	#sz = -1, tz = -1
		NS = len(states)
		below_fermi = range(0,N)
		above_fermi = range(N, NS)
		return np.asarray(states), below_fermi, above_fermi, NS
	
	def f(self,p,q):
		pass
	
	def h0(self,p,q):
		pass
	
	def assym(self,p,q,r,s):
		pass
	

class HEG:	#Homogenous Electron Gas
	def __init__(self):
		pass
	
	#box potential
	def makeStateSpace(self): #N = num of particles, NB = size of basis (number of energy levels)
		states = []
		for n2 in range(0,NB):	#sum over shells (energy levels)
			for nx in range(-n2-1,n2+1):
				for ny in range(-n2-1,n2+1):
					for nz in range(-n2-1,n2+1):
						if (nx**2 + ny**2 + nz**2 == n2):
							states.append(np.asarray([n2,nx,ny,nz,1]))
							states.append(np.asarray([n2,nx,ny,nz,-1]))
		NS = len(states)
		below_fermi = range(0,N)
		above_fermi = range(N, NS)
		return np.asarray(states), below_fermi, above_fermi, NS
	
	def f(self,p,q):
		returnVal = self.h0(p,q)
		for i in range(0,N):
			returnVal += self.assym(p,i,q,i)
			
		return returnVal
	
	def h0(self,p,q):
		if p == q:
			energy = states[p][0]
			return energy*(2*pi)**2/(2*m*L**2)
		else:
			return 0
	
	def assym(self,p,q,r,s):
		kp = states[p][1:4]; sp = states[p][4]
		kq = states[q][1:4]; sq = states[q][4]
		kr = states[r][1:4]; sr = states[r][4]
		ks = states[s][1:4]; ss = states[s][4]
		returnVal = 0
	
		if self.kd3(kp+kq, kr+ks):#these deltas might be length-, rather than element-wise
			if sp==sr and sq==ss:
				if self.kd3(kp,kr) == 0:
					returnVal += 1./(np.linalg.norm(kr-kp))**2
			if sp==ss and sq==sr:
				if self.kd3(kp,ks) == 0:
					returnVal -= 1./(np.linalg.norm(ks-kp))**2
		
		return returnVal/(L*pi)
		
	def kd(self,S1,S2):
		return (S1==S2)
		
	def kd3(self,V1, V2):
		counter = 0
		for i in range(0,3):
			if V1[i] == V2[i]:
				counter += 1
		return (counter == 3)



class PM:	#Pairing Model
	def __init__(self):
			pass
		
	#linear in energy (not box potential)
	def makeStateSpace(self): #N = num of particles, NB = size of basis (levels, each level has two spin states)
		states = []
		for i in range(1,NB/2+1):
			states.append(np.asarray([i,1]))
			states.append(np.asarray([i,-1]))
		NS = len(states)
		below_fermi = range(0,N)
		above_fermi = range(N, NS)
		return np.asarray(states), below_fermi, above_fermi, NS
	
	def h0(self,p,q):
		if p == q:
			p1, s1 = states[p]
			return (p1 - 1)
		else:
			return 0
	
	def f(self,p,q):
		#if p == q:
		#	return 0
		s = self.h0(p,q)
		for i in range(0,N):
			s += self.assym(p,i,q,i)
			#print s,p,q,i
		return s
	
	def assym(self,p,q,r,s):
		p1, s1 = states[p]
		p2, s2 = states[q]
		p3, s3 = states[r]
		p4, s4 = states[s]
	
		if p1 != p2 or p3 != p4:
			return 0
		if s1 == s2 or s3 == s4:
			return 0
		if s1 == s3 and s2 == s4:
			return -g/2.
		if s1 == s4 and s2 == s3:
			return g/2.
	
	def eps(self,holes, particles):
		E = 0
		for h in holes:
			p, s = states[h]
			E += (p-1)
		for p in particles:
			p, s = states[p]
			E -= (p-1)
		return E
		
	def FCI(self):	#IT LIIIIIVEEEEES (defining the few lines below took me half a day...)
		#matSize = 0
		#
		##divide by 2 since we work in spin pairs, and below_fermi/above_fermi are given in single-particle states
		#NFb = len(below_fermi)
		#NFa = len(above_fermi)
		##print NFb, NFa
		#determinants = 1
		#subtractor = 0
		#matSize = 0
		#
		#
		#matSize += determinants	
		#determinants *= NFb*NFa		#this is the number of possible singles, doubles, triples, etc, (pair) excitations
		#subtractor += 1
		#
		#while determinants != 0:
		#	matSize += determinants	
		#	determinants *= (NFb - subtractor)*(NFa - subtractor)/4.		#this is the number of possible singles, doubles, triples, etc, (pair) excitations
		#	determinants = int(determinants)
		#	subtractor += 1
		#	
		#H = np.zeros((matSize,matSize))
		
		g_val = g
		
		#I can't help it if this matrix looks terrible in editors other than geany
		H1 = matrix([[2-g_val 	, -g_val/2.	,  -g_val/2., -g_val/2.	, -g_val/2.	,  0		],
					 [-g_val/2.	,  4-g_val	,  -g_val/2., -g_val/2.	,  0.		, -g_val/2.	],
					 [-g_val/2.	, -g_val/2.	,   6-g_val	,  0		, -g_val/2.	, -g_val/2.	],
					 [-g_val/2.	, -g_val/2.	,   0		,  6-g_val	, -g_val/2.	, -g_val/2.	],
					 [-g_val/2.	,  0		,  -g_val/2., -g_val/2.	,  8-g_val	, -g_val/2.	],
					 [0    		, -g_val/2.	,  -g_val/2., -g_val/2.	, -g_val/2.	,  10-g_val	]])
		
		u1, v1 = linalg.eig(H1)
		return min(u1)

#END OF CLASSES

#we use natural units
N 	= 4							#number of particles
NB 	= 8							#number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 can at max have N=14)
rs 	= 1.							#Wigner Seitz radius
rb 	= 268						#Bohr radius [MeV^-1]
m 	= 0.511						#electron mass [MeV]
L 	= rb*rs*(4*pi/3)**(1./3)	#box length
g 	= 1							#pairing strength


system = PM()
PM_on = true

states, below_fermi, above_fermi, NS = system.makeStateSpace()

#print states

def makeFockMatrix():
	fockMat = np.zeros((NS,NS))
	for p in range(0,NS):
		for q in range(0,NS):
			fockMat[p,q] = system.f(p,q)
	return fockMat
    
def makeIntMatrix():
	intMat = np.zeros((NS,NS,NS,NS))
	for p in range(0,NS):
		for q in range(0,NS):
			for r in range(0,NS):
				for s in range(0,NS):
					intMat[p,q,r,s] = system.assym(p,q,r,s)
	return intMat

V = makeIntMatrix()
F = makeFockMatrix()

#print V[0:N,0:N,0:N,0:N]
#print V[0:N,0:N,N:NS,N:NS]
#print V[N:NS,N:NS,N:NS,N:NS]
###############################################################################################################



#Diagrams (audun)
def D3(i,j,a,b):	#L_a
	factor 		= 0.5
	returnVal 	= 	 np.einsum('cd,cd', V[a,b,N:NS,N:NS], T[i,j,:,:])
	return factor*returnVal
	
def D4(i,j,a,b):	#L_b
	factor 		= 0.5
	returnVal 	= 	 np.einsum('kl,kl', V[0:N,0:N,i,j], T[:,:,a-N,b-N])
	return factor*returnVal
	
def D5(i,j,a,b):	#L_c
	factor 		= -1
	returnVal 	= 	 np.einsum('kc,kc', V[a,0:N,N:NS,j], T[i,:,:,b-N]) \
				   - np.einsum('kc,kc', V[b,0:N,N:NS,j], T[i,:,:,a-N]) \
				   - np.einsum('kc,kc', V[a,0:N,N:NS,i], T[j,:,:,b-N]) \
				   + np.einsum('kc,kc', V[b,0:N,N:NS,i], T[j,:,:,a-N])
	return factor*returnVal
	
def D6(i,j,a,b):	#Q_a
	factor 		= 0.25
	returnVal 	= 	 np.einsum('cd,klcd,kl', T[i,j,:,:], V[0:N,0:N,N:NS,N:NS], T[:,:,a-N,b-N])
	return factor*returnVal
	
def D7(i,j,a,b):	#Q_b
	factor 		= 0.5
	returnVal 	=  	 np.einsum('kc,klcd,ld', T[i,:,:,a-N], V[0:N,0:N,N:NS,N:NS], T[j,:,:,b-N]) \
				   - np.einsum('kc,klcd,ld', T[j,:,:,a-N], V[0:N,0:N,N:NS,N:NS], T[i,:,:,b-N]) \
				   - np.einsum('kc,klcd,ld', T[i,:,:,b-N], V[0:N,0:N,N:NS,N:NS], T[j,:,:,a-N]) \
				   + np.einsum('kc,klcd,ld', T[j,:,:,b-N], V[0:N,0:N,N:NS,N:NS], T[i,:,:,a-N])
	return factor*returnVal
	
def D8(i,j,a,b):	#Q_c
	factor 		= -0.5
	returnVal 	=  	 np.einsum('klc,klcd,d', T[:,:,:,a-N], V[0:N,0:N,N:NS,N:NS], T[i,j,:,b-N]) \
				   - np.einsum('klc,klcd,d', T[:,:,:,b-N], V[0:N,0:N,N:NS,N:NS], T[i,j,:,a-N])
	return factor*returnVal
	
def D9(i,j,a,b):	#Q_d
	factor 		= -0.5
	returnVal 	=    np.einsum('kcd,klcd,l', T[i,:,:,:], V[0:N,0:N,N:NS,N:NS], T[j,:,a-N,b-N]) \
				   - np.einsum('kcd,klcd,l', T[j,:,:,:], V[0:N,0:N,N:NS,N:NS], T[i,:,a-N,b-N])
	return factor*returnVal
	
def diagrams(i,j,a,b):	#sum over all L_i and Q_i diagrams
	return D3(i,j,a,b) + D4(i,j,a,b) + D5(i,j,a,b) + D6(i,j,a,b) + D7(i,j,a,b) + D8(i,j,a,b) + D9(i,j,a,b)



###############################################################################################################


#Diagrams (bok)	
def D3b(i,j,a,b):	#L_a
	factor 		= 0.5
	returnVal 	= 	 np.einsum('cd,cd', V[a,b,N:NS,N:NS], T[i,j,:,:])
	return factor*returnVal
	
def D4b(i,j,a,b):	#L_b
	factor 		= 0.5
	returnVal 	= 	 np.einsum('kl,kl', V[0:N,0:N,i,j], T[:,:,a-N,b-N])
	return factor*returnVal
	
def D5b(i,j,a,b):	#L_c
	factor 		= 1
	returnVal 	= 	 np.einsum('kc,kc', V[0:N,b,N:NS,j], T[i,:,a-N,:]) \
				   - np.einsum('kc,kc', V[0:N,b,N:NS,i], T[j,:,a-N,:]) \
				   - np.einsum('kc,kc', V[0:N,a,N:NS,j], T[i,:,b-N,:]) \
				   + np.einsum('kc,kc', V[0:N,a,N:NS,i], T[j,:,b-N,:])
	return factor*returnVal
	
def D6b(i,j,a,b):	#Q_a
	factor 		= 0.25
	returnVal 	= 	 np.einsum('cd,klcd,kl', T[i,j,:,:], V[0:N,0:N,N:NS,N:NS], T[:,:,a-N,b-N])
	return factor*returnVal
	
def D7b(i,j,a,b):	#Q_b
	factor 		= 1
	returnVal 	=  	 np.einsum('kc,klcd,ld', T[i,:,a-N,:], V[0:N,0:N,N:NS,N:NS], T[j,:,b-N,:]) \
				   - np.einsum('kc,klcd,ld', T[j,:,a-N,:], V[0:N,0:N,N:NS,N:NS], T[i,:,b-N,:])
	return factor*returnVal
	
def D8b(i,j,a,b):	#Q_c
	factor 		= -0.5
	returnVal 	=  	 np.einsum('kdc,klcd,l', T[i,:,:,:], V[0:N,0:N,N:NS,N:NS], T[:,j,a-N,b-N]) \
				   - np.einsum('kdc,klcd,l', T[j,:,:,:], V[0:N,0:N,N:NS,N:NS], T[:,i,a-N,b-N])
	return factor*returnVal
	
def D9b(i,j,a,b):	#Q_d
	factor 		= -0.5
	returnVal 	=    np.einsum('lkc,klcd,d', T[:,:,a-N,:], V[0:N,0:N,N:NS,N:NS], T[i,j,:,b-N]) \
				   - np.einsum('lkc,klcd,d', T[:,:,b-N,:], V[0:N,0:N,N:NS,N:NS], T[i,j,:,a-N])
	return factor*returnVal
	
def diagramsb(i,j,a,b):	#sum over all L_i and Q_i diagrams
	return D3b(i,j,a,b) + D4b(i,j,a,b) + D5b(i,j,a,b) + D6b(i,j,a,b) + D7b(i,j,a,b) + D8b(i,j,a,b) + D9b(i,j,a,b)



###############################################################################################################


eps = 1e-10
conFac = 1

#set T from MBPT2
T = np.zeros((N, N, NS-N, NS-N))	#hole-indices first, then particle-indices
for i in range(0,N):
	for j in range(0,N):
		for a in range(N,NS):
			for b in range(N,NS):
				if abs(F[j,j] + F[i,i] - F[a,a] - F[b,b]) > eps:
					T[i,j,a-N,b-N] = V[a,b,i,j]/(F[j,j] + F[i,i] - F[a,a] - F[b,b])
					
#print (np.einsum('klcd,klcd',V[0:N,0:N,N:NS,N:NS],T)/(14*27.211383))*1e6

#iterate for T
ECCD = 0.25*np.einsum('klcd,klcd',V[0:N,0:N,N:NS,N:NS],T)
while conFac > eps:
	T_new = np.zeros((N, N, NS-N, NS-N))
	for i in range(0,N):
		for j in range(0,N):
			for a in range(N,NS):
				for b in range(N,NS):
					if abs(F[j,j] + F[i,i] - F[a,a] - F[b,b]) > eps:
						T_new[i,j,a-N,b-N] = (V[a,b,i,j] + diagrams(i,j,a,b))/(F[j,j] + F[i,i] - F[a,a] - F[b,b])
					#print T[0,1,0,-1]
					#print T[i,j,a-N,b-N]

	T = T_new
	ECCD_new = 0.25*np.einsum('klcd,klcd',V[0:N,0:N,N:NS,N:NS],T)
	
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
print "CCD correlation:", ECCD

