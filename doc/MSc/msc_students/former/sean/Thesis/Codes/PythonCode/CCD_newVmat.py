from sympy import *
from pylab import *
import numpy as np
import math
import matplotlib.pyplot as plt

#CLASSES

#class MP:	#Minnesota Potential
#	def __init__(self):
#		pass
#		
#	def makeStateSpace(self):
#		states = []
#		for n2 in range(0,NB):	#sum over shells
#			for nx in range(-n2,n2):
#				for ny in range(-n2,n2):
#					for nz in range(-n2,n2):
#						if (nx**2 + ny**2 + nz**2 == n2):
#							states.append(np.asarray([n2,nx,ny,nz,1,1]))	#sz =  1, tz =  1
#							states.append(np.asarray([n2,nx,ny,nz,1,-1]))	#sz =  1, tz = -1
#							states.append(np.asarray([n2,nx,ny,nz,-1,1]))	#sz = -1, tz =  1
#							states.append(np.asarray([n2,nx,ny,nz,-1,-1]))	#sz = -1, tz = -1
#		NS = len(states)
#		below_fermi = range(0,N)
#		above_fermi = range(N, NS)
#		return np.asarray(states), below_fermi, above_fermi, NS
#	
#	def f(self,p,q):
#		returnVal = self.h0(p,q)
#		for i in range(0,N):
#			returnVal += self.assym(p,i,q,i)
#		return returnVal
#	
#	def h0(self,p,q):
#		if p == q:
#			energy = states[p][0]
#			return energy*(2*pi)**2/(2*m*L2)
#		else:
#			return 0
#	
#	def assym(self,p,q,r,s):
#		#kp = states[p][1:4]; sp = states[p][4]; tp = states[p][5]
#		#kq = states[q][1:4]; sq = states[q][4]; tq = states[q][5]
#		#kr = states[r][1:4]; sr = states[r][4]; tr = states[r][5]
#		#ks = states[s][1:4]; ss = states[s][4]; ts = states[s][5]
#		#returnVal = 0
#		#
#		#if self.kd3(kp+kq, kr+ks):
#		#	q = 
#		#	returnVal = V_a*(pi/alpga)**(3./2)*math.exp()
#		pass
#	
#	def kd3(self,V1, V2):
#		counter = 0
#		for i in range(0,3):
#			if V1[i] == V2[i]:
#				counter += 1
#		return (counter == 3)

class HEG:	#Homogenous Electron Gas
	def __init__(self):
		pass
	
	#box potential
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
	
	def f(self,p,q):
		returnVal = self.h0(p,q)
		for i in range(0,N):
			returnVal += self.assym(p,i,q,i)
		#	if p==q:
		#		print self.assym(p,i,q,i), p,i
		#if p==q:
		#	print returnVal,p,i
		return returnVal
	
	def h0(self,p,q):
		if p == q:
			energy = states[p][0]
			return energy*2*pi**2/(m*L2)
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
				if self.kd3(kp,kr) != 1.0:
					returnVal += 1./(np.linalg.norm(kr-kp))**2
			if sp==ss and sq==sr:
				if self.kd3(kp,ks) != 1.0:
					returnVal -= 1./(np.linalg.norm(ks-kp))**2
		return returnVal/(L1*pi)
		
	def kd(self,S1,S2):
		return (S1==S2)
		
	def kd3(self,V1, V2):
		d = 1.0
		for i in range(len(V1)):
			d*=(V1[i]==V2[i])
		return d



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
		returnVal = self.h0(p,q)
		for i in range(0,N):
			returnVal += self.assym(p,i,q,i)
		return returnVal
	
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
N 	= 14							#number of particles
NB 	= 2							#number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 can at max have N=14)
rs 	= 1.6						#Wigner Seitz radius
rb 	= 1.#2.68#268*1e-6/27.211			#Bohr radius [MeV^-1]
m 	= 1.#0.511#*1e6/21.211			#electron mass [MeV]
L3 	= 4*pi*N*rs**3/3.				#box length
L2 	= L3**(2./3)
L1 	= L3**(1./3)
g 	= 1.							#pairing strength

system = HEG()
PM_on = false

states, below_fermi, above_fermi, NS = system.makeStateSpace()

#print states

def makeFockMatrix():
	fockMat = np.zeros((NS,NS))
	for p in range(0,NS):
		fockMat[p,p] = system.f(p,p)
	return fockMat
    
def makeIntMatrix():
	intMat_hhhh = np.zeros((N,N,N,N))			
	intMat_hhpp = np.zeros((N,N,NS-N,NS-N))		
	intMat_hphp = np.zeros((N,NS-N,N,NS-N))
	intMat_pppp = np.zeros((NS-N,NS-N,NS-N,NS-N))
	
	for i in range(0,N):
		for j in range(0,N):
			
			for k in range(0,N):
				for l in range(0,N):
					intMat_hhhh[i,j,k,l] = system.assym(i,j,k,l)
			
			for a in range(N,NS):
				for b in range(N,NS):
					intMat_hhpp[i,j,a-N,b-N] = system.assym(i,j,a,b)
					intMat_hphp[i,a-N,j,b-N] = system.assym(i,a,j,b)	
	
	for a in range(N,NS):
		for b in range(N,NS):			
			for c in range(N,NS):
				for d in range(N,NS):
					intMat_pppp[a-N,b-N,c-N,d-N] = system.assym(a,b,c,d)
	
	return intMat_hhhh, intMat_hhpp, intMat_hphp, intMat_pppp
print NS
Vhhhh, Vhhpp, Vhphp, Vpppp = makeIntMatrix()
F = makeFockMatrix()
#print Vhhhh
#print Vhhpp
#print Vpppp

###############################################################################################################

#Diagrams (audun)
def D3():	#L_a
	factor 		= 0.5
	returnVal 	= np.einsum('cdab,ijcd->ijab', Vpppp, T)
	return factor*returnVal
	
def D4():	#L_b
	factor 		= 0.5
	returnVal 	= np.einsum('klij,klab->ijab', Vhhhh, T)
	return factor*returnVal
	
def D5():	#L_c
	factor 		= -1
	Val 		=  np.einsum('kajc,ikcb->ijab', Vhphp, T)
	returnVal 	=  Val - np.einsum('ijab->jiab',Val) \
					   - np.einsum('ijab->ijba',Val) \
					   + np.einsum('ijab->jiba',Val)
	return factor*returnVal
	
def D6():	#Q_a
	factor 		= 0.25
	returnVal 	= np.einsum('ijcd,klcd,klab->ijab', T, Vhhpp, T)
	return factor*returnVal
	
def D7():	#Q_b
	factor 		= 0.5
	Val 		= np.einsum('ikca,klcd,jldb->ijab', T, Vhhpp, T)
	returnVal 	= Val - np.einsum('ijab->jiab', Val) \
					  - np.einsum('ijab->ijba', Val) \
					  + np.einsum('ijab->jiba', Val)
	return factor*returnVal
	
def D8():	#Q_c
	factor 		= -0.5
	Val 		=  np.einsum('klca,klcd,ijdb->ijab', T, Vhhpp, T)
	returnVal 	=  Val - np.einsum('ijab->ijba', Val)
	return factor*returnVal
	
def D9():	#Q_d
	factor 		= -0.5
	Val 		=  np.einsum('ikcd,klcd,jlab->ijab', T, Vhhpp, T)
	returnVal 	=  Val - np.einsum('ijab->jiab', Val)
	return factor*returnVal
	
def diagrams():	#sum over all L_i and Q_i diagrams
	#print sum(D3()), sum(D4()), sum(D5()), sum(D6()), sum(D7()), sum(D8()), sum(D9())
	returnVal = D3() + D4() + D5() + D6() + D7() + D8() + D9()
	return returnVal



###############################################################################################################


#Diagrams (bok)	


###############################################################################################################


eps = 1e-10
conFac = 1

#set T with diagrams giving zero contribution
T = np.zeros((N, N, NS-N, NS-N))	#hole-indices first, then particle-indices
for i in range(0,N):
	for j in range(0,N):
		for a in range(N,NS):
			for b in range(N,NS):
				if abs(F[j,j] + F[i,i] - F[a,a] - F[b,b]) > eps:
					T[i,j,a-N,b-N] = Vhhpp[i,j,a-N,b-N]/(F[j,j] + F[i,i] - F[a,a] - F[b,b])
					
print (np.einsum('klcd,klcd',Vhhpp,T)) #/(N*27.211383))*1e6

print "Python code"
#iterate for T
ECCD = 0.25*np.einsum('klcd,klcd',Vhhpp,T)
while conFac > eps:
	T_new = np.zeros((N, N, NS-N, NS-N))
	D = diagrams()
	for i in range(0,N):
		for j in range(0,N):
			for a in range(N,NS):
				for b in range(N,NS):
					if abs(F[j,j] + F[i,i] - F[a,a] - F[b,b]) > eps:
						T_new[i,j,a-N,b-N] = (Vhhpp[i,j,a-N,b-N] + D[i,j,a-N,b-N])/(F[j,j] + F[i,i] - F[a,a] - F[b,b])
						#T_new[i,j,a-N,b-N] = (Vhhpp[i,j,a-N,b-N])/(F[j,j] + F[i,i] - F[a,a] - F[b,b])

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

