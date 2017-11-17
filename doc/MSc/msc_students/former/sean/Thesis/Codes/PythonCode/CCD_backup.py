from sympy import *
from pylab import *
import numpy as np
import matplotlib.pyplot as plt

#we use natural units
N = 2						#number of particles
NB = 4						#number of closed-shells (n^2=0, n^2=1, n^2=2, etc..., for NB=2, can at max have N=14)
rs = 2						#Wigner Seitz radius
rb = 268					#Bohr radius [MeV^-1]
m = 0.511					#electron mass [MeV]
L = rb*rs*(4*pi/3)**(1./3)	#box length

#box potential
def makeStateSpace(N,NB): #N = num of particles, NB = size of basis (number of energy levels)
	states = []
	for n2 in range(0,NB):	#sum over shells
		for nx in range(-n2-1,n2+1):
			for ny in range(-n2-1,n2+1):
				for nz in range(-n2-1,n2+1):
					if (nx**2 + ny**2 + nz**2 == n2):
						states.append(np.asarray([n2,nx,ny,nz,1]))
						states.append(np.asarray([n2,nx,ny,nz,-1]))
	below_fermi = range(0,N/2)
	above_fermi = range(N/2, NB)
	return np.asarray(states), below_fermi, above_fermi

states, below_fermi, above_fermi = makeStateSpace(N,NB)
NS = len(states)		#number of single-particle states
print NS

def f(p,q):
	returnVal = 0
	if p == q:
		returnVal += states[p][0]*(2*pi)**2/(2*m*L**2)
	
	for i in below_fermi:
		returnVal += assym(p,i,q,i)
		
	return returnVal

def kd3(V1, V2):
	counter = 0
	for i in range(0,3):
		if V1[i] == V2[i]:
			counter += 1
	return (counter == 3)

def assym(p,q,r,s):
    #n21, nx1, ny1, nz1, s1 = states[p]
    #n22, nx2, ny2, nz2, s2 = states[q]
    #n23, nx3, ny3, nz3, s3 = states[r]
    #n24, nx4, ny4, nz4, s4 = states[s]
    
    kp = states[p][1:4]; sp = states[p][4]
    kq = states[q][1:4]; sq = states[q][4]
    kr = states[r][1:4]; sr = states[r][4]
    ks = states[s][1:4]; ss = states[s][4]
    returnVal = 0

    if kd3(kp+kq, kr+ks):
		if sp==sr and sq==ss:
			if kd3(kp,kr) != 1:
				returnVal += 1/(np.inner(kr-kp, kr-kp))
		if sp==ss and sq==sr:
			if kd3(kp,ks) != 1:
				returnVal -= 1/(np.inner(ks-kp, ks-kp))
	
    return returnVal/(L*pi)

    #if kd3([nx1 + nx2, ny1 + ny2, nz1 + nz2], [nx3 + nx4, ny3 + ny4, nz3 + nz4]):
	#	if s1==s3 and s2==s4:
	#		if not kd3([nx1,ny1,nz1],[nx3,ny3,nz3]):
	#			returnVal += 1/((nx3-nx1)**2 + (ny3-ny1)**2 + (nz3-nz1)**2)
	#	if s1==s4 and s2==s3:
	#		if not kd3([nx1,ny1,nz1],[nx4,ny4,nz4]):
	#			returnVal -= 1/((nx4-nx1)**2 + (ny4-ny1)**2 + (nz4-nz1)**2)
	#
    #return returnVal
    
    
def makeFockMatrix():
	fockMat = np.zeros((NS,NS))
	for p in range(0,NS):
		for q in range(0,NS):
			fockMat[p,q] = f(p,q)
	return fockMat
    
def makeIntMatrix():
	intMat = np.zeros((NS,NS,NS,NS))	#array
	for p in range(0,NS):
		for q in range(0,NS):
			for r in range(0,NS):
				for s in range(0,NS):
					intMat[p,q,r,s] = assym(p,q,r,s)
	return intMat

#before I learned I should just store a 4D matrix
#	intMat = np.zeros(NS**4).reshape(NS**2, NS**2)
#	indices = []
#	for p in range(0,NS):
#		for q in range(0,NS):
#			indices.append([p,q])	#this is pointless, might as well use einsum with 4d tensors
#	#print indices
#	counter1 = 0
#	counter2 = 0
#	for h1 in indices:
#		for h2 in indices:
#			p, q = h1
#			r, s = h2
#			intMat[counter1][counter2] = assym(p,q,r,s)
#			counter2 = counter2 + 1
#		counter2 = 0
#		counter1 = counter1 + 1
#	return intMat

V = makeIntMatrix()
F = makeFockMatrix()

#Diagrams
def D1(i,j,a,b):
	factor 		= -1
	returnVal 	= 	 np.einsum('c,c', F[a,N:NS], T[i,j,:,b-N]) \
				   - np.einsum('c,c', F[b,N:NS], T[i,j,:,a-N])
	return factor*returnVal

def D2(i,j,a,b):
	factor 		= 1
	returnVal 	=	 np.einsum('k,k', F[0:N,j], T[i,:,a-N,b-N]) \
				   - np.einsum('k,k', F[0:N,i], T[j,:,a-N,b-N])
	return factor*returnVal
	
def D3(i,j,a,b):
	factor 		= 0.5
	returnVal 	= 	 np.einsum('cd,cd', V[a,b,N:NS,N:NS], T[i,j,:,:])
	return factor*returnVal
	
def D4(i,j,a,b):
	factor 		= 0.5
	returnVal 	= 	 np.einsum('kl,kl', V[0:N,0:N,i,j], T[:,:,a-N,b-N])
	return factor*returnVal
	
def D5(i,j,a,b):
	factor 		= -1
	returnVal 	= 	 np.einsum('kc,kc', V[a,0:N,N:NS,j], T[i,:,:,b-N]) \
				   - np.einsum('kc,kc', V[b,0:N,N:NS,j], T[i,:,:,a-N]) \
				   - np.einsum('kc,kc', V[a,0:N,N:NS,i], T[j,:,:,b-N]) \
				   + np.einsum('kc,kc', V[b,0:N,N:NS,i], T[j,:,:,a-N])
	return factor*returnVal
	
def D6(i,j,a,b):
	factor 		= 0.25
	returnVal 	= 	 np.einsum('cd,klcd,kl', T[i,j,:,:], V[0:N,0:N,N:NS,N:NS], T[:,:,a-N,b-N])
	return factor*returnVal
	
def D7(i,j,a,b):
	factor 		= 0.5
	returnVal 	=  	 np.einsum('kc,klcd,ld', T[i,:,:,a-N], V[0:N,0:N,N:NS,N:NS], T[j,:,:,b-N]) \
				   - np.einsum('kc,klcd,ld', T[j,:,:,a-N], V[0:N,0:N,N:NS,N:NS], T[i,:,:,b-N]) \
				   - np.einsum('kc,klcd,ld', T[i,:,:,b-N], V[0:N,0:N,N:NS,N:NS], T[j,:,:,a-N]) \
				   + np.einsum('kc,klcd,ld', T[j,:,:,b-N], V[0:N,0:N,N:NS,N:NS], T[i,:,:,a-N])
	return factor*returnVal
	
def D8(i,j,a,b):
	factor 		= -0.5
	returnVal 	=  	 np.einsum('klc,klcd,d', T[:,:,:,a-N], V[0:N,0:N,N:NS,N:NS], T[i,j,:,b-N]) \
				   - np.einsum('klc,klcd,d', T[:,:,:,b-N], V[0:N,0:N,N:NS,N:NS], T[i,j,:,a-N])
	return factor*returnVal
	
def D9(i,j,a,b):
	factor 		= -0.5
	returnVal 	=    np.einsum('kcd,klcd,l', T[i,:,:,:], V[0:N,0:N,N:NS,N:NS], T[j,:,a-N,b-N]) \
				   - np.einsum('kcd,klcd,l', T[j,:,:,:], V[0:N,0:N,N:NS,N:NS], T[i,:,a-N,b-N])
	return factor*returnVal
	
def diagrams(i,j,a,b):
	#if abs(D1(i,j,a,b) + D2(i,j,a,b) + D3(i,j,a,b) + D4(i,j,a,b) + D5(i,j,a,b) + D6(i,j,a,b) + D7(i,j,a,b) + D8(i,j,a,b) + D9(i,j,a,b) ) > 1e-1:
	#	print i,j,a,b, D1(i,j,a,b), D2(i,j,a,b), D3(i,j,a,b), D4(i,j,a,b), D5(i,j,a,b), D6(i,j,a,b), D7(i,j,a,b), D8(i,j,a,b), D9(i,j,a,b) 
	return D1(i,j,a,b) + D2(i,j,a,b) + D3(i,j,a,b) + D4(i,j,a,b) + D5(i,j,a,b) + D6(i,j,a,b) + D7(i,j,a,b) + D8(i,j,a,b) + D9(i,j,a,b)



eps = 1e-10
conFac = 1

#set T from MBPT2
T = np.zeros((N, N, NS-N, NS-N))	#hole-indices first, then particle-indices
for i in range(0,N):
	for j in range(0,N):
		for a in range(N,NS):
			for b in range(N,NS):
				if (F[j,j] + F[i,i] - F[a,a] - F[b,b]) > eps:
					T[i,j,a-N,b-N] = V[a,b,i,j]/(F[j,j] + F[i,i] - F[a,a] - F[b,b])
				else:
					T[i,j,a-N,b-N] = 0
					
#print (np.einsum('klcd,klcd',V[0:N,0:N,N:NS,N:NS],T)/(14*27.211383))*1e6

#iterate for T
ECCD = np.einsum('klcd,klcd',V[0:N,0:N,N:NS,N:NS],T)
while conFac > eps:
	T_new = np.zeros((N, N, NS-N, NS-N))
	for i in range(0,N):
		for j in range(0,N):
			for a in range(N,NS):
				for b in range(N,NS):
					if (F[j,j] + F[i,i] - F[a,a] - F[b,b]) != 0:
						T_new[i,j,a-N,b-N] = (V[a,b,i,j] + diagrams(i,j,a,b))/(F[j,j] + F[i,i] - F[a,a] - F[b,b])
					#print T[0,1,0,-1]
					#print T[i,j,a-N,b-N]

	T = T_new
	ECCD_new = np.einsum('klcd,klcd',V[0:N,0:N,N:NS,N:NS],T)
	
	#print T
	#print ECCD_new
	
	conFac = abs(ECCD_new - ECCD)
	ECCD = ECCD_new
	print ECCD
