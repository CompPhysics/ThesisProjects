from sympy import *
#from pylab import *
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
	
	#def assym(self,p,q,r,s):
	#	kp = states[p][1:4]; sp = states[p][4]
	#	kq = states[q][1:4]; sq = states[q][4]
	#	kr = states[r][1:4]; sr = states[r][4]
	#	ks = states[s][1:4]; ss = states[s][4]
	#	
	#	if np.prod(kp+kq==kr+ks):
	#		returnVal = 0
	#		if not np.prod(kp==kr):
	#			returnVal += (sp==sr)*(sq==ss)/(np.linalg.norm(kr-kp))**2
	#		if not np.prod(kp==ks):
	#			returnVal -= (sp==ss)*(sq==sr)/(np.linalg.norm(ks-kp))**2
	#		return returnVal/(L1*pi)
	#	else:
	#		return 0
	
	def assym(self, kp,kq,kr,ks):	#CURRENTLY, THIS RETURNS A 4D MATRIX --> NOT GOOD
		p  = states[kp,1:4][:,None,None,None]
		pm =   states[kp,4][:,None,None,None]
		q  = states[kq,1:4][None,:,None,None]
		qm =   states[kq,4][None,:,None,None]
		r  = states[kr,1:4][None,None,:,None]
		rm =   states[kr,4][None,None,:,None]
		s  = states[ks,1:4][None,None,None,:]
		sm =   states[ks,4][None,None,None,:]
		kdplus = 1.0*np.prod(p+q==r+s,4)
		kspin1 = (pm==rm)*(qm==sm)*(1.-np.prod((p==r),4))/np.sum((r-p)**2,4)
		kspin1[kspin1!=kspin1] = 0
		kspin2 = (pm==sm)*(qm==rm)*(1.-np.prod((p==s),4))/np.sum((s-p)**2,4)
		kspin2[kspin2!=kspin2] = 0
		
		ret = kdplus*(kspin1-kspin2)/float(L1*pi)
		
		#ret is 4d, t2 is the 2d version of it: t2[i,j] = ret[i,i,j,j]
		rows = np.arange(len(kp))
		cols = np.arange(len(kr))
		t1 = ret[rows,rows,:,:]
		t2 = t1[:,cols,cols]
		
		return t2
	
	def kUnique(self,i):#i is index i in states, not momentum itself
		k = states[i][1:-1]
		dk = 2*NB + 1
		kuni = k[0] + k[1]*dk + k[2]*dk**2
		return kuni
		
	def kUnique2(self,k,p):
		#kuni = self.kUnique(k) + self.kUnique(p)
		km = states[k][1:]
		pm = states[p][1:]
		mom = km+pm
		dk = 2*max(mom) + 1
		kuni = mom[0] + mom[1]*dk + mom[2]*dk**2 + mom[3]*dk**3
		#print kuni
		return kuni
		
	def restrictions(self):
		pass

#END OF CLASSES

#we use natural units
N 	= 14							#number of particles
NB 	= 3							#number of closed-shells (n^2=0, n^2=1, n^2=2, etc... For NB=2 can at max have N=14)
rs 	= 1						#Wigner Seitz radius
rb 	= 1.#2.68#268*1e-6/27.211			#Bohr radius [MeV^-1]
m 	= 1.#0.511#*1e6/21.211			#electron mass [MeV]
L3 	= 4*pi*N*rs**3/3.				#box length
L2 	= L3**(2./3)
L1 	= L3**(1./3)
g 	= 1.							#pairing strength
eps = 1e-14
conFac = 1
print L3
system = HEG()

states, below_fermi, above_fermi, NS = system.makeStateSpace()
    
#merely a shell until i vectorize my interaction
#golden rule, don't optimize before you know it works and all that
def temp1(inMat):
	dim = np.shape(inMat)[0]
	vMat = np.zeros((dim,dim))
	kp = np.asarray(inMat[:,0])
	kq = np.asarray(inMat[:,1])
	kr = np.asarray(inMat[:,0])
	ks = np.asarray(inMat[:,1])
	vMat = system.assym(kp,kq,kr,ks)
	#print vMat
	return vMat
    
def temp2(inMat1,inMat2):	#i suggest inMat1 contains rows while inMat2 contains columns, convention is VERY important
	dim1 = np.shape(inMat1)[0]
	dim2 = np.shape(inMat2)[0]
	vMat = np.zeros((dim1,dim2))
	kp = np.asarray(inMat1[:,0])
	kq = np.asarray(inMat1[:,1])
	kr = np.asarray(inMat2[:,0])
	ks = np.asarray(inMat2[:,1])
	vMat = system.assym(kp,kq,kr,ks)
	#print vMat
	return vMat

#def tempFock(inMat1,inMat2):
#	dim1 = np.shape(inMat1)[0]
#	dim2 = np.shape(inMat2)[0]
#	fMat = np.zeros((dim1,dim2))
#	for row in range(0,dim1):
#		for col in range(0,dim2):
#			temp = system.f(inMat1[row,0]) + system.f(inMat1[row,1]) - system.f(inMat2[col,0]) - system.f(inMat2[col,1])
#			if temp > eps:
#				fMat[row,col] = temp
#	return fMat



#########################################################################################

#in this function, I've just c/p-ed what I used to have into pp and hh because it *seemingly* works
#to whomever might want to learn from this function, don't. It's difficult to understand all the steps
#and i suggest you rather find someone to explain how blocking works and how to do it, and then try to make your own
#seriously, even i have problems re-reading this, but it's written as it is simply because it's the fastest method i
#can think of. sorry if that's a bummer for you
def makeIntMatrix():	#just fyi, this function took me half a day, I'm not fluent in numpy functions
	blockArrays_pp_temp = np.zeros(((NS-N)**2,3),dtype=np.int)						#i really need to start using Np = NS-N
	blockArrays_hh_temp = np.zeros((N**2,3),dtype=np.int)
	#dtype=np.int is important as these are indices, and few functions are forgiving with double-type indices
	
	#only need two indices due to matrix diagonal symmetry
	index = 0
	for a in range(N,NS):
		for b in range(N,NS):
			blockArrays_pp_temp[index] = [system.kUnique2(a,b),a,b] 	#[k_u, a, b], i can probably use a+b*Np instead of "index", but whatevs
			index += 1
	index = 0
	for i in range(0,N):
		for j in range(0,N):
			blockArrays_hh_temp[index] = [system.kUnique2(i,j),i,j]
			index += 1
	print "made blockArrays_hh/pp"
	#print blockArrays_pp_temp.shape
	
	#blockArrays is now sorted in ascending orders of k_u
	blockArrays_pp = blockArrays_pp_temp[blockArrays_pp_temp[:,0].argsort()]		#we got this from the internets
	blockArrays_hh = blockArrays_hh_temp[blockArrays_hh_temp[:,0].argsort()]
	print "sorted blockArrays_hh/pp"
	
	#print blockArrays_pp
	
	#number of distinct identification numbers
	distinct_ku_pp = np.shape(np.unique(blockArrays_pp[:,0]))[0]
	distinct_ku_hh = np.shape(np.unique(blockArrays_hh[:,0]))[0]
	
	#array with the last index for each k_u, which seperates the blocks
	lastIndices_pp = np.zeros((distinct_ku_pp),dtype=np.int)
	lastIndices_hh = np.zeros((distinct_ku_hh),dtype=np.int)
	
	#BEWARE: for-loop notation might be wrong here
	index = 0
	for val in np.unique(blockArrays_pp[:,0]):									#val now takes all different values of k_u
		lastIndices_pp[index] = np.where(blockArrays_pp[:,0] == val)[0][-1]		#you get back a tuple where the first element is an array with indices for blockArrays==val)
		index += 1
	index = 0
	for val in np.unique(blockArrays_hh[:,0]):
		lastIndices_hh[index] = np.where(blockArrays_hh[:,0] == val)[0][-1]
		index += 1
	print "made lastIndices_hh/pp"
	
	#i think list is best, due to the varying size of the elements
	blockArraysNew_pp = []
	blockArraysNew_hh = []
	prev_i = 0
	#if i only wanted (a,b) instead of (k_u,a,b), then i put "1:" at the end of blockArraysNew_pp[-1] = blockArrays_pp[prev_i:i,:]
	#and 3 to 2 in blockArraysNew_pp.append(np.zeros((i+1-prev_i,3)))
	#however, i want k_u to do the if-test in the loop for V_hhpp
	for i in lastIndices_pp:
		blockArraysNew_pp.append(np.zeros((i+1-prev_i,3),dtype=np.int))
		blockArraysNew_pp[-1] = blockArrays_pp[prev_i:i+1,:]	#be careful with this indexing, i'm not too sure about it
		prev_i = i+1
	prev_i = 0
	for i in lastIndices_hh:
		blockArraysNew_hh.append(np.zeros((i+1-prev_i,3),dtype=np.int))
		blockArraysNew_hh[-1] = blockArrays_hh[prev_i:i+1,:]	#same here
		prev_i = i+1
	print "made blockArraysNew_hh/pp"
	
	intMat_pppp 	= []
	identify_pppp	= []
	intMat_hhpp 	= []
	identify_hhpp	= []
	FockMat 		= []
	identify_fock	= []
	#these loops ought to disappear once assym is vectorized (ps: you really need to vectorize, just sayin')
	
	for i in range(0,distinct_ku_pp):
		intMat_pppp.append(temp1(blockArraysNew_pp[i][:,1:]))
		identify_pppp.append(blockArraysNew_pp[i][0,0])
		
	for i in range(0,distinct_ku_hh):
		for j in range(0,distinct_ku_pp):
			if blockArraysNew_hh[i][0,0] == blockArraysNew_pp[j][0,0]:	#we need k_u to be the same, and i know of no better way atm
				intMat_hhpp.append(temp2(blockArraysNew_hh[i][:,1:],blockArraysNew_pp[j][:,1:]))
				identify_hhpp.append(blockArraysNew_hh[i][0,0])
				#FockMat.append(tempFock(blockArraysNew_hh[i][:,1:],blockArraysNew_pp[j][:,1:]))
				#identify_fock.append(blockArraysNew_hh[i][0,0])
	
	#intMat_pppp = temp1(blockArraysNew_pp[:][:,1:])
	#identify_pppp = blockArraysNew_pp[:][0,0]
	#
	#intMat_hhpp = temp2(blockArraysNew_hh[:][:,1:],blockArraysNew_pp[:][:,1:])
	#identify_hhpp = blockArraysNew_hh[:][0,0]
	#
	##FockMat.append(tempFock(blockArraysNew_hh[i][:,1:],blockArraysNew_pp[j][:,1:]))
	##identify_fock.append(blockArraysNew_hh[i][0,0])
	
	print "made intMat-s"
	
	return intMat_pppp, identify_pppp, intMat_hhpp, identify_hhpp#, FockMat, identify_fock, blockArraysNew_pp


#################################################################################



print NS
Vpppp, idpppp, Vhhpp, idhhpp = makeIntMatrix() #, F, idF, tempRet = makeIntMatrix()
print idpppp
#counter = 0
#for i in range(0,len(Vpppp)):
#	dim = np.shape(Vpppp[i])[0]
#	for j in range(0,dim):
#		for k in range(0,dim):
#			if Vpppp[i][j,k] != 0.:
#				print Vpppp[i][j,k], tempRet[i][j,1:], tempRet[i][k,1:]
#				counter += 1
#print counter
#Diagrams
def D3(i,j):	#L_a
	factor 		= 0.5
	#print np.shape(T[i])
	#print np.shape(np.transpose(Vpppp[j]))
	returnVal 	=  np.dot(T[i],np.transpose(Vpppp[j]))	#np.einsum('abcd,ijcd->ijab', Vpppp, T)
	return factor*returnVal
	
def diagrams():	#sum over all L_i and Q_i diagrams
	returnVal = []
	#print T[1][1]
	
	for i,x in enumerate(idT):
		for j,y in enumerate(idpppp):
			if x == y:
				#print x,y
				returnVal.append(D3(i,j))
	return returnVal

#########################################################################

T 		 = []
idT 	 = []
ECCD 	 = 0
ECCD_new = 0
for i in range(0,len(Vhhpp)):
	T.append(Vhhpp[i])	#need to figure out the division by F
	idT.append(idhhpp[i])
	#print np.dot(np.transpose(Vhhpp[i]),T[i])
	ECCD += 0.25*np.einsum('ab,ab',Vhhpp[i],T[i])
	#print ECCD

#for i in range(0,len(Vpppp)):
#	print idpppp[i]

#print len(T)
#print len(Vhhpp)
#ECCD = 0.25*np.einsum('klcd,klcd',Vhhpp,T)

while conFac > eps:
	D = diagrams()
	T_new = [] #this might be dodgy, considering numpy stuff
	for i in range(0,len(T)):	#HERE COMES THE PROBLEM OF ENSURING EQUAL K_U FOR Vhhpp AND Vpppp
		#print D[i]
		T_new.append(Vhhpp[i]+D[i])		#*denMat
	T = T_new

	for i in range(0,len(T)):
		ECCD_new += 0.25*np.einsum('ab,ab',Vhhpp[i],T[i])
		
	#ECCD_new = 0.25*np.einsum('klcd,klcd',Vhhpp,T)
	
	#print T
	#print ECCD_new
	
	conFac = abs(ECCD_new - ECCD)
	ECCD = ECCD_new
	ECCD_new = 0
	print ECCD

#print "CCD correlation:", ECCD/N, "Hatrees"
