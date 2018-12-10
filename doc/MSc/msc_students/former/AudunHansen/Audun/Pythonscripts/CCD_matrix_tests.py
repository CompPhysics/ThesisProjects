from numpy import *
from time import *

Nh = 10
Np = 7
Ns = Nh+Np

class tempbase():
    def __init__(self, Np, Nh):
        self.nstates = Np+Nh
        self.nparticles = Np
        self.nholes = Nh

class CCD():
    def __init__(self, bs):
        self.nstates = bs.nstates      #total number of states
        self.Np = bs.nparticles        #number of particle states
        self.Nh = self.nholes #number of hole states
        
        #set these properly up later
        self.Vhhhh = random.uniform(0,1,(Nh**2, Nh**2))
        self.Vhhpp = random.uniform(0,1,(Nh**2, Np**2))
        self.Vphhp = random.uniform(0,1,(Nh*Np, Nh*Np))
        self.Vhpph = random.uniform(0,1,(Nh*Np, Nh*Np))
        self.Vpppp = random.uniform(0,1,(Np**2, Np**2))

        self.Tpphh = random.uniform(0,1,(Np**2, Nh**2))
        
    #setup functions for each diagram in the CCD amplitude equation (L = linear, Q=quadratic in amplitude)
    def advance(self):
        t0 = clock()
        self.sL1()
        self.sL2()
        self.sL3()
        self.sQ1()
        self.sQ2()
        self.sQ3()
        self.sQ4()
        self.PL3 = self.P_ij(self.P_ab(self.L3))
        self.PQ2 = self.P_ij(self.Q2)
        self.PQ3 = self.P_ij(self.Q3)
        self.PQ4 = self.P_ab(self.Q4)
        
        self.Tpphh = .5*(self.L1+self.L2-self.PQ3-self.PQ4) + self.PL3 + .25 * self.Q1 + self.PQ2 
        #remember to divide elementwise by eps^ab_ij from the sp-energy
        print "Time spent on one iteration:", clock()-t0
    def sL1(self):
        #setup L1
        self.L1 = dot(self.Vpppp,self.Tpphh)
    def sL2(self):
        self.L2 = dot(self.Tpphh,self.Vhhhh)
    def sL3(self):
        self.L3 = self.TrL3(self.Vhpph, self.Tpphh)
    
    def sQ1(self):
        self.Q1 = dot(self.Tpphh,dot(self.Vhhpp, self.Tpphh))
    
    def sQ2(self):
        self.Q2 = self.TrQ2(self.Tpphh, self.Vhhpp)
        
    def sQ3(self):
        self.Q3 = self.TrQ3(Tpphh, Vhhpp)
    
    def sQ4(self):
        self.Q4 = self.TrQ4(Tpphh, Vhhpp)[a+b*Np, i + j*Nh] 
    
    #supporting functions
    def TrL3(self,V,T):
        TL3 = zeros((Nh*Np, Nh*Np))
        VL3 = zeros((Nh*Np, Nh*Np))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        TL3[a + i*Np, b + j*Np] = T[a + b*Np, i + j*Nh]
                        VL3[a + i*Np, b + j*Np] = V[i + b*Nh, a + j*Np]
        L3_ = dot(TL3, VL3)
        L3 = zeros((Np**2, Nh**2))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        L3[a + b*Np, i + j*Nh] = L3_[a + i*Np, b + j*Np]
        return L3  
          
    def TrQ2(self,T,V):
        #support for Q2 calculation
        TQ21 = zeros((Np*Nh, Np*Nh))
        VQ2  = zeros((Np*Nh, Np*Nh))
        TQ22 = zeros((Np*Nh, Np*Nh))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        TQ21[a + i*Np, b + j*Np] = T[a + b*Np, i + j*Nh]
                        TQ22[b + j*Np, a + i*Np] = T[a + b*Np, i + j*Nh]
                        VQ2[a + i*Np, b + j*Np]  = V[i + j*Nh, a + b*Np]
        Q2_ = dot(TQ21,dot(VQ2, TQ22))
        Q2  = zeros((Np**2, Nh**2))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        Q2[a + b*Np, i + j*Nh] = Q2_[a + i*Np, b + j*Np]
        
        #return TQ21, VQ2, TQ22
        return Q2  
    def TrQ3(self,T,V):
        TQ31 = zeros((Nh, Nh*Np**2))
        TQ32 = zeros((Nh, Nh*Np**2))
        VQ3  = zeros((Nh*Np**2, Nh))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        TQ31[i,j+ a*Nh + b *Np*Nh] = T[a + b*Np, i + j*Nh]
                        TQ32[i, j + a*Nh + b *Np*Nh] = T[a + b*Np, i + j*Nh]
                        VQ3[i + a*Nh + b *Np*Nh, j]  = V[i+ j*Nh, b + a*Np]
        Q3_ = dot(TQ31,dot(VQ3, TQ32))
        Q3  = zeros((Np**2, Nh**2))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        Q3[a + b*Np, i + j*Nh] = Q3_[i, j + a*Nh + b*Nh*Np]    
        return Q3 
    def TrQ4(self,T,V):
        TQ41 = zeros((Np, Np*Nh**2))
        TQ42 = zeros((Np, Np*Nh**2))
        VQ4  = zeros((Np*Nh**2, Np))     
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        TQ41[a,b+ j*Np + i *Nh*Np] = T[a + b*Np, j + i*Nh]
                        #TQ42[b, i + j*Nh + a *Nh*Nh] = T[b + a*Np, i + j*Nh]
                        #VQ4[i + j*Nh + a *Nh*Nh, b]  = V[i+ j*Nh, a + b*Np]
                        TQ42[b, a + i*Np + j *Nh*Np] = T[b + a*Np, i + j*Nh]
                        VQ4[b + j*Np + i *Nh*Np, a]  = V[i+ j*Nh, b + a*Np]
        Q4_ = dot(TQ41,dot(VQ4, TQ42))
        Q4  = zeros((Np**2, Nh**2))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        Q4[a + b*Np, i + j*Nh] = Q4_[a, b + i*Np + j*Np*Nh]    
        return Q4 
        
    #Permuting functions       
    def P_ab(self,M):
        #permute P(a,b)
        Z = zeros((Np**2, Nh**2))
        for i in range(self.Nh):
            for j in range(self.Nh): 
                for a in range(self.Np):
                    for b in range(self.Np):
                        Z[a+b*Np, i+j*Nh] = M[a+b*Np, i+j*Nh] - M[b+a*Np, i+j*Nh] 
        return Z
    def P_ij(M):
        #permute P(i,j)
        Z = zeros((Np**2, Nh**2))
        for i in range(self.Nh):
            for j in range(self.Nh): 
                for a in range(self.Np):
                    for b in range(self.Np):
                        Z[a+b*Np, i+j*Nh] = M[a+b*Np, i+j*Nh] - M[a+b*Np, j+i*Nh] 
        return Z



Vhhhh = random.uniform(0,1,(Nh**2, Nh**2))
Vhhpp = random.uniform(0,1,(Nh**2, Np**2))
Vphhp = random.uniform(0,1,(Nh*Np, Nh*Np))
Vhpph = random.uniform(0,1,(Nh*Np, Nh*Np))
Vpppp = random.uniform(0,1,(Np**2, Np**2))

Tpphh = random.uniform(0,1,(Np**2, Nh**2))
#a,b,i,j = 1,1,1,0
a,b,i,j =1,4,2,4
v = 0.0
t0 = clock()
for c in range(Np):
    for d in range(Np):
        #print a,b,c,d,i,j,k,l
        #print Vpppp[a+b*Np,c+d*Np]
        #print Tpphh[c+d*Np, i+j*Nh]
        v += Vpppp[a+b*Np,c+d*Np]*Tpphh[c+d*Np, i+j*Nh]
print "time 1:", clock()-t0

t0 = clock()
L1 = dot(Vpppp,Tpphh)[a+b*Np, i+j*Nh] #L1 contribution
print "time 2:", clock()-t0
print "L1 (target):", v
print "L1 (dot)   :", L1
print "================="


#
# L2 Contrib
#

v = 0
vv = zeros(Nh**2)
for k in range(Nh):
    for l in range(Nh):
        v += Vhhhh[k + l*Nh, i + j*Nh]*Tpphh[a+b*Np, k + l*Nh]
        vv[l + k*Nh] = Vhhhh[k + l*Nh, i + j*Nh]*Tpphh[a+b*Np, k + l*Nh]

print "L2 vv:", sum(vv)
print "L2 (target):", v
print "L2 (dot)   :", dot(Tpphh,Vhhhh)[a+b*Np, i+j*Nh]
print "================="

#L3 contrib

#Transformation
def TrL3_(M):
    Z = zeros((len(M), len(M[0])))
    #Z = zeros((Nh*Np, Nh*Np))
    for k in range(Nh):
        for c in range(Np):
            for j in range(Nh):
                for b in range(Np):
                    Z[c+k*Np,b + j*Np] = M[k+b*Nh,c+j*Np] #k+b*Nh,c+j*Np

    return Z

def TrT3(M):
    #for amplitudes in L3 contributions
    #Z = zeros((len(M), len(M[0])))
    Z = zeros((Np*Nh, Np*Nh))
    for k in range(Nh):
        for c in range(Np):
            for j in range(Nh):
                for b in range(Np):
                    Z[b + j*Np,k+c*Nh] = M[c+b*Np,i+j*Nh] #Tpphh[a+c*Np, i+k*Nh]

    return Z

def Ittr(M):
    Z = zeros((len(M), len(M[0])))
    for k in range(Nh):
        for c in range(Np):
            for j in range(Nh):
                for b in range(Np):
                    Z[k+b*Nh,c+j*Np] = M[k+c*Nh, b + j*Np]
    return Z               


v = 0
vv = zeros(Np*Nh)
vvv = 0
#VV = TrL3(Vhpph)
#TT = TrT3(Tpphh)
tt = zeros(Nh*Np)
NN = 0
for k in range(Nh):
    for c in range(Np):
        #print a,b,c,d,i,j,k,l
        #+= <kb||cj>t^ac_ik

        v += Vhpph[k+b*Nh,c+j*Np]
        vvv += Vhpph[k+b*Nh,c+j*Np]*Tpphh[a+c*Np, i+k*Nh]
        #vv[NN] = Vhpph[k+b*Nh,c+j*Np]
        vv[NN] = Vhpph[k+b*Nh,c+j*Np]
        tt[NN] = Tpphh[a+c*Np, i+k*Nh]
        NN+=1
        """
        v += Vhpph[k+b*Np,c+j*Nh]
        vvv += Vhpph[k+b*Np,c+j*Nh]*Tpphh[a+c*Np, i+k*Nh]
        vv[NN] = Vhpph[k+b*Np,c+j*Nh]
        tt[NN] = Tpphh[a+c*Np, i+k*Nh]
        NN+=1
        """
def TrL3(V,T):
    TL3 = zeros((Nh*Np, Nh*Np))
    VL3 = zeros((Nh*Np, Nh*Np))
    for i in range(Nh):
        for j in range(Nh): 
            for a in range(Np):
                for b in range(Np):
                    TL3[a + i*Np, b + j*Np] = T[a + b*Np, i + j*Nh]
                    VL3[a + i*Np, b + j*Np] = V[i + b*Nh, a + j*Np]
    L3_ = dot(TL3, VL3)
    L3 = zeros((Np**2, Nh**2))
    for i in range(Nh):
        for j in range(Nh): 
            for a in range(Np):
                for b in range(Np):
                    L3[a + b*Np, i + j*Nh] = L3_[a + i*Np, b + j*Np]
    return L3
                        

print " ...."


#vv = sum(TV[:,b+j*Np])
"""
print "tt:", tt, sum(tt)
print "..."
print "TT:", TT, sum(TT[:,a+i*Np])
print "..."
print "vv:", vv, sum(vv)
print "..."
print "VV:", VV, VV[:,b+j*Np]
print "TV:", sum(VV[:,b+j*Np]*TT[:,a+i*Np])
"""
print "TARGET(L3):", vvv
#print "Dot:", dot(VV,TT)[a+b*Np,i+j*Nh]
#print "Dot for L3:", dot(VV.T, TT)[i+j*Nh,a+b*Np]
print "Dot for L3:", TrL3(Vhpph, Tpphh)[a+b*Np,i+j*Nh]

#L3 = Ittr(dot(Ttr(Vhpph), Tpphh))
#print L3
#print "L3:", v, vv,  L3[a+b*Np, i+j*Nh]
print "...."

##
## Testing Q-terms
##

v2 = zeros((Nh**2, Np**2))

for c in range(Np):
    for d in range(Np):
        v2
v1  = 0
for k in range(Nh):
    for l in range(Nh): 
        for c in range(Np):
            for d in range(Np):
                v1 += Tpphh[a+b*Np, k+l*Nh]*Vhhpp[k+l*Nh,c+d*Np]*Tpphh[c+d*Np, i+j*Nh]

#VTT = dot(dot(Tpphh,Vhhpp), Tpphh)[a+b*Np, i+j*Nh]
Q1 = dot(Tpphh,dot(Vhhpp, Tpphh))[a+b*Np, i+j*Nh]  #the contributions from Q1
print "Q1(target) :", v1
print "Q1 (dot)   :", Q1


#Testing Q2
print "==="
v = 0
for k in range(Nh):
    for l in range(Nh): 
        for c in range(Np):
            for d in range(Np):
                v += Vhhpp[k+l*Nh,c+d*Np]*Tpphh[a+c*Np, i + k*Nh]*Tpphh[b+d*Np, j+l*Nh]

#transforming matrices

def TrQ2(T,V):
    #support for Q2 calculation
    TQ21 = zeros((Np*Nh, Np*Nh))
    VQ2  = zeros((Np*Nh, Np*Nh))
    TQ22 = zeros((Np*Nh, Np*Nh))
    for i in range(Nh):
        for j in range(Nh): 
            for a in range(Np):
                for b in range(Np):
                    TQ21[a + i*Np, b + j*Np] = T[a + b*Np, i + j*Nh]
                    TQ22[b + j*Np, a + i*Np] = T[a + b*Np, i + j*Nh]
                    VQ2[a + i*Np, b + j*Np]  = V[i + j*Nh, a + b*Np]
    Q2_ = dot(TQ21,dot(VQ2, TQ22))
    Q2  = zeros((Np**2, Nh**2))
    for i in range(Nh):
        for j in range(Nh): 
            for a in range(Np):
                for b in range(Np):
                    Q2[a + b*Np, i + j*Nh] = Q2_[a + i*Np, b + j*Np]
    
    #return TQ21, VQ2, TQ22
    return Q2

print "Q2 (target):", v
#T1, Q, T2 = TrQ2(Tpphh, Vhhpp)
print "Q2 (dot)   :", TrQ2(Tpphh, Vhhpp)[a+b*Np, i+j*Nh]
print "==="


#Testing Q3
v =0.0
for k in range(Nh):
    for l in range(Nh): 
        for c in range(Np):
            for d in range(Np):
                v += Vhhpp[k+l*Nh,c+d*Np]*Tpphh[d+c*Np, i + k*Nh]*Tpphh[a+b*Np, l+j*Nh]

def TrQ3(T,V):
    TQ31 = zeros((Nh, Nh*Np**2))
    TQ32 = zeros((Nh, Nh*Np**2))
    VQ3  = zeros((Nh*Np**2, Nh))
    for i in range(Nh):
        for j in range(Nh): 
            for a in range(Np):
                for b in range(Np):
                    TQ31[i,j+ a*Nh + b *Np*Nh] = T[a + b*Np, i + j*Nh]
                    TQ32[i, j + a*Nh + b *Np*Nh] = T[a + b*Np, i + j*Nh]
                    VQ3[i + a*Nh + b *Np*Nh, j]  = V[i+ j*Nh, b + a*Np]
                    
                    #TQ31[i, b + j*Np + a *Np*Nh] = T[a + b*Np, i + j*Nh]
                    #TQ32[i, b +a*Np + j *Np*Nh] = T[a + b*Np, i + j*Nh]
                    #VQ3[b + j*Np + a*Np*Nh, i]  = V[i+ j*Nh, a + b*Np]
                    
                    #TQ31[i, a + j*Np + b *Np*Nh] = T[a + b*Np, i + j*Nh]
                    #TQ32[i, a + j*Np + b *Np*Nh] = T[a + b*Np, i + j*Nh]
                    #VQ3[a + j*Nh + b *Np*Nh, i]  = V[i+ j*Nh, a + b*Np]
    Q3_ = dot(TQ31,dot(VQ3, TQ32))
    Q3  = zeros((Np**2, Nh**2))
    for i in range(Nh):
        for j in range(Nh): 
            for a in range(Np):
                for b in range(Np):
                    Q3[a + b*Np, i + j*Nh] = Q3_[i, j + a*Nh + b*Nh*Np]    
    return Q3

print "Q3 (target):", v
print "Q3 (dot)   :", TrQ3(Tpphh, Vhhpp)[i+j*Nh, a+b*Np]
print "Q3 (dot)   :", TrQ3(Tpphh, Vhhpp)[a+b*Np, i+j*Nh]


#test Q4
print "==="

v = 0.0
for k in range(Nh):
    for l in range(Nh): 
        for c in range(Np):
            for d in range(Np):
                v += Vhhpp[k+l*Nh,c+d*Np]*Tpphh[a+c*Np, l + k*Nh]*Tpphh[d+b*Np, i+j*Nh]
 
def TrQ4(T,V):
    TQ41 = zeros((Np, Np*Nh**2))
    TQ42 = zeros((Np, Np*Nh**2))
    VQ4  = zeros((Np*Nh**2, Np))     
    for i in range(Nh):
        for j in range(Nh): 
            for a in range(Np):
                for b in range(Np):
                    TQ41[a,b+ j*Np + i *Nh*Np] = T[a + b*Np, j + i*Nh]
                    #TQ42[b, i + j*Nh + a *Nh*Nh] = T[b + a*Np, i + j*Nh]
                    #VQ4[i + j*Nh + a *Nh*Nh, b]  = V[i+ j*Nh, a + b*Np]
                    TQ42[b, a + i*Np + j *Nh*Np] = T[b + a*Np, i + j*Nh]
                    VQ4[b + j*Np + i *Nh*Np, a]  = V[i+ j*Nh, b + a*Np]
    Q4_ = dot(TQ41,dot(VQ4, TQ42))
    Q4  = zeros((Np**2, Nh**2))
    for i in range(Nh):
        for j in range(Nh): 
            for a in range(Np):
                for b in range(Np):
                    Q4[a + b*Np, i + j*Nh] = Q4_[a, b + i*Np + j*Np*Nh]    
    return Q4        
print "Q4 (target):", v
print "Q4 (dot)   :", TrQ4(Tpphh, Vhhpp)[a+b*Np, i + j*Nh]    