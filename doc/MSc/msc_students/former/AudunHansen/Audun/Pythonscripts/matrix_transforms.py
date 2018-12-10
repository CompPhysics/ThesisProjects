from numpy import *

class transforms():
    def __init__(self, Np, Nh):
        self.Np = Np
        self.Nh = Nh
        
        
        self.Tab = zeros((Np**2), dtype = int)
        for b in range(Np):
            self.Tab[range(b*Np,(b+1)*Np)] = range(b,Np**2, Np)

        self.Tai = zeros((Np*Nh), dtype = int)
        for a in range(Nh):
            self.Tai[range(a*Np,(a+1)*Np)] = range(b,Np*Nh, Nh)
        
        self.Tij = zeros((Nh**2), dtype = int)
        for i in range(Nh):
            self.Tij[range(i*Nh,(i+1)*Nh)] = range(i,Nh**2,Nh)
    def abij2baji(self,M):
        return M[self.Tab][:, self.Tij]
    def abij2baji_(self,M): 
        Z = zeros((Np*Np,Nh*Nh))
        for i in range(Nh):
            for j in range(Nh):
                for a in range(Np):
                    for b in range(Np):
                        Z[a + b*Np, i + j * Nh] = M[b + a*Np, j + i*Nh]
        return Z
    def abij2aibj_(self,M):
        Z = zeros((Nh*Np,Nh*Np))
        Z2 = zeros((Nh*Np,Nh*Np))
        for i in range(self.Nh):
            for j in range(self.Nh):
                for a in range(self.Np):
                    for b in range(self.Np):
                        #TL3[a + i*Np, b + j*Np] = T[a + b*Np, i + j*Nh]
                        Z[a + i*self.Np, b + j*self.Np] = M[a + b*Np, i + j*Nh]
        return Z
        
    def abij2aibj(self,M):
        Z = zeros((Nh*Np,Nh*Np))
        print M[self.Tai, :]
        #return M[self.Tai][:, self.Tai]       

class transformer():
    def __init__(self, Np, Nh):
        self.Np = Np
        self.Nh = Nh
        self.init_index_reorganization()
    def init_index_reorganization(self):
        Np = self.Np
        Nh = self.Nh
        self.Pai_bjA = zeros((Np*Nh,Nh*Np), dtype = int)
        self.Pai_bjB = zeros((Np*Nh,Nh*Np), dtype = int)
        
        self.L3_A = zeros((Np**2,Nh**2), dtype = int)
        self.L3_B = zeros((Np**2,Nh**2), dtype = int)
        
        self.Pbj_aiA = zeros((Np*Nh,Nh*Np), dtype = int)
        self.Pbj_aiB = zeros((Np*Nh,Nh*Np), dtype = int)
        
        self.Pi_jabA = zeros((Nh,Nh*Np**2), dtype = int)
        self.Pi_jabB = zeros((Nh,Nh*Np**2), dtype = int)
        
        self.Pa_bjiA = zeros((Np,Np*Nh**2), dtype = int)
        self.Pa_bjiB = zeros((Np,Np*Nh**2), dtype = int)
        
        self.Pb_aijA = zeros((Np,Np*Nh**2), dtype = int)
        self.Pb_aijB = zeros((Np,Np*Nh**2), dtype = int)
        
        self.Q2_A = zeros((Np**2, Nh**2), dtype = int)
        self.Q2_B = zeros((Np**2, Nh**2), dtype = int)
        
        self.Q3_A = zeros((Np**2, Nh**2), dtype = int)
        self.Q3_B = zeros((Np**2, Nh**2), dtype = int)
        
        self.Q4_A = zeros((Np**2, Nh**2), dtype = int)
        self.Q4_B = zeros((Np**2, Nh**2), dtype = int)
        
        
        for i in range(Nh):
            for j in range(Nh):
                for a in range(Np):
                    for b in range(Np):
                        self.Pai_bjA[a + i*Np, b + j*Np] = a + b*Np
                        self.Pai_bjB[a + i*Np, b + j*Np] = i + j*Nh
                        
                        self.L3_A[a + b*Np, i + j*Nh] = a + i*Np
                        self.L3_B[a + b*Np, i + j*Nh] = b + j*Np
                        
                        self.Pbj_aiA[b + j*Np, a + i*Np] = a + b*Np
                        self.Pbj_aiB[b + j*Np, a + i*Np] = i + j*Nh
                        
                        self.Pi_jabA[i,j+ a*Nh + b *Np*Nh] = a + b*Np #TQ31&TQ32
                        self.Pi_jabB[i,j+ a*Nh + b *Np*Nh] = i + j*Nh
                        
                        self.Pa_bjiA[a,b+ j*Np + i *Nh*Np] = a + b*Np
                        self.Pa_bjiB[a,b+ j*Np + i *Nh*Np] = j + i*Nh
                        
                        self.Pb_aijA[b, a + i*Np + j *Nh*Np] = a + b*Np
                        self.Pb_aijB[b, a + i*Np + j *Nh*Np] = i + j*Nh
                        
                        self.Q2_A[a + b*Np, i + j*Nh] = a + i*Np
                        self.Q2_B[a + b*Np, i + j*Nh] = b + j*Np
                        
                        self.Q3_A[a + b*Np, i + j*Nh] = i
                        self.Q3_B[a + b*Np, i + j*Nh] = j + a*Nh + b*Nh*Np
    
                        self.Q4_A[a + b*Np, i + j*Nh] = a
                        self.Q4_B[a + b*Np, i + j*Nh] = b + i*Np + j*Np*Nh
                        
                        #Q3[a + b*Np, i + j*Nh] = Q3_[i, j + a*Nh + b*Nh*Np]
                        #Q2[a + b*Np, i + j*Nh] = Q2_[a + i*Np, b + j*Np]

        
        
Np = 3
Nh = 2

h = transformer(Np,Nh)

pphh = random.randint(0,10,(Np**2,Nh**2))
vv = random.randint(0,10,(Np*Nh,Nh*Np))
#pphh = zeros((Np**2,Nh**2))

#n = 0
#for a in range(Np**2):
#    for i in range(Nh**2):
#        pphh[a,i] = n
#        n+= 1
Z = zeros((Np,Np*Nh**2), dtype = int)
for i in range(Nh):
    for j in range(Nh):
        for a in range(Np):
            for b in range(Np):
               Z[a,b+ j*Np + i *Nh*Np] = pphh[a + b*Np, j + i*Nh]
                
def PermuteIndexes(array, perm):
    return array[ix_(*(perm[:s] for s in array.shape))]

print pphh

print "==="
print Z
print "==="
print pphh[h.Pa_bjiA, h.Pa_bjiB]
#print "PP:", pphh[PaibjA, PaibjB]
a,b,i,j = 0,2,1,1
L3_t = 0.0
for k in range(Nh):
    for c in range(Np):
        L3_t += vv[k + b*Nh, c + j*Np]*pphh[a + c*Np, i + k*Nh]
print "L3_linear:", L3_t
#print pphh[pphh[0]]

ktrA = zeros((Nh*Np, Nh*Np), dtype = int)
ktrB = zeros((Nh*Np, Nh*Np), dtype = int)
ktr = zeros((Nh*Np, Nh*Np))

ptrA = zeros((Nh*Np, Nh*Np), dtype = int)
ptrB = zeros((Nh*Np, Nh*Np), dtype = int)
ptr = zeros((Nh*Np, Nh*Np))

L3A = zeros((Np**2, Nh**2), dtype = int)
L3B = zeros((Np**2, Nh**2), dtype = int)

for i_ in range(Nh):
    for j_ in range(Nh):
        for a_ in range(Np):
            for b_ in range(Np):
                ktrA[j_ + b_*Nh, a_ + i_*Np] = j_ + a_*Nh #L3 transformation indices
                ktrB[j_ + b_*Nh, a_ + i_*Np] = b_ + i_*Np
                ktr[j_ + b_*Nh, a_ + i_*Np] = vv[j_ + a_*Nh, b_ + i_*Np]
                
                ptrA[a_ + i_*Np, j_ + b_*Nh] = a_ + b_*Np
                ptrB[a_ + i_*Np, j_ + b_*Nh] = i_ + j_*Nh
                ptr[a_ + i_*Np, j_ + b_*Nh] = pphh[a_ + b_*Np, i_+j_*Nh]
                
                L3A[a_ + b_*Np, i_ + j_*Nh] = a_ + i_*Np
                L3B[a_ + b_*Np, i_ + j_*Nh] = b_ + j_*Np
print "==="
print ktr
KT = vv[ktrA, ktrB]
print "==="
print ptr
PT =  pphh[ptrA, ptrB]

L3_ = dot(PT,KT)
L3 = zeros((Np**2, Nh**2))


for i_ in range(Nh):
    for j_ in range(Nh):
        for a_ in range(Np):
            for b_ in range(Np):
                L3[a_ + b_*Np, i_ + j_*Nh] = L3_[a_ + i_*Np, b_ + j_*Np]
                
print "L3_linear:", L3_t
print L3
print L3_[L3A, L3B][a + b*Np, i + j*Nh]


#print pphh[range(0,3)]
#print h.abij2aibj_(pphh)
#print h.abij2aibj(pphh)
#print h.abij2aibj(pphh)


"""
TL3[a + i*Np, b + j*Np] = T[a + b*Np, i + j*Nh]
L3[a + b*Np, i + j*Nh] = L3_[a + i*Np, b + j*Np]
TQ21[a + i*Np, b + j*Np] = T[a + b*Np, i + j*Nh]
TQ22[b + j*Np, a + i*Np] = T[a + b*Np, i + j*Nh]
TQ31[i,j+ a*Nh + b *Np*Nh] = T[a + b*Np, i + j*Nh]
TQ32[i,j+ a*Nh + b *Np*Nh] = T[a + b*Np, i + j*Nh]
TQ41[a,b+ j*Np + i *Nh*Np] = T[a + b*Np, j + i*Nh]
TQ42[b, a + i*Np + j *Nh*Np] = T[b + a*Np, i + j*Nh]

Q4[a + b*Np, i + j*Nh] = Q4_[a, b + i*Np + j*Np*Nh]
"""