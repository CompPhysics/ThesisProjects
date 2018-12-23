from numpy import *

Np = 3
Nh = 2

T = random.randint(0,10,(Np**2, Nh**2))
V = random.randint(0,10,(Nh**2, Np**2))

Tt = zeros((Np*Nh, Np*Nh))

TA = zeros((Nh*Nh*Np, Np))
TB = zeros((Nh*Nh*Np, Np))

T0A = zeros((Np**2, Nh**2))
T0B = zeros((Np**2, Nh**2))

for i in range(Nh):
    for j in range(Nh):
        for a in range(Np):
            for b in range(Np):
                TA[b + j*Np + i *Nh*Np, a] = i+ j*Nh
                TB[b + j*Np + i *Nh*Np, a] = b + a*Np
                

                #self.VQ4[b + j*Np + i *Nh*Np, a]  = self.Vhhpp[i+ j*Nh, b + a*Np]
                
                
                #self.L3_A[a + b*Np, i + j*Nh] = a + i*Np
                #self.L3_B[a + b*Np, i + j*Nh] = b + j*Np 
                
                #self.Pb_aijA[b, a + i*Np + j *Nh*Np] = a + b*Np
                #self.Pb_aijB[b, a + i*Np + j *Nh*Np] = i + j*Nh
                #T0A[a + b*Np, j + i*Nh] = a + b*Np
                #T0B[a + b*Np, j + i*Nh] = i + j*Nh
                
         
#print outer(range(10), ones(9, dtype = int))
#unaltered indexing

#Optimization madness. Do not freak out! Or even better: do not read!
ij2ij = outer(ones(Np**2, dtype = int), range(Nh**2))
ab2ab = outer(range(Np**2), ones(Nh**2, dtype = int))

#unaltered
ab_ij2ab_ij = [outer(range(Np**2), ones(Nh**2, dtype = int)), outer(ones(Np**2, dtype = int), range(Nh**2))]

#self.Pba_ijA[b + a*Np, i + j*Nh] = a + b*Np
#self.Pba_ijB[b + a*Np, i + j*Nh] = i + j*Nh
abba = arange(0,Np**3).reshape(Np**2,Np)[:,0]%(Np**2-1)
abba[Np**2 - 1] = Np**2 - 1
e1 = kron(ones((1,Nh**2), dtype = int).T, abba).T  
e2 =  arange(0,Nh**3).reshape(Nh**2,Nh)[:,0]%(Nh**2- 1)
e2[Nh**2 - 1] = Nh**2 - 1
e2 =  kron(ones((1,Np**2), dtype = int).T, e2)
ab_ij2ba_ji = [e1,e2]  
#self.Pba_jiA[b + a*Np, j + i*Nh] = a + b*Np
#self.Pba_jiB[b + a*Np, j + i*Nh] = i + j*Nh

ijji = arange(0,Nh**3).reshape(Nh**2,Nh)[:,0]%(Nh**2-1)
ijji[Nh**2 - 1] = Nh**2 - 1
ab_ij2ab_ji = [ab2ab, outer(ones(Np**2, dtype = int), ijji)]        
#self.Pab_jiA[a + b*Np, j + i*Nh] = a + b*Np
#self.Pab_jiB[a + b*Np, j + i*Nh] = i + j*Nh

abba = arange(0,Np**3).reshape(Np**2,Np)[:,0]%(Np**2-1)
abba[Np**2 - 1] = Np**2 - 1
ab_ij2ba_ij = [outer(abba, ones(Nh**2, dtype = int)),ij2ij]

ab_ij2aj_bi =[kron(ones((Nh,Nh), dtype = int),arange(0,Np**2).reshape(Np,Np).T), kron(arange(0,Nh**2).reshape(Nh,Nh),ones((Nh,Nh), dtype = int))]
ab_ij2ai_bj = [kron(ones((Nh,Nh), dtype = int), arange(0,Np**2).reshape(Np,Np).T), kron(arange(0,Nh**2).reshape(Nh,Nh).T,ones((Np,Np), dtype = int))]
#ab_ij2ai_bj = []        
#self.Pai_bjA[a + i*Np, b + j*Np] = a + b*Np
#self.Pai_bjB[a + i*Np, b + j*Np] = i + j*Nh



e1= kron(kron(ones((Nh), dtype =int).T, arange(0,Np**2).reshape(Np,Np).T), ones((Nh,1), dtype = int))
#print " "
e2 = kron(ones((Np,1), dtype = int),kron(arange(0,Nh**2).reshape(Nh,Nh).T, ones((Np),dtype = int)))
#print kron(arange(0,Nh**2).reshape(Nh,Nh).T, ones((Np),dtype = int))
ab_ij2ia_bj = [e1,e2]
#self.Pia_bjA[i + a*Nh, b + j*Np] = a + b*Np
#self.Pia_bjB[i + a*Nh, b + j*Np] = i + j*Nh

ab_ij2ia_jb = [kron(kron(arange(0,Np**2).reshape(Np,Np).T, ones((Nh),dtype = int)),ones((Nh,1), dtype = int)),kron(ones((Np,1), dtype = int),kron(ones((Np), dtype =int).T, arange(0,Nh**2).reshape(Nh,Nh).T))]

#self.Pia_bjA[i + a*Nh, b + j*Np] = a + b*Np
ab_ij2ai_jb = [kron(ones((Nh,1), dtype = int),kron(arange(0,Np**2).reshape(Np,Np).T, ones((Nh),dtype = int))), kron(ones((Nh**2,1), dtype = int).T,kron(arange(0,Nh**2).reshape(Nh,Nh).T, ones((1,Np),dtype = int).T))]
#self.Pai_jbA[a + i*Np, j + b*Nh] = a + b*Np
#self.Pai_jbB[a + i*Np, j + b*Nh] = i + j*Nh

ab_ij2bj_ai = [kron(ones((Nh,Nh), dtype = int), arange(0,Np**2).reshape(Np,Np).T).T,kron(arange(0,Nh**2).reshape(Nh,Nh).T,ones((Np,Np), dtype = int)).T]
#self.Pbj_aiA[b + j*Np, a + i*Np] = a + b*Np
#self.Pbj_aiB[b + j*Np, a + i*Np] = i + j*Nh

ab_ij2i_jab = [kron(arange(0,Np**2),ones((Nh,Nh), dtype =int)),kron(ones((1,Np**2), dtype = int), arange(0,Nh**2).reshape(Nh,Nh).T) ]        
#self.Pi_jabA[i,j+ a*Nh + b *Np*Nh] = a + b*Np #TQ31&TQ32
#self.Pi_jabB[i,j+ a*Nh + b *Np*Nh] = i + j*Nh

e1 =  kron(ones((1,Nh**2), dtype = int), arange(0,Np**2).reshape(Np,Np).T)
e2 = kron(arange(0,Nh**2),ones((Np,Np), dtype =int))
ab_ij2a_bji = [e1,e2]  
#self.Pa_bjiA[a,b+ j*Np + i *Nh*Np] = a + b*Np
#self.Pa_bjiB[a,b+ j*Np + i *Nh*Np] = j + i*Nh

e1 =  kron(ones((1,Nh**2), dtype = int), arange(0,Np**2).reshape(Np,Np))
e2 = arange(0,Nh**3).reshape(Nh**2,Nh)[:,0]%(Nh**2-1)
e2[Nh**2 - 1] = Nh**2 - 1 
e2 = kron(e2,ones((Np,Np), dtype =int))
ab_ij2b_aji = [e1,e2]        
#self.Pb_aijA[b, a + i*Np + j *Nh*Np] = a + b*Np
#self.Pb_aijB[b, a + i*Np + j *Nh*Np] = i + j*Nh

e1 =  kron(ones((1,Nh**2), dtype = int), arange(0,Np**2).reshape(Np,Np))
e2 = arange(0,Nh**3).reshape(Nh**2,Nh)[:,0]%(Nh**2-1)
e2[Nh**2 - 1] = Nh**2 - 1 
e2 = arange(0,Nh**2)
e2 = kron(e2,ones((Np,Np), dtype =int))
#print e1
#print " "
#print e2
ab_ij2b_aij = [e1,e2]  

e1 =  kron(ones((Np**2,Nh), dtype = int), arange(0,Nh).reshape(Nh))
e2 = arange(0,Nh*Np**2).reshape(Np**2,Nh) #[:,0]%(Nh**2-1)
e2 = kron(e2,ones((1,Nh), dtype =int))
i_jab2ab_ij = [e1,e2]        
#self.Q3_A[a + b*Np, i + j*Nh] = i
#self.Q3_B[a + b*Np, i + j*Nh] = j + a*Nh + b*Nh*Np

e1 =  kron(ones((Nh**2, Np), dtype = int), arange(0,Np).reshape(Np)).T
e2 = arange(0,Np*Nh**2).reshape(Nh**2,Np)
e2 = kron(e2,ones((1,Np), dtype =int)).T
a_bij2ab_ij = [e1,e2]
#self.Q4_A[a + b*Np, i + j*Nh] = a
#self.Q4_B[a + b*Np, i + j*Nh] = b + i*Np + j*Np*Nh

e1 = kron(ones((Np, Nh)), arange(Np*Nh).reshape(Nh,Np).T)
e2 = kron(arange(0,Np*Nh).reshape(Nh,Np).T,ones((Np,Nh),dtype = int))
ai_bj2ab_ij = [e1,e2]
#self.L3_A[a + b*Np, i + j*Nh] = a + i*Np
#self.L3_B[a + b*Np, i + j*Nh] = b + j*Np     
#self.Q2_A[a + b*Np, i + j*Nh] = a + i*Np
#self.Q2_B[a + b*Np, i + j*Nh] = b + j*Np

e1 = kron( arange(Np**2), ones((1, Nh**2), dtype = float).T).T
e2 = arange(0,Nh**3).reshape(Nh**2,Nh)[:,0]%(Nh**2-1)
e2[Nh**2 - 1] = Nh**2 - 1 
e2 = outer(ones((Np**2), dtype = int), e2)

ab_ji2ab_ij = [e1,e2]


e1 =arange(0,Np*Nh).reshape(Np,Nh).T
e1 = kron(e1,ones((1,Np), dtype = int).T)
e1 = kron(ones((Nh), dtype = int),e1)
e2 = arange(0,Np*Nh).reshape(Nh,Np).T
e2 = kron(ones((1, Nh), dtype = int),kron(e2,ones((1,Np), dtype = int)).T).T
vl3 = [e1,e2]




vq2 = [kron(arange(0,Nh**2).reshape(Nh,Nh).T,ones((Np,Np), dtype = int)),kron(ones((Nh,Nh), dtype = int), arange(0,Np**2).reshape(Np,Np).T)]



e1 =  kron(ones((Np*Np), dtype = int), arange(0,Nh**2).reshape(Nh,Nh)).T
e2 = arange(0,Np**3).reshape(Np**2,Np)[:,0]%(Np**2-1)
e2[Np**2 - 1] = Np**2 - 1 
e2 = kron(kron(ones((Nh,1), dtype = int),e2),ones((Nh),dtype =int)).T
vq3 = [e1,e2]


e1 = arange(0,Nh**3).reshape(Nh**2,Nh)[:,0]%(Nh**2-1)
e1[Nh**2 - 1] = Nh**2 - 1 
e1 =  kron(e1,ones((Np,Np), dtype = int)).T
e2 =  kron(ones((1,Nh**2), dtype = int).T, arange(0,Np**2).reshape(Np,Np).T)

vq4 = [e1,e2]  
print e1
print " "
print e2
print "===T0==="
print TA
print "..."
print TB
print "========"
#A = array([[0,3,6],[1,4,7],[2,5,8]])
g = True
if True:
    for i in range(Nh):
        for j in range(Nh):
            for a in range(Np):
                for b in range(Np):
                    if TA[b + j*Np + i *Nh*Np, a] != vq4[0][b + j*Np + i *Nh*Np, a]:
                        print "error"
                        g = False
                        break
                    if TB[b + j*Np + i *Nh*Np, a] != vq4[1][b + j*Np + i *Nh*Np, a]:
                        print "error"
                        g = False
                        break
                if g == False:
                    break
            if g == False:
                break
        if g == False:
            break




#print ravel(A).reshape(3,3)
#print linspace(0,Np**2 - 1,Np**2).reshape(Np,Np).T
#print kron(ones((Nh,Nh), dtype = int), arange(0,Np**2).reshape(Np,Np).T)
#print kron(arange(0,Nh**2).reshape(Nh,Nh).T,ones((Np,Np), dtype = int))
#print kron(ones((2,2)),A)
#print outer( A, array([[1,0],[0,1]]), array([[1,0],[0,1]]))