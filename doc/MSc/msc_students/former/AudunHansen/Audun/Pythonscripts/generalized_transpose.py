from numpy import *
from scipy.sparse import csr_matrix, coo_matrix
from matplotlib.pyplot import *

#a sparse matrix multiplication algorithm

class flexmat():
    def __init__(self, data, p,q,r,s, Np, Nq, Nr, Ns):
        self.data = data
        self.p = p
        self.q = q
        self.r = r
        self.s = s
        self.Np = Np
        self.Nq = Nq
        self.Nr = Nr
        self.Ns = Ns
    #Casting routines
    def pq_rs(self):
        try:
            return self.Mpq_rs
        except:
            #coo_matrix((M.data, (a+ b*self.Np, i + j*self.Nh)), shape=(self.Np**2,self.Nh**2)).tocsr() 
            self.Mpq_rs = coo_matrix((data, (self.p+self.q*self.Np, self.r+self.s*self.Nr)), shape = (self.Np*self.Nq, self.Nr*self.Ns)).tocsr()
            return self.Mpq_rs
    def pr_qs(self):
        try:
            return self.Mpr_qs
        except:
            #coo_matrix((M.data, (a+ b*self.Np, i + j*self.Nh)), shape=(self.Np**2,self.Nh**2)).tocsr() 
            self.Mpr_qs = coo_matrix((data, (self.p+self.r*self.Np, self.q+self.s*self.Nq)), shape = (self.Np*self.Nr, self.Nq*self.Ns)).tocsr()
            return self.Mpr_qs
    
    def p_qrs(self):
        pass
        
#NOTE: Why do we need the sp-mat in the first place? Is it possible to just perform the multiplication right on the "data object above?
#Does the armadillo-sp_mat-dot perform better than a iterative procedure? Why?
        

Np = 30
Nh = 14
Nv = 600

data = random.randint(0, 4, Nv)
a = random.randint(0,Np**2, Nv)%Np
b = random.randint(0,Np**2, Nv)//Np
i = random.randint(0,Nh**2, Nv)%Nh
j = random.randint(0,Nh**2, Nv)//Nh

T = flexmat(data, a,b,i,j,Np,Np,Nh,Nh)

def spmult(val1, val2, p1, p2,  q1, q2):
    N1 = len(val1)
    N2 = len(val2)
    p3 = p1
    q3 = q2
    val3 = zeros(len(p3))
    for i in range(N1):
        for e in range(N2):
            if q1[i]==p2[e]:
                val3[i] += val1[i]*val2[e]
    return val3, p3, q3

val, p, q = spmult(data, data, a+b*Np, i+j*Nh, a+b*Np, i+j*Nh) 

figure(1)
imshow(coo_matrix((data, (p,q)), shape = (Np**2, Np**2)).toarray())
show()

figure(2)
T2 = T.pq_rs().dot(T.pq_rs().T).toarray()
print T2
imshow(T2)
show()

#imshow(T.pr_qs().toarray())
#show()


        