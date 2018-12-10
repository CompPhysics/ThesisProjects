from numpy import *
from time import *
from matplotlib.pyplot import *
from scipy.sparse import csr_matrix, coo_matrix

class electronbasis():
    def __init__(self, N, rs, Nparticles):
        self.rs = rs
        self.states = []
        self.nstates = 0
        self.nparticles = Nparticles
        Nm = int(sqrt(N) + 1)
        self.Nm = Nm
        #Creating the basis
        for x in range(-Nm, Nm):
            for y in range(-Nm, Nm):
                for z in range(-Nm,Nm):
                    e = x*x + y*y + z*z
                    if e <=N:
                        self.states.append([e, x,y,z, 1])
                        self.states.append([e, x,y,z,-1])
                        self.nstates += 2
        self.states.sort() #Sorting the basis in increasing energy

        self.L3 = (4*pi*self.nparticles*self.rs**3)/3.0
        self.L2 = self.L3**(2/3.0)
        self.L = pow(self.L3, 1/3.0)
        
        for i in range(self.nstates):
            self.states[i][0] *= 2*(pi**2)/self.L**2 #Multiplying in the missing factors in the single particle energy
        self.states = array(self.states) #converting to array to utilize vectorized calculations
    def hfenergy(self, nParticles):
        #Calculating the HF-energy (reference energy)
        e0 = 0.0
        if nParticles<=self.nstates:
            for i in range(nParticles):
                e0 += self.h(i,i)
                for j in range(nParticles):
                    if j != i:
                        e0 += .5*self.v(i,j,i,j)
        else:
            #Safety for cases where nParticles exceeds size of basis
            print "Not enough basis states."
            
        return e0
                
    def h(self, p,q):
        #Return single particle energy
        return self.states[p,0]*(p==q)
    def veval(self, p,q,r,s):
        #A test for evaluating the two-body interaction
        val = ""
        if self.kdplus(p,q,r,s):
            val+= "kdplus "
            if self.kdspin(p,r):
                val += "Direct[kdspin_pr "
                if self.kdspin(q,s):
                    val += "kdspin_qs "
                    if self.kdwave(p,r) != 0:
                        val += "kdwave!=0 "
                        val += str(self.absdiff2(r,p))
                val += "] "

            if self.kdspin(p,s):
                val += "Exchange[kdspin_pr "
                if self.kdspin(q,r):
                    val += "kdspin_qs "
                    if self.kdwave(p,s) != 0:
                        val += "kdwave!=0 "
                        val += str(self.absdiff2(s,p))
                val += "] "
        return val
    def vevalHF(self, N):
        #Evaluation of all expressions of two-body contributions to the HF-energy
        for i in range(N):
            for j in range(N):
                if i!= j:
                    print "<",i,j,"|",i,j,"> =",self.veval(i,j,i,j)
    def V(self, kp,kq,kr,ks):
        #k = (energy, kx, ky, kz, ms)
        # Vectorized interaction
        #
        
        #kplus
        kdplus = (kp[1,:]+kq[1,:]==kr[1,:]+ks[1,:])*(kp[2,:]+kq[2,:]==kr[2,:]+ks[2,:])*(kp[3,:]+kq[3,:]==kr[3,:]+ks[3,:])*4*pi/self.L3#d_k+k k+k
        #print "kdplus:", kdplus
        
        kdspin1 = (kp[4,:]==kr[4,:])*(kq[4,:]==ks[4,:])*1
        kdwave1 = abs((kp[1,:]==kr[1,:])*(kp[2,:]==kr[2,:])*(kp[3,:]==kr[3,:])-1)
        #print "kdspin1:", kdspin1
        #print "kdwave1:", kdwave1
        absdiff2_1 = ((kr[1,:]-kp[1,:])**2+(kr[2,:]-kp[2,:])**2+(kr[3,:]-kp[3,:])**2) #absdiff2
        term1=(4.0*absdiff2_1*pi**2)/self.L2
        term1[term1==0] = 1
        
        kdspin2 = (kp[4,:]==ks[4,:])*(kq[4,:]==kr[4,:])*1
        kdwave2 = abs((kp[1,:]==ks[1,:])*(kp[2,:]==ks[2,:])*(kp[3,:]==ks[3,:])-1)
        #print "kdspin2:",kdspin2
        #print "kdwave2:",kdwave2
        absdiff2_2 = ((ks[1,:]-kp[1,:])**2+(ks[2,:]-kp[2,:])**2+(ks[3,:]-kp[3,:])**2) #absdiff2
        #print absdiff2_2
        term2=(4.0*absdiff2_2*pi**2)/self.L2
        term2[term2==0] = 1
        
        return kdplus*(kdspin1*kdwave1/term1 - kdspin2*kdwave2/term2)
        
    def v(self,p,q,r,s):
        #Two body interaction
        #To optimize bottleneck: vectorize this function ! (remove if-tests)
        val = 0
        terms = 0.0
        kdpl = self.kdplus(p,q,r,s)
        if kdpl != 0:
            val = 4*pi/self.L3
            term1 = 0.0
            term2 = 0.0
            if self.kdspin(p,r)*self.kdspin(q,s)==1:
                if self.kdwave(p,r) != 1.0:
                    term1=(4*self.absdiff2(r,p)*pi**2)/self.L2
                    terms += 1.0/term1
    
            if self.kdspin(p,s)*self.kdspin(q,r)==1:
                if self.kdwave(p,s) != 1.0:
                    term2=(4*self.absdiff2(s,p)*pi**2)/self.L2
                    terms -= 1.0/term2
        return val*terms

    
    #The following is a series of kroenecker deltas used in the two-body interactions. 
    #Run kd_integrity() to ensure that they work as intended.
    def kdi(self,a,b):
        #Kroenecker delta integer
        return 1.0*(a==b)
    def kda(self,a,b):
        #Kroenecker delta array
        d = 1.0
        #print a,b,
        for i in range(len(a)):
            d*=(a[i]==b[i])
        return d
    def kdfullplus(self,p,q,r,s):
        #Kroenecker delta wavenumber p+q,r+s
        return self.kda(self.states[p][1:5]+self.states[q][1:5],self.states[r][1:5]+self.states[s][1:5])
    def kdplus(self,p,q,r,s):
        #Kroenecker delta wavenumber p+q,r+s
        return self.kda(self.states[p][1:4]+self.states[q][1:4],self.states[r][1:4]+self.states[s][1:4])
    def kdspin(self,p,q):
        #Kroenecker delta spin
        return self.kdi(self.states[p][4], self.states[q][4])
    def kdwave(self,p,q):
        #Kroenecker delta wavenumber
        return self.kda(self.states[p][1:4],self.states[q][1:4])
    def absdiff2(self,p,q):
        val = 0.0
        for i in range(1,4):
            val += (self.states[p][i]-self.states[q][i])*(self.states[p][i]-self.states[q][i])
        #if val == 0:
        #    print "div0"
        return val
    def kd_integrity(self):
        #test integrity of kroenecker deltas
        print "Array KD            :", self.kda([0,1,2], [0,1,2]) == True
        print "Integer KD          :", self.kdi(1,1) == True
        print "Opposite spin       :", self.kdspin(0,1) == False
        print "Equal spin          :", self.kdspin(1,1) == True
        print "Wavenumber equal    :", self.kdwave(1,0) == True
        print "Wavenumber not equal:", self.kdwave(1,2) == False
    def liststates(self):
        for i in range(self.nstates):
            print self.states[i]


class tempbase():
    def __init__(self, Np, Nh):
        self.nstates = Np+Nh
        self.nparticles = Np
        self.nholes = Nh

class CCD():
    def __init__(self, bs):
        self.bs = bs
        self.nstates = bs.nstates            #total number of states
        self.Nh = bs.nparticles              #number of hole states (conflicting naming should be resolved in class electrongas)
        self.Np = self.nstates-bs.nparticles #number of particle states

        self.Vhhhh = csr_matrix((self.Nh**2, self.Nh**2))
        self.Vhhpp = csr_matrix((self.Nh**2, self.Np**2))
        self.Vphhp = csr_matrix((self.Nh*self.Np, self.Nh*self.Np))
        self.Vhpph = csr_matrix((self.Nh*self.Np, self.Nh*self.Np))
        self.Vpppp = csr_matrix((self.Np**2, self.Np**2))
        self.Vpphh = csr_matrix((self.Np**2, self.Nh**2))
        self.Tpphh = csr_matrix((self.Np**2, self.Nh**2))
        self.Epphh = zeros((self.Np**2, self.Nh**2))
        
        self.setup_matrices_optimized()
        
    ################################################
    ##
    ## MAIN PROGRAM ROUTINES 
    ##
    ################################################

    def setup_matrices_optimized(self):
        #Fill inn all matrices
        #This is probably the bottleneck right now, should apply symmetries to oprimize
        
        Nh = self.Nh
        Np = self.Np
        
        
        
        #alternate setup for Epphh
        E = self.bs.states[:,0]
        pp = arange(Np**2)
        hh = arange(Nh**2)
        a = pp%Np
        b = pp//Np
        i = hh%Nh
        j = hh//Nh

        ij = kron(ones((Np**2,1)), E[i] + E[j])
        ab = kron(ones((Nh**2,1)), E[a+Nh] + E[b+Nh])

        self.Epphh = ij - ab.T
       
        
        
        
        t0 = clock()
        
        
        
        """
        for i in range(Nh):
            for j in range(i,Nh):
                for a in range(Np):
                    for b in range(a,Np):

                        
                        val = self.bs.v(a+Nh,i,j,b+Nh)
                        if val != 0:
                            self.Vphhp[a + i*Np, j + b*Nh] = val
                            self.Vphhp[b + j*Np, i + a*Nh] = val 

                        val = self.bs.v(j,a+Nh,b+Nh,i)
                        if val != 0:
                            self.Vphhp[a + j*Np, i + b*Nh] = val
                            self.Vhpph[j + a*Nh, b + i*Np] = val
                            self.Vphhp[b + i*Np, j + a*Nh] = val
                            self.Vhpph[i + b*Nh, a + j*Np] = val


                        val = self.bs.v(a+Nh,b+Nh,i,j)
                        #eps = self.bs.h(i,i) + self.bs.h(j,j) -self.bs.h(a+Nh,a+Nh) - self.bs.h(b+Nh,b+Nh)
                        eps = self.Epphh[a + b*Np, i + j*Nh]
                        #if self.Epphh[a + b*Np, i +j*Nh] != val:
                        #    #print val, self.Epphh[a + b*Np, i +j*Np]
                        #    self.Epphh[a + b*Np, i + j*Nh] = eps
                        #    self.Epphh[a + b*Np, j + i*Nh] = eps
                        #    self.Epphh[b + a*Np, i + j*Nh] = eps
                        #    self.Epphh[b + a*Np, j + i*Nh] = eps
                        if val != 0:

                            self.Vpphh[a + b*Np, i + j*Nh] = val
                            self.Vpphh[a + b*Np, j + i*Nh] = -val
                            self.Vpphh[b + a*Np, i + j*Nh] = -val
                            self.Vpphh[b + a*Np, j + i*Nh] = val                        
                            self.Vhhpp[i + j*Nh, a + b*Np] = val
                            self.Vhhpp[j + i*Nh, b + a*Np] = val
                            self.Vhhpp[j + i*Nh, a + b*Np] = -val
                            self.Vhhpp[i + j*Nh, b + a*Np] = -val
                    
                        
                            self.Tpphh[a + b*Np, i + j*Nh] = val/eps
                            self.Tpphh[a + b*Np, j + i*Nh] = -val/eps
                            self.Tpphh[b + a*Np, i + j*Nh] = -val/eps
                            self.Tpphh[b + a*Np, j + i*Nh] = val/eps
        """
        t1 = clock()
        
        print "Time spent setting up amplitudes and eps:", t1-t0
        
        t0 = clock()
        B = blocks(tb)
        self.Vhhhh = B.Vhhhh
        self.Vpppp = B.Vpppp
        self.Vhhpp = B.Vhhpp
        self.Vpphh = B.Vpphh
        self.Vhpph = B.Vhpph
        self.Vphhp = B.Vphhp
        t1 = clock()
        
        self.Tpphh = csr_matrix(self.Vpphh/self.Epphh)
        
        print "Time spent setting up interactions:", t1-t0
        """


        optiv = optimV(self.bs)
        self.Vhhhh = csr_matrix(optiv.Vhhhh)
        self.Vpppp = csr_matrix(optiv.Vpppp)
        t2 = clock()
        print "Time spent on setting up hhpp terms:", t1-t0
        print "Time spent on setting up pppp and hhhh terms:", t2-t1
        """
        
        
        """

        t0 = clock()
        
        for i in range(Nh):
            for j in range(i,Nh):
                for k in range(Nh):
                    for l in range(k,Nh):
                        val = self.bs.v(i,j,k,l)
                        if val!=0:
                            self.Vhhhh[i + j*Nh, k+ l*Nh] = val
                            self.Vhhhh[j + i*Nh, l+ k*Nh] = val
                            self.Vhhhh[j + i*Nh, k+ l*Nh] = -val
                            self.Vhhhh[i + j*Nh, l+ k*Nh] = -val
        t1 = clock()

        for a in range(Np):
            for b in range(a,Np):
                for c in range(Np):
                    for d in range(c,Np):
                        val = self.bs.v(a+Nh,b+Nh,c+Nh,d+Nh)
                        if val!= 0:
                            self.Vpppp[a + b*Np, c+ d*Np] = val 
                            self.Vpppp[b + a*Np, d+ c*Np] = val
                            self.Vpppp[b + a*Np, c+ d*Np] = -val
                            self.Vpppp[a + b*Np, d+ c*Np] = -val                
        t2 = clock()
        print "Time spent setting up Vhhhh (iteratively):", t1-t0
        print "Time spent setting up Vpppp (iteratively):", t2-t1
        """

        #Aligned matrices for L3, Q2, Q3 and Q4 multiplications
        self.VL3 = self.perm_ind_ib_aj2ai_bj(self.Vhpph)       
        self.VQ2 = self.perm_ind_ij_ab2ai_bj(self.Vhhpp)
        self.VQ3 = self.perm_ind_ij_ba2iab_j(self.Vhhpp)
        self.VQ4 = self.perm_ind_ij_ba2bji_a(self.Vhhpp)
        

        
    def advance(self):
        #Main loop, run this to advance solution one iteration

        #setup linear contributions
        self.sL1()
        self.sL2()
        self.sL3()
        
        #setup quadratic contributions
        self.sQ1()
        self.sQ2()
        self.sQ3()
        self.sQ4()
        
        #permute contributions
        self.PL3 = self.L3 - self.perm_ind_ba_ij(self.L3) - self.perm_ind_ab_ji(self.L3) + self.perm_ind_ba_ji(self.L3)
        self.PQ2 = self.Q2 - self.perm_ind_ab_ji(self.Q2)
        self.PQ3 = self.Q3 - self.perm_ind_ab_ji(self.Q3)
        self.PQ4 = self.Q4 - self.perm_ind_ba_ij(self.Q4)
        
        #Sum all contributions
        
        self.Tpphh = (self.Vpphh + .5*(self.L1 + self.L2) + self.PL3 + .25*self.Q1 + self.PQ2 - .5*(self.PQ3 + self.PQ4))/self.Epphh
        
        #self.sp_epsdiv(self.Tpphh)

        #calculate energy
        self.energy()
        
        #Update UI
        print "        Correlation energy:", self.C_energy

        #Update amplitudes (have been temporarily dense due to division above)
        self.Tpphh = csr_matrix(self.Tpphh)

        
    def e0_(self):
        Np = self.Np
        Nh = self.Nh
        e0 = 0.0
        for i in range(Nh):
            for j in range(Nh):
                for a in range(Np):
                    for b in range(Np):    
                        e0 += self.Vhhpp[i+j*Nh, a+b*Np]*self.Tpphh[a + b*Np, i+j*Nh]   
        return e0     
    def energy(self):
        Np = self.Np
        Nh = self.Nh                    
        C = self.Vhhpp.dot(self.Tpphh)
        N = len(C)
        #self.C_energy = .25*sum(C.diagonal())
        self.C_energy = .25*sum(C[range(0,N), range(0,N)])
    
    def sp_epsdiv(self, M):
        #sparse matrix energy division
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        print self.bs.states[:,0][i] + self.bs.states[:,0][j] - self.bs.states[:,0][a] - self.bs.states[:,0][b]
        M.data/=(self.bs.states[:,0][i] + self.bs.states[:,0][j] - self.bs.states[:,0][a] - self.bs.states[:,0][b])
        
    #######################################  
    ##
    ##   SPARSE PERMUTATION ROUTINES
    ##   A set of functions that efficiently permutes and reshapes sparse matrix representations of rank 4 tensors
    ##
    #######################################

    def unpack_indptr(self,indptr):
        #Unpack row-compressed indices
        I =zeros(indptr[-1], dtype = int)
        for i in range(len(indptr)-1):
            I[indptr[i]:indptr[i+1]] = i
        return I
    
    def perm_ind_ai_bj(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        return coo_matrix((M.data, (a + i*self.Np, b + j*self.Np)), shape=(self.Np*self.Nh, self.Np*self.Nh)).tocsr()

    def perm_ind_ia_bj(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        return coo_matrix((M.data, (i + a*self.Nh, b + j*self.Np)), shape=(self.Np*self.Nh, self.Np*self.Nh)).tocsr()

    def perm_ind_bj_ai(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        return coo_matrix((M.data, (b + j*self.Np, a + i*self.Np)), shape=(self.Np*self.Nh, self.Np*self.Nh)).tocsr()

    def perm_ind_ai_jb(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        return coo_matrix((M.data, (a + i*self.Np, j + b*self.Nh)), shape=(self.Np*self.Nh, self.Np*self.Nh)).tocsr()

    def perm_ind_ba_ij(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        return coo_matrix((M.data, (b + a*self.Np, i + j*self.Nh)), shape=(self.Np**2, self.Nh**2)).tocsr()

    def perm_ind_ab_ji(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        return coo_matrix((M.data, (a + b*self.Np, j + i*self.Nh)), shape=(self.Np**2, self.Nh**2)).tocsr()
        
    def perm_ind_ba_ji(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        return coo_matrix((M.data, (b + a*self.Np, j + i*self.Nh)), shape=(self.Np**2, self.Nh**2)).tocsr()        

    def perm_ind_i_jab(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        return coo_matrix((M.data, (i, j + a*self.Nh+ b*self.Nh*self.Np)), shape=(self.Nh, self.Nh*self.Np**2)).tocsr()        

    def perm_ind_a_bji(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        return coo_matrix((M.data, (a, b + j*self.Np+ i*self.Nh*self.Np)), shape=(self.Np, self.Np*self.Nh**2)).tocsr()   

    def perm_ind_b_aji(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        return coo_matrix((M.data, (b, a + j*self.Np+ i*self.Nh*self.Np)), shape=(self.Np, self.Np*self.Nh**2)).tocsr()   
    
    def perm_ind_ij_ab2ai_bj(self,M):
        #Sparse permutations
        #print M.shape
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        i,j,a,b = rows%self.Nh, rows//self.Nh,cols%self.Np, cols//self.Np
        return coo_matrix((M.data, (a + i*self.Np, b + j*self.Np)), shape=(self.Np*self.Nh, self.Np*self.Nh)).tocsr()
     
    def perm_ind_ij_ba2iab_j(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        i,j,b,a = rows%self.Nh, rows//self.Nh,cols%self.Np, cols//self.Np
        return coo_matrix((M.data, (i + a*self.Nh + b*self.Nh*self.Np, j)), shape=(self.Np*self.Nh*self.Np, self.Nh)).tocsr()
        
    def perm_ind_ij_ba2bji_a(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        i,j,b,a = rows%self.Nh, rows//self.Nh,cols%self.Np, cols//self.Np
        return coo_matrix((M.data, (b + j*self.Np + i*self.Np*self.Nh, a)), shape=(self.Np*self.Nh**2, self.Np)).tocsr()     
     
    #def perm_ind_ai_bj2ab_ij(self,M):
    #    #Sparse permutations
    #    cols, rows = M.indices, self.unpack_indptr(M.indptr)
    #    a,i,b,j = rows%self.Np, rows//self.Np,cols%self.Np, cols//self.Np
    #    return coo_matrix((M.data, (a + b*self.Np,i + j*self.Nh)), shape=(self.Np**2, self.Nh**2)).tocsr()     
     
    def perm_ind_ai_bj2a_bji(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,i,b,j = rows%self.Np, rows//self.Np,cols%self.Np, cols//self.Np
        return coo_matrix((M.data, (a, b + j*self.Np + i*self.Np*self.Nh)), shape=(self.Np, self.Np*self.Nh**2)).tocsr()

    def perm_ind_ib_aj2ai_bj(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        i,b,a,j = rows%self.Nh, rows//self.Nh,cols%self.Np, cols//self.Np
        return coo_matrix((M.data, (a + i*self.Np, b + j*self.Np)), shape=(self.Np*self.Nh, self.Np*self.Nh)).tocsr()

    def perm_ind_ai_bj2ab_ij(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,i,b,j = rows%self.Np, rows//self.Np,cols%self.Np, cols//self.Np
        return coo_matrix((M.data, (a+ b*self.Np, i + j*self.Nh)), shape=(self.Np**2,self.Nh**2)).tocsr()        

    def perm_ind_a_bij2ab_ij(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a = rows
        b = cols%self.Np
        i = ((cols-b)/self.Np)%self.Nh
        j = ((cols-b)/self.Np)//self.Nh        
        return coo_matrix((M.data, (a+ b*self.Np, i + j*self.Nh)), shape=(self.Np**2,self.Nh**2)).tocsr()          

    def perm_ind_i_jab2ab_ij(self, M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        i = rows
        j = cols%self.Nh
        a = ((cols-j)/self.Nh)%self.Np
        b = ((cols-j)/self.Nh)//self.Np        
        return coo_matrix((M.data, (a+ b*self.Np, i + j*self.Nh)), shape=(self.Np**2,self.Nh**2)).tocsr()                        

    ##############################################
    ##
    ##  Contributions to the CCD amplitude
    ##  As in S-B, the contributions is defined as linear L (t) and quadratic Q (tt) 
    ##  The CCD amplitude equation then reads
    ##  Tpphh = (v + L1 + L2 + L3 + Q1 + Q2 + Q3 + Q4)/eps
    ##
    ##############################################
        
    def sL1(self):
        self.L1 = self.Vpppp.dot(self.Tpphh)
    def sL2(self):
        self.L2 = (self.Vhhhh.T.dot(self.Tpphh.T)).T
    def sL3(self):
        self.L3 = self.TL3()
    
    def sQ1(self):
        self.Q1 = ((self.Vhhpp.dot(self.Tpphh)).T.dot(self.Tpphh.T)).T
    def sQ2(self):
        self.Q2 = self.TQ2(self.Tpphh, self.Vhhpp)
    def sQ3(self):
        self.Q3 = self.TQ3(self.Tpphh, self.Vhhpp)
    def sQ4(self):
        self.Q4 = self.TQ4(self.Tpphh, self.Vhhpp)#[a+b*Np, i + j*Nh] 
                          
           
    def TL3(self):
        #The L3 Contribution
        self.TL3_ = self.perm_ind_ai_bj(self.Tpphh)
        L3_ = (self.VL3.T.dot(self.TL3_.T)).T
        return self.perm_ind_ai_bj2ab_ij(L3_)     
                              
    def TQ2(self,T,V):
        #The Q2 contrubution
        TQ21 = self.perm_ind_ai_bj(self.Tpphh)
        TQ22 = self.perm_ind_bj_ai(self.Tpphh)
        Q2_ = (self.VQ2.dot(TQ22).T.dot(TQ21.T)).T
        return  self.perm_ind_ai_bj2ab_ij(Q2_)

    def TQ3(self,T,V):
        #The Q3-contrubution
        TQ31 = self.perm_ind_i_jab(self.Tpphh)
        Q3_ = (self.VQ3.dot(TQ31).T.dot(TQ31.T)).T
        return self.perm_ind_i_jab2ab_ij(Q3_)
        
    def TQ4(self,T,V):
        #The Q4 contribution
        TQ41 = self.perm_ind_a_bji(self.Tpphh)
        Q4_ = (self.VQ4.dot(TQ41).T.dot(TQ41.T)).T
        return self.perm_ind_a_bij2ab_ij(Q4_)

class optimV():
    def __init__(self, bs):
        self.bs = bs
        self.Np = bs.nstates-bs.nparticles
        self.Nh = bs.nparticles
        self.Ns = bs.nstates
        self.Nm = self.bs.Nm #Max possible momentum  
        self.Nm2 = self.Nm**2
        self.Nm3 = self.Nm**3
        self.Vpppp = zeros((self.Np**2, self.Np**2))
        self.Vhhhh = zeros((self.Nh**2, self.Nh**2))
        self.Vhhpp = zeros((self.Nh**2, self.Np**2))
        self.Vhpph = zeros((self.Nh*self.Np, self.Nh*self.Np))
        #self.Vhhhh = zeros((self.Nh**2, self.Nh**2))
        self.setup_pppp()
        self.setup_hhhh()
        self.setup_hhpp()
        self.setup_hpph()
        #setup hp
        #seutp ph
        

                    
    def ident(self,v):
        #A unique identifying integer for the momentum combinations
        return v[0] + v[1]*self.Nm + v[2]*self.Nm2 + v[3]*self.Nm3                            
    def setup_pppp(self):
        t0 = clock()
        Np = self.Np
        
        combs_pp = 20000*ones((Np**2,Np**2), dtype = int) #arbitrary large number since identifier will include zeros
        idents = zeros((Np**2))
        for p in range(Np):
            for q in range(Np):
                v = self.bs.states[p+self.Nh][1:5]+self.bs.states[q+self.Nh][1:5]
                iv =  self.ident(v)
                combs_pp[p + q*Np, :] = iv  #this one should not be zero, as most elements in array is already zero, or?
                idents[p+q*Np] = iv
        spectrum = unique(idents)
        combs_pp[combs_pp!=combs_pp.T]=20000 #identify each pair of quantum numbers sharing the same added momentum
        self.combs_pp = combs_pp
        t1 = clock()
        print "Time spent determining unique sortings:", t1-t0
        self.setup_Vpppp()
        
    def setup_Vpppp(self):
        for P in range(self.Np**2):
            for Q in range(P,self.Np**2):
                if self.combs_pp[P,Q] != 20000:
                    a,b = P%self.Np, P//self.Np
                    c,d = Q%self.Np, Q//self.Np
                    #if self.ident(self.bs.states[a][1:5]+self.bs.states[b][1:5])==self.ident(self.bs.states[c][1:5]+self.bs.states[d][1:5]):
                    #if product(self.bs.states[a+Nh][1:5]+self.bs.states[b+Nh][1:5])==product(self.bs.states[c+Nh][1:5]+self.bs.states[d+Nh][1:5]):
                    val = self.bs.v(a+self.Nh,b+self.Nh,c+self.Nh,d+self.Nh)
                    self.Vpppp[P,Q] = val
                    self.Vpppp[Q,P] = val
                    
    def setup_hhhh(self):
        Nh = self.Nh
        
        combs_hh = 20000*ones((Nh**2,Nh**2), dtype = int) #arbitrary large number since identifier will include zeros
        idents = zeros((Nh**2))
        for p in range(Nh):
            for q in range(Nh):
                v = self.bs.states[p][1:5]+self.bs.states[q][1:5]
                iv =  self.ident(v)
                combs_hh[p + q*Nh, :] = iv  #this one should not be zero, as most elements in array is already zero, or?
                idents[p+q*Nh] = iv
        spectrum = unique(idents)
        combs_hh[combs_hh!=combs_hh.T]=20000 #identify each pair of quantum numbers sharing the same added momentum
        self.combs_hh = combs_hh   
        self.setup_Vhhhh()  
           
    def setup_Vhhhh(self):
        for P in range(self.Nh**2):
            for Q in range(P,self.Nh**2):
                if self.combs_pp[P,Q] != 20000:
                    i,j = P%self.Nh, P//self.Nh
                    k,l = Q%self.Nh, Q//self.Nh
                    #if self.ident(self.bs.states[a][1:5]+self.bs.states[b][1:5])==self.ident(self.bs.states[c][1:5]+self.bs.states[d][1:5]):
                    #if product(self.bs.states[a+Nh][1:5]+self.bs.states[b+Nh][1:5])==product(self.bs.states[c+Nh][1:5]+self.bs.states[d+Nh][1:5]):
                    val = self.bs.v(i,j,k,l)
                    self.Vpppp[P,Q] = val
                    self.Vpppp[Q,P] = val    
    def setup_hhpp(self):
        Nh = self.Nh
        Np = self.Np
        combs_hh = 20000*ones((Nh**2,Np**2), dtype = int) #arbitrary large number since identifier will include zeros
        combs_pp = 20000*ones((Nh**2,Np**2), dtype = int) #arbitrary large number since identifier will include zeros
        #idents = zeros((Nh*Np))
        for p in range(Np):
            for q in range(Np):
                v = self.bs.states[p+Nh][1:5]+self.bs.states[Nh + q][1:5]
                iv =  self.ident(v)
                combs_pp[:,p + q*Np] = iv  #this one should not be zero, as most elements in array is already zero, or?
                #idents[p+q*Np] = iv
        for p in range(Nh):
            for q in range(Nh):
                v = self.bs.states[p][1:5]+self.bs.states[q][1:5]
                iv =  self.ident(v)
                combs_hh[p + q*Nh, :] = iv  #this one should not be zero, as most elements in array is already zero, or?
                #idents[p+q*Np] = iv
        #spectrum = unique(idents)
        combs_hh[combs_pp!=combs_hh]=20000 #identify each pair of quantum numbers sharing the same added momentum
        self.combs_hp = combs_hh   
        self.setup_Vhhpp()  
           
    def setup_Vhhpp(self):
        for P in range(self.Nh**2):
            for Q in range(self.Np**2):
                if self.combs_hp[P,Q] != 20000:
                    #Run trough common setup routine here
                    i,j = P%self.Nh, P//self.Nh
                    a,b = Q%self.Np, Q//self.Np
                    val = self.bs.v(i,j,a+self.Nh,b+self.Nh)
                    self.Vhhpp[P,Q] = val
                    #self.Vpppp[Q,P] = val  
        self.Vpphh = self.Vhhpp.T
        
    def setup_hpph(self):
        Nh = self.Nh
        Np = self.Np
        combs_hp = 20000*ones((Nh*Np,Nh*Np), dtype = int) #arbitrary large number since identifier will include zeros
        combs_ph = 20000*ones((Nh*Np,Nh*Np), dtype = int) #arbitrary large number since identifier will include zeros
        idents = zeros((Nh**2))
        for p in range(Nh):
            for q in range(Np):
                v = self.bs.states[p][1:5]+self.bs.states[q+Nh][1:5]
                iv =  self.ident(v)
                combs_hp[p + q*Nh, :] = iv  #this one should not be zero, as most elements in array is already zero, or?
                
                combs_ph[:,q + p*Np] = iv  #this one should not be zero, as most elements in array is already zero, or?
                #idents[p+q*Nh] = iv
        #spectrum = unique(idents)
        combs_hp[combs_hp!=combs_ph]=20000 #identify each pair of quantum numbers sharing the same added momentum
        self.combs_hpph = combs_hp 
        self.setup_Vhpph()  
           
    def setup_Vhpph(self):
        for P in range(self.Nh*self.Np):
            for Q in range(self.Np*self.Nh):
                if self.combs_hpph[P,Q] != 20000:
                    i,a = P%self.Nh, P//self.Nh
                    b,j = Q%self.Np, Q//self.Np
                    #if self.ident(self.bs.states[a][1:5]+self.bs.states[b][1:5])==self.ident(self.bs.states[c][1:5]+self.bs.states[d][1:5]):
                    #if product(self.bs.states[a+Nh][1:5]+self.bs.states[b+Nh][1:5])==product(self.bs.states[c+Nh][1:5]+self.bs.states[d+Nh][1:5]):
                    val = self.bs.v(i,a,b,j)
                    self.Vhpph[P,Q] = val
                    #self.Vhpph[Q,P] = val            
        
class blocks():
    def __init__(self, bs):
        self.bs = bs
        self.Np = bs.nstates-bs.nparticles
        self.Nh = bs.nparticles
        self.Ns = bs.nstates
        self.Nm = self.bs.Nm #Max possible momentum
        
        self.Vhhhh = zeros((self.Nh**2, self.Nh**2))
        self.Vhhpp = zeros((self.Nh**2, self.Np**2))
        self.Vphhp = zeros((self.Nh*self.Np, self.Nh*self.Np))
        self.Vhpph = zeros((self.Nh*self.Np, self.Nh*self.Np))
        self.Vpppp = zeros((self.Np**2, self.Np**2))
        self.Vpphh = zeros((self.Np**2, self.Nh**2))

        self.Tpphh = zeros((self.Np**2, self.Nh**2))
        self.Epphh = zeros((self.Np**2, self.Nh**2))
        
        #self.setup_matrices_optimized()
        #self.Tpphh = random.uniform(0,1,(self.Np**2, self.Nh**2))
        self.setup_pppp()
        self.setup_hhhh()
        self.setup_hhpp()
        self.setup_hpph()
    def ident(self,v):
        #A unique identifying integer for the momentum combinations
        return v[0] + v[1]*self.bs.Nm + v[2]*self.bs.Nm**2 + v[3]*self.bs.Nm**3
        
    def setup_pppp(self):
        Np = self.Np
        combs_pp = 20000*ones((Np**2,Np**2), dtype = int) #arbitrary large number since identifier will include zeros
        idents = zeros((Np**2))
        for p in range(Np):
            for q in range(Np):
                v = self.bs.states[p+self.Nh][1:5]+self.bs.states[q+self.Nh][1:5]
                iv =  self.ident(v)
                combs_pp[p + q*Np, :] = iv  #this one should not be zero, as most elements in array is already zero, or?
                idents[p+q*Np] = iv
        combs_pp[combs_pp!=combs_pp.T]=20000 #identify each pair of quantum numbers sharing the same added momentum
        t = where(combs_pp!=2000)
        a = self.bs.states[t[0]%Np  + self.Nh].T
        b = self.bs.states[t[0]//Np + self.Nh].T
        c = self.bs.states[t[1]%Np  + self.Nh].T
        d = self.bs.states[t[1]//Np + self.Nh].T

        data = self.bs.V(a,b,c,d)
        #print data[data!=0]
        self.Vpppp = coo_matrix((data, (t[0], t[1])), shape=(self.Np**2, self.Np**2)).tocsr()
    def setup_hhhh(self):
        Np = self.Np
        Nh = self.Nh
        combs_hh = 20000*ones((Nh**2,Nh**2), dtype = int) #arbitrary large number since identifier will include zeros
        #idents = zeros((Nh**2))
        for p in range(Nh):
            for q in range(Nh):
                v = self.bs.states[p][1:5]+self.bs.states[q][1:5]
                iv =  self.ident(v)
                combs_hh[p + q*Nh, :] = iv  #this one should not be zero, as most elements in array is already zero, or?
                #idents[p+q*Nh] = iv
        combs_hh[combs_hh!=combs_hh.T]=20000 #identify each pair of quantum numbers sharing the same added momentum
        t = where(combs_hh!=2000)
        a = self.bs.states[t[0]%Nh ].T
        b = self.bs.states[t[0]//Nh].T
        c = self.bs.states[t[1]%Nh ].T
        d = self.bs.states[t[1]//Nh].T

        data = self.bs.V(a,b,c,d)
        #print data[data!=0]
        self.Vhhhh = coo_matrix((data, (t[0], t[1])), shape=(self.Nh**2, self.Nh**2)).tocsr()

    def setup_hpph(self):
        Np = self.Np
        Nh = self.Nh
        combs_hp = 20000*ones((Nh*Np,Nh*Np), dtype = int) #arbitrary large number since identifier will include zeros
        combs_ph = 20000*ones((Nh*Np,Nh*Np), dtype = int) #arbitrary large number since identifier will include zeros
        
        #idents = zeros((Nh**2))
        for p in range(Nh):
            for q in range(Np):
                v = self.bs.states[p][1:5]+self.bs.states[q+Nh][1:5]
                iv =  self.ident(v)
                combs_hp[p + q*Nh, :] = iv  #this one should not be zero, as most elements in array is already zero, or?
                combs_ph[:,q + p*Np ] = iv
                #idents[p+q*Nh] = iv
        combs_hp[combs_hp!=combs_ph.T]=20000 #identify each pair of quantum numbers sharing the same added momentum
        t = where(combs_hp!=2000)
        i = self.bs.states[t[0]%Nh ].T
        a = self.bs.states[t[0]//Nh + Nh].T
        b = self.bs.states[t[1]%Np  + Nh].T
        j = self.bs.states[t[1]//Np].T

        data = self.bs.V(i,a,b,j)
        #print data[data!=0]
        self.Vhpph = coo_matrix((data, (t[0], t[1])), shape=(self.Nh*Np, self.Nh*Np)).tocsr()
        self.Vphhp = self.Vhpph.T
   
    def setup_hhpp(self):
        Np = self.Np
        Nh = self.Nh
        combs_hh = 20000*ones((Nh*Nh,Np*Np), dtype = int) #arbitrary large number since identifier will include zeros
        combs_pp = 20000*ones((Nh*Nh,Np*Np), dtype = int) #arbitrary large number since identifier will include zeros
        
        #idents = zeros((Nh**2))
        for p in range(Nh):
            for q in range(Nh):
                v = self.bs.states[p][1:5]+self.bs.states[q][1:5]
                iv =  self.ident(v)
                combs_hh[p + q*Nh, :] = iv  #this one should not be zero, as most elements in array is already zero, or?
                #combs_ph[:,q + p*Np ] = iv
                #idents[p+q*Nh] = iv
        #idents = zeros((Nh**2))
        for p in range(Np):
            for q in range(Np):
                v = self.bs.states[p+ Nh][1:5]+self.bs.states[q+Nh][1:5]
                iv =  self.ident(v)
                combs_pp[:,p + q*Np] = iv  #this one should not be zero, as most elements in array is already zero, or?
                #combs_ph[:,q + p*Np ] = iv
                #idents[p+q*Nh] = iv                
        
        combs_hh[combs_hh!=combs_pp]=20000 #identify each pair of quantum numbers sharing the same added momentum
        t = where(combs_hh!=2000)
        i = self.bs.states[t[0]%Nh ].T
        j = self.bs.states[t[0]//Nh].T
        a = self.bs.states[t[1]%Np  + Nh].T
        b = self.bs.states[t[1]//Np + Nh].T

        data = self.bs.V(i,j,a,b)
        #print data[data!=0]
        self.Vhhpp = coo_matrix((data, (t[0], t[1])), shape=(self.Nh**2, self.Np**2)).tocsr()
        self.Vpphh = self.Vhhpp.T


def compare(Z1,Z2):
    Nx = len(Z1)
    Ny = len(Z1[0])
    EQ = True
    NE = 0
    toter = 0
    er = 0
    try:
        for i in range(Nx):
            for e in range(Ny):
                if Z1[i,e]!=Z2[i,e]:
                    #print Z1[i,e],Z2[i,e]
                    er = abs(Z1[i,e]-Z2[i,e])
                    if er>toter:
                        toter = er
                    NE = 1
                    
    except:
        print "NOT EQUAL, total failure"
        NE = 1
    return NE, toter

def Vpppp_check(Z1, bs):
    Np = bs.nstates-bs.nparticles
    Nh = bs.nparticles
    for a in range(Np):
        for b in range(Np):
            for c in range(Np):
                for d in range(Np):
                    if Z1[a + b*Np, c+ d*Np] != bs.v(a+Nh,b+Nh,c+Nh,d+Nh):
                        print a,b,c,d, Z1[a + b*Np, c+ d*Np], bs.v(a+Nh,b+Nh,c+Nh,d+Nh)

def Vhpph_check(Z1, bs):
    Np = bs.nstates-bs.nparticles
    Nh = bs.nparticles
    for i in range(Nh):
        for a in range(Np):
            for b in range(Np):
                for j in range(Nh):
                    if Z1[i + a*Nh, b+ j*Np] != bs.v(i,a+Nh,b+Nh,j):
                        print i,a,b,j, Z1[i + a*Nh, b+ j*Np], bs.v(i,a+Nh,b+Nh,j)
                        
def Vphhp_check(Z1, bs):
    Np = bs.nstates-bs.nparticles
    Nh = bs.nparticles
    for i in range(Nh):
        for a in range(Np):
            for b in range(Np):
                for j in range(Nh):
                    if Z1[a + i*Np, j+ b*Nh] != bs.v(a+Nh,i,j,b+Nh):
                        print i,a,b,j, Z1[a + i*Np, j+ b*Nh], bs.v(a+Nh,i,j,b+Nh)
def Vhhpp_check(Z1, bs):
    Np = bs.nstates-bs.nparticles
    Nh = bs.nparticles
    for i in range(Nh):
        for a in range(Np):
            for b in range(Np):
                for j in range(Nh):
                    if Z1[i + j*Nh, a+ b*Np] != bs.v(i,j,a+Nh,b+Nh):
                        print i,j,a,b, Z1[i + j*Nh, a+ b*Np], bs.v(i,j,a+Nh,b+Nh)
                                                           
t0 = clock()                        
tb = electronbasis(2,1.0,14)
t1 = clock()

 
print "Time spent on initializing basis:", t1-t0
print "====="
print "Number of states   :", tb.nstates
print "Number of particles:", tb.nparticles
print "====="
t0 = clock()
Q = CCD(tb)
t1 = clock()
print "Time spent on initializing solver:", t1-t0
#B = optimV(tb)
t0 = clock()
"""
B = blocks(tb)
t2 = clock()
B.setup_hhpp()
print "Time spent initializing vectorized interaction:", t2-t0
Vhhpp_check(B.Vhhpp.toarray(), tb)
"""
#Q.Vpppp = B.Vpppp


for i in range(20):
    Q.advance()



"""
print "pppp:", compare(Q.Vpppp.toarray(), B.Vpppp.toarray())
print compare(Q.Vhhpp.toarray(), B.Vhhpp)
print compare(Q.Vpphh.toarray(), B.Vhhpp.T)
print compare(Q.Vhpph.toarray(), B.Vhpph)
#Q.Vpppp = csr_matrix(B.Vpppp)

figure(1)
imshow(B.Vpppp.toarray())
show()

figure(2)
imshow(Q.Vpppp.toarray())
show()




"""