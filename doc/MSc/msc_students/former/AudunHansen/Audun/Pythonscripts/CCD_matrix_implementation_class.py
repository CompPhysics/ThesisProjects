from numpy import *
from time import *


class electronbasis():
    #A class for setting up a basis for the Homogeneous Electron Gas and interactions
    def __init__(self, N, rs, Nparticles):
        self.rs = rs
        self.states = []
        self.nstates = 0
        self.nparticles = Nparticles
        Nm = int(sqrt(N) + 1)
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
        
    def v(self,p,q,r,s):
        #Two body interaction
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

class CCD2():
    def __init__(self, bs):
        self.bs = bs
        self.nstates = bs.nstates            #total number of states
        self.Nh = bs.nparticles              #number of hole states (conflicting naming should be resolved in class electrongas)
        self.Np = self.nstates-bs.nparticles #number of particle states
        
        self.Vhhhh = zeros((self.Nh**2, self.Nh**2))
        self.Vhhpp = zeros((self.Nh**2, self.Np**2))
        self.Vphhp = zeros((self.Nh*self.Np, self.Nh*self.Np))
        self.Vhpph = zeros((self.Nh*self.Np, self.Nh*self.Np))
        self.Vpppp = zeros((self.Np**2, self.Np**2))
        self.Vpphh = zeros((self.Np**2, self.Nh**2))

        self.Tpphh = zeros((self.Np**2, self.Nh**2))
        self.Epphh = zeros((self.Np**2, self.Nh**2))

        self.setup_matrices_optimized()
        self.init_index_reorg()

        
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
        C = dot(self.Vhhpp,self.Tpphh)
        N = len(C)
        self.C_energy = .25*sum(C[range(0,N), range(0,N)])

    def test_matrices(self):
        t0 = clock()

        Nh = self.Nh
        Np = self.Np
        for i in range(Nh):
            for j in range(Nh):
                for a in range(Np):
                    for b in range(Np):
                        if self.Vhhpp[i + j*Nh, a + b*Np] != self.bs.v(i,j,a+Nh,b+Nh):
                            print "Vhhpp", i,j,a,b
                        
                        if self.Vphhp[a + i*Np, j + b*Nh] != self.bs.v(a+Nh,i,j,b+Nh):
                            print "Vphhp", a,i,j,b, self.Vphhp[a + i*Np, j + b*Nh], self.bs.v(a+Nh,i,j,b+Nh)
                        
                        if self.Vhpph[i + a*Nh, b + j*Np] != self.bs.v(i,a+Nh,b+Nh,j):
                            print "Vhpph", i,a,b,j
                        
                        if self.Vpphh[a + b*Np, i + j*Nh] != self.bs.v(a+Nh,b+Nh,i,j):
                            print "Vpphh", a,b,i,j

        t1 = clock()    
        for i in range(Nh):
            for j in range(Nh):
                for k in range(Nh):
                    for l in range(Nh):
                        val = self.bs.v(i,j,k,l)
                        if self.Vhhhh[i + j*Nh, k+ l*Nh] != val:
                            print "Vhhhh:", i,j,k,l,val,self.Vhhhh[i + j*Nh, k+ l*Nh]

        t2 = clock()
        for a in range(Np):
            for b in range(Np):
                for c in range(Np):
                    for d in range(Np):
                        val = self.bs.v(a+Nh,b+Nh,c+Nh,d+Nh)
                        if val!=self.Vpppp[a + b*Np, c+ d*Np]:
                            print "Vpppp:", a,b,c,d,val,self.Vpppp[a + b*Np, c+ d*Np]

                             
                        
        t3 = clock()
        print "matrices scanned in:", t1-t0,t2-t1,t3-t2
        print "                  :", t3-t0
    def setup_matrices_optimized(self):
        t0 = clock()

        Nh = self.Nh
        Np = self.Np
        for i in range(Nh):
            for j in range(i,Nh):
                for a in range(Np):
                    for b in range(a,Np):
                        #print self.bs.v(i,j,a+Nh,b+Nh), self.bs.v(a+Nh,i,j,b+Nh), self.bs.v(j,a+Nh,b+Nh,i), self.bs.v(a+Nh,b+Nh,i,j)

                        
                        #print self.bs.v(a+Nh,i,j,b+Nh) == self.bs.v(j,a+Nh,b+Nh,i)
                        
                        val = self.bs.v(a+Nh,i,j,b+Nh)
                        self.Vphhp[a + i*Np, j + b*Nh] = val
                        self.Vphhp[b + j*Np, i + a*Nh] = val 
                        
                        #self.Vhpph[i + a*Nh, b + j*Np] = val           
                        #self.Vhpph[j + b*Nh, a + i*Np] = val
                        
                        val = self.bs.v(j,a+Nh,b+Nh,i)
                        self.Vphhp[a + j*Np, i + b*Nh] = val
                        #self.Vhpph[j + a*Nh, b + i*Np] = val
                        self.Vphhp[b + i*Np, j + a*Nh] = val
                        #self.Vhpph[i + b*Nh, a + j*Np] = val
                        
                        
                        
                        
                        #val = self.bs.v(i,a+Nh,b+Nh,j)
                        #self.Vhpph[i + a*Nh, b + j*Np] = self.bs.v(i,a+Nh,b+Nh,j)
                        

                        
                        val = self.bs.v(a+Nh,b+Nh,i,j)
                        #val = self.bs.v(i,j,a+Nh,b+Nh)
                        #self.Vhhpp[i + j*Nh, a + b*Np] = val
                        #self.Vhhpp[j + i*Nh, b + a*Np] = val
                        #self.Vhhpp[j + i*Nh, a + b*Np] = -val
                        #self.Vhhpp[i + j*Nh, b + a*Np] = -val
                        
                        
                        #self.Vpphh[a + b*Np, i + j*Nh] = self.bs.v(a+Nh,b+Nh,i,j)
                        self.Vpphh[a + b*Np, i + j*Nh] = val
                        self.Vpphh[a + b*Np, j + i*Nh] = -val
                        self.Vpphh[b + a*Np, i + j*Nh] = -val
                        self.Vpphh[b + a*Np, j + i*Nh] = val                        
                        
                        
                        eps = self.bs.h(i,i) + self.bs.h(j,j) -self.bs.h(a+Nh,a+Nh) - self.bs.h(b+Nh,b+Nh)
                        
                        self.Tpphh[a + b*Np, i + j*Nh] = val/eps
                        self.Epphh[a + b*Np, i + j*Nh] = eps
                        
                        self.Tpphh[a + b*Np, j + i*Nh] = -val/eps
                        self.Epphh[a + b*Np, j + i*Nh] = eps
                        
                        self.Tpphh[b + a*Np, i + j*Nh] = -val/eps
                        self.Epphh[b + a*Np, i + j*Nh] = eps
                        
                        self.Tpphh[b + a*Np, j + i*Nh] = val/eps
                        self.Epphh[b + a*Np, j + i*Nh] = eps
        self.Vhhpp = self.Vpphh.T
        self.Vhpph = self.Vphhp.T
        t1 = clock()    
        for i in range(Nh):
            for j in range(i,Nh):
                for k in range(Nh):
                    for l in range(k,Nh):
                        val = self.bs.v(i,j,k,l)
                        self.Vhhhh[i + j*Nh, k+ l*Nh] = val
                        self.Vhhhh[j + i*Nh, l+ k*Nh] = val
                        self.Vhhhh[j + i*Nh, k+ l*Nh] = -val
                        self.Vhhhh[i + j*Nh, l+ k*Nh] = -val
                        
                        #self.Vhhhh[k + l*Nh, i+ j*Nh] = val
                        
                        
                        #self.Vhhhh[j + i*Nh, l+ k*Nh] = val
                        #self.Vhhhh[j + i*Nh, k+ l*Nh] = -val
                        #self.Vhhhh[i + j*Nh, l+ k*Nh] = -val
                        
                        
        t2 = clock()
        for a in range(Np):
            for b in range(a,Np):
                for c in range(Np):
                    for d in range(c,Np):
                        val = self.bs.v(a+Nh,b+Nh,c+Nh,d+Nh)
                        self.Vpppp[a + b*Np, c+ d*Np] = val 
                        self.Vpppp[b + a*Np, d+ c*Np] = val
                        self.Vpppp[b + a*Np, c+ d*Np] = -val
                        self.Vpppp[a + b*Np, d+ c*Np] = -val                

        t3 = clock()
        print "matrices set up in:", t1-t0,t2-t1,t3-t2
        print "                  :", t3-t0
    #setup functions for each diagram in the CCD amplitude equation (L = linear, Q=quadratic in amplitude)
    def dref(self, p,N):
        return p%N, p/N
    def test_L3(self):
        a,b,i,j = 23,10,12,9
        self.sL3()
        print a,b,i,j
        val = self.testfunction(a,b,i,j)
    def advance(self):
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
        self.PL3 = self.L3 - self.P_ab(self.L3) - self.P_ij(self.L3) + self.P_ij_ab(self.L3)
        self.PQ2 = self.Q2 - self.P_ij(self.Q2)
        self.PQ3 = self.Q3 - self.P_ij(self.Q3)
        self.PQ4 = self.Q4 - self.P_ab(self.Q4)
        
        #update amplitudes
        self.Tpphh = (self.Vpphh + .5*(self.L1 + self.L2) + self.PL3 + .25*self.Q1 + self.PQ2 - .5*(self.PQ3 + self.PQ4))/self.Epphh
        #calculate energy
        self.energy()
        print "        Correlation energy:", self.C_energy
        #print "                          :", self.e0_()  
    def testadvance(self):
        t0 = clock()
        
        #t1 = clock()
        self.sL1()
        self.sL2()
        self.sL3()
        self.sQ1()
        self.sQ2()
        self.sQ3()
        self.sQ4()
        
        
        t2 = clock()
        #self.PL3 = self.P_ij(self.P_ab(self.L3))
        self.PL3 = self.L3 - self.P_ab(self.L3) - self.P_ij(self.L3) + self.P_ij_ab(self.L3) #+#+ self.P_ij_ab(self.L3)
        self.PQ2 = self.Q2 - self.P_ij(self.Q2)
        self.PQ3 = self.Q3 - self.P_ij(self.Q3)
        self.PQ4 = self.Q4 - self.P_ab(self.Q4)
        t3 = clock()
        #self.LT = .5*(self.L1+self.L2) + self.PL3
        #self.QT = .5*(-self.PQ3-self.PQ4) + .25 * self.Q1 + self.PQ2 
        
        
        
        #a,b,i,j = 23,10,12,9
        for a in range(self.Np):
            for b in range(self.Np):
                for i in range(self.Nh):
                    for j in range(self.Nh):
                        l1,l2,l3,q1,q2,q3,q4,nnn = self.compare(a,b,i,j)
                        if nnn != True:
                            print "FOUND DISCREPANCY at a,b,i,j:", a,b,i,j
                            print "L1(vector):", self.L1[a + b*self.Np, i + j*self.Nh]-l1
                            print "L1(linear):", l1
                            print "L2(vector):", self.L2[a + b*self.Np, i + j*self.Nh]-l2
                            print "L2(linear):", l2
                            print "L3(vector):", self.PL3[a + b*self.Np, i + j*self.Nh]-l3
                            print "L3(linear):", l3
                            print "Q1(vector):", self.Q1[a + b*self.Np, i + j*self.Nh]-q1
                            print "Q1(linear):", q1
                            print "Q2(vector):", self.PQ2[a + b*self.Np, i + j*self.Nh]-q2
                            print "Q2(linear):", q2
                            print "Q3(vector):", self.PQ3[a + b*self.Np, i + j*self.Nh]-q3
                            print "Q3(linear):", q3
                            print "Q4(vector):", self.PQ4[a + b*self.Np, i + j*self.Nh]-q4
                            print "Q4(linear):", q4
    def compare(self, a,b,i,j):
        l1,l2,l3,q1,q2,q3,q4 = self.testfunction(a,b,i,j)
    
        n1 = self.L1[a + b*self.Np, i + j*self.Nh] == l1

        n2 = self.L2[a + b*self.Np, i + j*self.Nh] == l2

        n3 = self.PL3[a + b*self.Np, i + j*self.Nh]==l3

        n4 = self.Q1[a + b*self.Np, i + j*self.Nh]==q1

        n5 = self.PQ2[a + b*self.Np, i + j*self.Nh]==q2

        n6 = self.PQ3[a + b*self.Np, i + j*self.Nh]==q3

        n7 = self.PQ4[a + b*self.Np, i + j*self.Nh]==q4

        
        #print "Linear   :", self.LT[a + b*self.Np, i + j*self.Nh]
        #print "Quadratic:", self.QT[a + b*self.Np, i + j*self.Nh]
        #print "Full     :",(self.LT+self.QT)[a + b*self.Np, i + j*self.Nh]
       
        
        #self.Tpphh = self.LT +self.QT
        #self.Tpphh /= self.Epphh
        #remember to divide result by sp-energy eps^ab_ij
        #print "Time spent on first iteration:", clock()-t0
        #print "Time spent setting up contributions, ex. permutations:", t2-t0
        #print "Time spent on setting up contributions, incl. permutations:", t3-t0
        #print "Time spent on permuting V:", t1-t0
        #self.energy()
        #print "Correlation energy:", self.C_energy
        return l1,l2,l3,q1,q2,q3,q4,n1*n2*n3*n4*n5*n6*n7
    
    def sL1(self):
        #setup L1
        self.L1 = dot(self.Vpppp,self.Tpphh)
    def sL2(self):
        self.L2 = dot(self.Tpphh,self.Vhhhh)
    def sL3(self):
        #self.L3lin = self.TrL3(self.Vhpph, self.Tpphh)
        self.L3 = self.TL3()
    
    def sQ1(self):
        self.Q1 = dot(self.Tpphh,dot(self.Vhhpp, self.Tpphh))
    
    def sQ2(self):
        self.Q2 = self.TQ2(self.Tpphh, self.Vhhpp)
        
    def sQ3(self):
        self.Q3 = self.TQ3(self.Tpphh, self.Vhhpp)
    
    def sQ4(self):
        self.Q4 = self.TQ4(self.Tpphh, self.Vhhpp) #[a+b*Np, i + j*Nh] 
    
    #supporting functions
    def test_indexing(self):
        Np = self.Np
        Nh = self.Nh
        for i in range(Nh):
            for j in range(Nh):
                for a in range(Np):
                    for b in range(Np):
                        if self.test_reorg(a,b,i,j,self.Pba_ijA,self.ab_ij2ba_ij[0]):
                            print "1 - error"
                        if self.test_reorg(a,b,i,j,self.Pba_ijB,self.ab_ij2ba_ij[1]):
                            print "2 - error"
                        
                        if self.test_reorg(a,b,i,j,self.Pab_jiA,self.ab_ij2ab_ji[0]):
                            print "3 - error"
                        if self.test_reorg(a,b,i,j,self.Pab_jiB,self.ab_ij2ab_ji[1]):
                            print "4 - error"
                        
                        if self.Pai_bjA[a + i*Np, b + j*Np] != self.ab_ij2ai_bj[0][a + i*Np, b + j*Np]:
                            print "5 - error"
                        if self.Pai_bjB[a + i*Np, b + j*Np] != self.ab_ij2ai_bj[1][a + i*Np, b + j*Np]:
                            print "6 - error"
                        
                        if self.Pia_bjA[i + a*Nh, b + j*Np] != self.ab_ij2ia_bj[0][i + a*Nh, b + j*Np]:
                            print "7 - error"
                        if self.Pia_bjB[i + a*Nh, b + j*Np] != self.ab_ij2ia_bj[1][i + a*Nh, b + j*Np]:
                            print "8 - error"
                        
                        if self.Pai_jbA[a + i*Np, j + b*Nh] != self.ab_ij2ai_jb[0][a + i*Np, j + b*Nh]:
                            print "9 - error"
                        if self.Pai_jbB[a + i*Np, j + b*Nh] != self.ab_ij2ai_jb[1][a + i*Np, j + b*Nh]:
                            print "10 - error"
                            
                        if self.Pbj_aiA[b + j*Np, a + i*Np] != self.ab_ij2bj_ai[0][b + j*Np, a + i*Np]:
                            print "11 - error"
                        if self.Pbj_aiB[b + j*Np, a + i*Np] != self.ab_ij2bj_ai[1][b + j*Np, a + i*Np]:
                            print "12 - error"
                      
                        if self.Pi_jabA[i,j+ a*Nh + b *Np*Nh] != self.ab_ij2i_jab[0][i,j+ a*Nh + b *Np*Nh]:
                            print "13 - error"
                        if self.Pi_jabB[i,j+ a*Nh + b *Np*Nh] != self.ab_ij2i_jab[1][i,j+ a*Nh + b *Np*Nh]:
                            print "14 - error"
                        
                        if self.Pa_bjiA[a,b+ j*Np + i *Nh*Np] != self.ab_ij2a_bji[0][a,b+ j*Np + i *Nh*Np]:
                            print "15 - error"
                        if self.Pa_bjiB[a,b+ j*Np + i *Nh*Np] != self.ab_ij2a_bji[1][a,b+ j*Np + i *Nh*Np]:
                            print "16 - error"
                        
                        if self.Pb_aijA[b,a+ i*Np + j *Nh*Np] != self.ab_ij2b_aij[0][b,a+ i*Np + j *Nh*Np]:
                            print "17 - error"
                        if self.Pb_aijB[b,a+ i*Np + j *Nh*Np] != self.ab_ij2b_aij[1][b,a+ i*Np + j *Nh*Np]:
                            print "18 - error"
                        
                        if self.Q2_A[a + b*Np, i + j*Nh] != self.ai_bj2ab_ij[0][a + b*Np, i + j*Nh]:
                            print "23 - error"
                        if self.Q2_B[a + b*Np, i + j*Nh] != self.ai_bj2ab_ij[1][a + b*Np, i + j*Nh]:
                            print "24 - error"
                        
                        if self.Q3_A[a + b*Np, i + j*Nh] != self.i_jab2ab_ij[0][a + b*Np, i + j*Nh]:
                            print "19 - error"
                        if self.Q3_B[a + b*Np, i + j*Nh] != self.i_jab2ab_ij[1][a + b*Np, i + j*Nh]:
                            print "20 - error"
                        
                        if self.Q4_A[a + b*Np, i + j*Nh] != self.a_bij2ab_ij[0][a + b*Np, i + j*Nh]:
                            print "21 - error"
                        if self.Q4_B[a + b*Np, i + j*Nh] != self.a_bij2ab_ij[1][a + b*Np, i + j*Nh]:
                            print "22 - error"
                        
                        #self.ai_bj2ab_ij = [e1,e2]  

                            
    def test_reorg(self,a,b,i,j,P1,P2):
        return P1[a + b*self.Np, i + j*self.Nh] != P2[a + b*self.Np, i +j*self.Nh]
    def init_index_reorg(self):
        #Optimization madness. NB! Do not freak out! Or even better: do not read.
        Np = self.Np
        Nh = self.Nh
        self.ij2ij = outer(ones(Np**2, dtype = int), range(Nh**2))
        self.ab2ab = outer(range(Np**2), ones(Nh**2, dtype = int))
        
        #unaltered
        self.ab_ij2ab_ij = [outer(range(Np**2), ones(Nh**2, dtype = int)), outer(ones(Np**2, dtype = int), range(Nh**2))]     
        
        abba = arange(0,Np**3).reshape(Np**2,Np)[:,0]%(Np**2-1)
        abba[Np**2 - 1] = Np**2 - 1
        e1 = kron(ones((1,Nh**2), dtype = int).T, abba).T  
        e2 =  arange(0,Nh**3).reshape(Nh**2,Nh)[:,0]%(Nh**2- 1)
        e2[Nh**2 - 1] = Nh**2 - 1
        e2 =  kron(ones((1,Np**2), dtype = int).T, e2)
        self.ab_ij2ba_ji = [e1,e2]  
        
        ijji = arange(0,Nh**3).reshape(Nh**2,Nh)[:,0]%(Nh**2-1)
        ijji[Nh**2 - 1] = Nh**2 - 1
        self.ab_ij2ab_ji = [self.ab2ab, outer(ones(Np**2, dtype = int), ijji)]        
        
        abba = arange(0,Np**3).reshape(Np**2,Np)[:,0]%(Np**2-1)
        abba[Np**2 - 1] = Np**2 - 1
        self.ab_ij2ba_ij = [outer(abba, ones(Nh**2, dtype = int)),self.ij2ij]
        
        self.ab_ij2aj_bi =[kron(ones((Nh,Nh), dtype = int),arange(0,Np**2).reshape(Np,Np).T), kron(arange(0,Nh**2).reshape(Nh,Nh),ones((Nh,Nh), dtype = int))]
        self.ab_ij2ai_bj = [kron(ones((Nh,Nh), dtype = int), arange(0,Np**2).reshape(Np,Np).T), kron(arange(0,Nh**2).reshape(Nh,Nh).T,ones((Np,Np), dtype = int))]

        self.ab_ij2ia_bj = [kron(kron(ones((Nh), dtype =int).T, arange(0,Np**2).reshape(Np,Np).T), ones((Nh,1), dtype = int)), kron(ones((Np,1), dtype = int),kron(arange(0,Nh**2).reshape(Nh,Nh).T, ones((Np),dtype = int))) ]
        
        self.ab_ij2ia_jb = [kron(kron(arange(0,Np**2).reshape(Np,Np).T, ones((Nh),dtype = int)),ones((Nh,1), dtype = int)),kron(ones((Np,1), dtype = int),kron(ones((Np), dtype =int).T, arange(0,Nh**2).reshape(Nh,Nh).T))]

        self.ab_ij2ai_jb = [kron(ones((Nh,1), dtype = int),kron(arange(0,Np**2).reshape(Np,Np).T, ones((Nh),dtype = int))), kron(ones((Nh**2,1), dtype = int).T,kron(arange(0,Nh**2).reshape(Nh,Nh).T, ones((1,Np),dtype = int).T))]

        
        self.ab_ij2bj_ai = [kron(ones((Nh,Nh), dtype = int), arange(0,Np**2).reshape(Np,Np).T).T,kron(arange(0,Nh**2).reshape(Nh,Nh).T,ones((Np,Np), dtype = int)).T]
        
        self.ab_ij2i_jab = [kron(arange(0,Np**2),ones((Nh,Nh), dtype =int)),kron(ones((1,Np**2), dtype = int), arange(0,Nh**2).reshape(Nh,Nh).T) ]        
        
        e1 =  kron(ones((1,Nh**2), dtype = int), arange(0,Np**2).reshape(Np,Np).T)
        e2 = kron(arange(0,Nh**2),ones((Np,Np), dtype =int))
        self.ab_ij2a_bji = [e1,e2]    
        
        e1 =  kron(ones((1,Nh**2), dtype = int), arange(0,Np**2).reshape(Np,Np))
        e2 = arange(0,Nh**3).reshape(Nh**2,Nh)[:,0]%(Nh**2-1)
        e2[Nh**2 - 1] = Nh**2 - 1 
        e2 = arange(0,Nh**2)
        e2 = kron(e2,ones((Np,Np), dtype =int))
        self.ab_ij2b_aij = [e1,e2]  
        
        e1 =  kron(ones((Np**2,Nh), dtype = int), arange(0,Nh).reshape(Nh))
        e2 = arange(0,Nh*Np**2).reshape(Np**2,Nh) 
        e2 = kron(e2,ones((1,Nh), dtype =int))
        self.i_jab2ab_ij = [e1,e2]        

        
        e1 =  kron(ones((Nh**2, Np), dtype = int), arange(0,Np).reshape(Np)).T
        e2 = arange(0,Np*Nh**2).reshape(Nh**2,Np)
        e2 = kron(e2,ones((1,Np), dtype =int)).T
        self.a_bij2ab_ij = [e1,e2]

        e1 = kron(ones((Np, Nh), dtype = int), arange(Np*Nh).reshape(Nh,Np).T)
        e2 = kron(arange(0,Np*Nh).reshape(Nh,Np).T,ones((Np,Nh),dtype = int))
        self.ai_bj2ab_ij = [e1,e2]     
        
        e1 = kron( arange(Np**2), ones((1, Nh**2), dtype = int).T).T
        e2 = arange(0,Nh**3).reshape(Nh**2,Nh)[:,0]%(Nh**2-1)
        e2[Nh**2 - 1] = Nh**2 - 1 
        e2 = outer(ones((Np**2), dtype = int), e2)
        self.ab_ji2ab_ij = [e1,e2]
        
        #V-permutations
        e1 =arange(0,Np*Nh).reshape(Np,Nh).T
        e1 = kron(e1,ones((1,Np), dtype = int).T)
        e1 = kron(ones((Nh), dtype = int),e1)
        e2 = arange(0,Np*Nh).reshape(Nh,Np).T
        e2 = kron(ones((1, Nh), dtype = int),kron(e2,ones((1,Np), dtype = int)).T).T
        self.vl3 = [e1,e2]
        self.VL3 = self.Vhpph[self.vl3]
        
        self.vq2 = [kron(arange(0,Nh**2).reshape(Nh,Nh).T,ones((Np,Np), dtype = int)),kron(ones((Nh,Nh), dtype = int), arange(0,Np**2).reshape(Np,Np).T)]
        self.VQ2 = self.Vhhpp[self.vq2]
        
        e1 =  kron(ones((Np*Np), dtype = int), arange(0,Nh**2).reshape(Nh,Nh)).T
        e2 = arange(0,Np**3).reshape(Np**2,Np)[:,0]%(Np**2-1)
        e2[Np**2 - 1] = Np**2 - 1 
        e2 = kron(kron(ones((Nh,1), dtype = int),e2),ones((Nh),dtype =int)).T
        self.vq3 = [e1,e2]
        self.VQ3 = self.Vhhpp[self.vq3]
        
        e1 = arange(0,Nh**3).reshape(Nh**2,Nh)[:,0]%(Nh**2-1)
        e1[Nh**2 - 1] = Nh**2 - 1 
        e1 =  kron(e1,ones((Np,Np), dtype = int)).T
        e2 =  kron(ones((1,Nh**2), dtype = int).T, arange(0,Np**2).reshape(Np,Np).T)
        self.vq4 = [e1,e2]  
        self.VQ4 = self.Vhhpp[self.vq4] 
           
    def init_index_reorganization(self):
        #initialize all index reorganization matrices
        Np = self.Np
        Nh = self.Nh
        """
        self.Pai_bjA = zeros((Np*Nh,Nh*Np), dtype = int)
        self.Pai_bjB = zeros((Np*Nh,Nh*Np), dtype = int)
        
        self.Pba_ijA = zeros((Np**2,Nh**2), dtype = int)
        self.Pba_ijB = zeros((Np**2,Nh**2), dtype = int)
        
        self.Pba_jiA = zeros((Np**2,Nh**2), dtype = int)
        self.Pba_jiB = zeros((Np**2,Nh**2), dtype = int)
        
        self.Pab_jiA = zeros((Np**2,Nh**2), dtype = int)
        self.Pab_jiB = zeros((Np**2,Nh**2), dtype = int)
       
        
        self.Pai_jbA = zeros((Np*Nh,Nh*Np), dtype = int)
        self.Pai_jbB = zeros((Np*Nh,Nh*Np), dtype = int)
        
        self.Pia_bjA = zeros((Np*Nh,Nh*Np), dtype = int)
        self.Pia_bjB = zeros((Np*Nh,Nh*Np), dtype = int)
        
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
        """
        
        self.VL3 = zeros((Nh*Np, Nh*Np))
        self.VQ2  = zeros((self.Np*self.Nh, self.Np*self.Nh))
        self.VQ3  = zeros((self.Nh*self.Np**2, self.Nh))
        self.VQ4  = zeros((self.Np*self.Nh**2, self.Np))
        
        for i in range(Nh):
            for j in range(Nh):
                for a in range(Np):
                    for b in range(Np):
                        self.VL3[a + i*Np, b + j*Np]      = self.Vhpph[i + b*Nh, a + j*Np]
                        self.VQ2[a + i*Np, b + j*Np]      = self.Vhhpp[i + j*Nh, a + b*Np] 
                        self.VQ3[i + a*Nh + b *Np*Nh, j]  = self.Vhhpp[i+ j*Nh, b + a*Np]
                        self.VQ4[b + j*Np + i *Nh*Np, a]  = self.Vhhpp[i+ j*Nh, b + a*Np]

                        """
                        self.Pba_ijA[b + a*Np, i + j*Nh] = a + b*Np
                        self.Pba_ijB[b + a*Np, i + j*Nh] = i + j*Nh
                        
                        self.Pba_jiA[b + a*Np, j + i*Nh] = a + b*Np
                        self.Pba_jiB[b + a*Np, j + i*Nh] = i + j*Nh
                        
                        self.Pab_jiA[a + b*Np, j + i*Nh] = a + b*Np
                        self.Pab_jiB[a + b*Np, j + i*Nh] = i + j*Nh
                        
                        self.Pba_jiA[b + a*Np, j + i*Nh] = a + b*Np
                        self.Pba_jiB[b + a*Np, j + i*Nh] = i + j*Nh
                        
                        self.Pai_bjA[a + i*Np, b + j*Np] = a + b*Np
                        self.Pai_bjB[a + i*Np, b + j*Np] = i + j*Nh

                        self.Pia_bjA[i + a*Nh, b + j*Np] = a + b*Np
                        self.Pia_bjB[i + a*Nh, b + j*Np] = i + j*Nh
                        
                        self.Pai_jbA[a + i*Np, j + b*Nh] = a + b*Np
                        self.Pai_jbB[a + i*Np, j + b*Nh] = i + j*Nh
                        #ptrA[a_ + i_*Np, j_ + b_*Nh] = a_ + b_*Np
                        #ptrB[a_ + i_*Np, j_ + b_*Nh] = i_ + j_*Nh
                        
                        #L3A[a_ + b_*Np, i_ + j_*Nh] = a_ + i_*Np
                        #L3B[a_ + b_*Np, i_ + j_*Nh] = b_ + j_*Np
                        self.L3_A[a + b*Np, i + j*Nh] = a + i*Np
                        self.L3_B[a + b*Np, i + j*Nh] = b + j*Np
                        
                        #self.L3_A[a + b*Np, j + i*Nh] = a + j*Np
                        #self.L3_B[a + b*Np, j + i*Nh] = b + i*Np
                        #L3[a + b*Np, i + j*Nh] = L3_[a + i*Np, b + j*Np]
                        
                        self.Pbj_aiA[b + j*Np, a + i*Np] = a + b*Np
                        self.Pbj_aiB[b + j*Np, a + i*Np] = i + j*Nh
                        
                        self.Pi_jabA[i,j+ a*Nh + b *Np*Nh] = a + b*Np #TQ31&TQ32
                        self.Pi_jabB[i,j+ a*Nh + b *Np*Nh] = i + j*Nh
                        
                        self.Pa_bjiA[a,b+ j*Np + i *Nh*Np] = a + b*Np
                        self.Pa_bjiB[a,b+ j*Np + i *Nh*Np] = j + i*Nh ############# Check this one
                        
                        self.Pb_aijA[b, a + i*Np + j *Nh*Np] = a + b*Np
                        self.Pb_aijB[b, a + i*Np + j *Nh*Np] = i + j*Nh
                        
                        self.Q2_A[a + b*Np, i + j*Nh] = a + i*Np
                        self.Q2_B[a + b*Np, i + j*Nh] = b + j*Np
                        
                        self.Q3_A[a + b*Np, i + j*Nh] = i
                        self.Q3_B[a + b*Np, i + j*Nh] = j + a*Nh + b*Nh*Np
    
                        self.Q4_A[a + b*Np, i + j*Nh] = a
                        self.Q4_B[a + b*Np, i + j*Nh] = b + i*Np + j*Np*Nh
                        """
                        
    def TrV(self):
        Nh = self.Nh
        Np = self.Np
        self.VL3 = zeros((Nh*Np, Nh*Np))
        #self.VL3 = self.Vhpph[self.Pai_bjA, self.Pai_bjB]
        #self.VQ2 = self.Vhhpp[]
        self.VQ2  = zeros((self.Np*self.Nh, self.Np*self.Nh))
        self.VQ3  = zeros((self.Nh*self.Np**2, self.Nh))
        self.VQ4  = zeros((self.Np*self.Nh**2, self.Np))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        self.VL3[a + i*Np, b + j*Np] = self.Vhpph[i + b*Nh, a + j*Np]
                        self.VQ2[a + i*Np, b + j*Np]  = self.Vhhpp[i + j*Nh, a + b*Np] 
                        self.VQ3[i + a*Nh + b *Np*Nh, j]  = self.Vhhpp[i+ j*Nh, b + a*Np]
                        self.VQ4[b + j*Np + i *Nh*Np, a]  = self.Vhhpp[i+ j*Nh, b + a*Np]
                        
                        #    VL3[a + i*Np, b + j*Np] =          V[i + b*Nh, a + j*Np]
                        #self.VL3[b + j*Np, a + i*Np] = self.Vhpph[i + a*Nh, b + j*Nh]
                        
                        #VL3[a + i*Np, b + j*Np] = V[i + b*Nh, a + j*Np]
    def TL3(self):
        Nh = self.Nh
        Np = self.Np
        #TL3 = zeros((Nh*Np, Nh*Np))
        
        #TL3[a + i*Np, b + j*Np] = T[a + b*Np, i + j*Nh]
        #VL3[a + i*Np, b + j*Np] = V[i + b*Nh, a + j*Np]
        #L3[a + b*Np, i + j*Nh] = L3_[a + i*Np, b + j*Np]
        
        #TL3 = self.Tpphh[self.Pai_bjA, self.Pai_bjB]
        
        TL3 = self.Tpphh[self.ab_ij2ai_bj]
        L3_ = dot(TL3, self.VL3)
        #return L3 = L3_[self.L3_A, self.L3_B]
        return L3_[self.ai_bj2ab_ij]

    def TrL3(self,V,T):
        #Vhpph -> Vphph
        Nh = self.Nh
        Np = self.Np
        TL3 = zeros((Nh*Np, Nh*Np))

        VL3 = zeros((Nh*Np, Nh*Np))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        TL3[a + i*Np, b + j*Np] = T[a + b*Np, i + j*Nh]

                        VL3[a + i*Np, b + j*Np] = V[i + b*Nh, a + j*Np]
                        #self.VL3[a + i*Np, b + j*Np] = V[i + b*Nh, a + j*Np]
        #print TL3
        #print VL3
        #print self.VL3
        L3_ = dot(TL3,VL3)
        L3 = zeros((Np**2, Nh**2))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        L3[a + b*Np, i + j*Nh] = L3_[a + i*Np, b + j*Np]
        return L3            
                              
    def TQ2(self,T,V):
        #support for Q2 calculation
        Nh = self.Nh
        Np = self.Np

        #TQ21 = self.Tpphh[self.Pai_bjA, self.Pai_bjB]
        #TQ22 = self.Tpphh[self.Pbj_aiA, self.Pbj_aiB]
        TQ21 = self.Tpphh[self.ab_ij2ai_bj]
        TQ22 = self.Tpphh[self.ab_ij2bj_ai]
        Q2_ = dot(TQ21,dot(self.VQ2, TQ22))
        return Q2_[self.ai_bj2ab_ij]
        
        #return Q2_[self.Q2_A, self.Q2_B]

    def TQ3(self,T,V):
        Nh = self.Nh
        Np = self.Np
        #TQ31 = T[self.Pi_jabA, self.Pi_jabB]
        TQ31 = T[self.ab_ij2i_jab]
        Q3_ = dot(TQ31, dot(self.VQ3, TQ31))
        
        return Q3_[self.i_jab2ab_ij]

    def TQ4(self,T,V):
        TQ41 = T[self.ab_ij2a_bji]
        Q4_ = dot(TQ41,dot(self.VQ4, TQ41))
        return Q4_[self.a_bij2ab_ij]


          
    def TrQ2(self,T,V):
        #support for Q2 calculation
        Nh = self.Nh
        Np = self.Np
        TQ21 = zeros((self.Np*self.Nh, self.Np*self.Nh))
        #VQ2  = zeros((self.Np*self.Nh, self.Np*self.Nh))
        TQ22 = zeros((self.Np*self.Nh, self.Np*self.Nh))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        TQ21[a + i*Np, b + j*Np] = T[a + b*Np, i + j*Nh]
                        TQ22[b + j*Np, a + i*Np] = T[a + b*Np, i + j*Nh]
                        #VQ2[a + i*Np, b + j*Np]  = V[i + j*Nh, a + b*Np]
        Q2_ = dot(TQ21,dot(self.VQ2, TQ22))
        Q2  = zeros((Np**2, Nh**2))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        Q2[a + b*Np, i + j*Nh] = Q2_[a + i*Np, b + j*Np]
        
        #return TQ21, VQ2, TQ22
        return Q2  
    def TrQ3(self,T,V):
        Nh = self.Nh
        Np = self.Np
        TQ31 = zeros((Nh, Nh*Np**2))
        TQ32 = zeros((Nh, Nh*Np**2))
        #VQ3  = zeros((Nh*Np**2, Nh))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        TQ31[i,j+ a*Nh + b *Np*Nh] = T[a + b*Np, i + j*Nh]
                        TQ32[i, j + a*Nh + b *Np*Nh] = T[a + b*Np, i + j*Nh]
                        #VQ3[i + a*Nh + b *Np*Nh, j]  = V[i+ j*Nh, b + a*Np]
        Q3_ = dot(TQ31,dot(self.VQ3, TQ32))
        Q3  = zeros((Np**2, Nh**2))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        Q3[a + b*Np, i + j*Nh] = Q3_[i, j + a*Nh + b*Nh*Np]    
        return Q3 
    def TrQ4(self,T,V):
        Nh = self.Nh
        Np = self.Np
        TQ41 = zeros((Np, Np*Nh**2))
        TQ42 = zeros((Np, Np*Nh**2))
        #VQ4  = zeros((Np*Nh**2, Np))     
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        TQ41[a,b+ j*Np + i *Nh*Np] = T[a + b*Np, j + i*Nh]
                        TQ42[b, a + i*Np + j *Nh*Np] = T[b + a*Np, i + j*Nh]
                        #VQ4[b + j*Np + i *Nh*Np, a]  = V[i+ j*Nh, b + a*Np]
        Q4_ = dot(TQ41,dot(self.VQ4, TQ42))
        Q4  = zeros((Np**2, Nh**2))
        for i in range(Nh):
            for j in range(Nh): 
                for a in range(Np):
                    for b in range(Np):
                        Q4[a + b*Np, i + j*Nh] = Q4_[a, b + i*Np + j*Np*Nh]    
        return Q4 
        
    #Permuting functions       
    def P_ab(self,M):
        return M[self.ab_ij2ba_ij]

    def P_ij(self,M):
        return M[self.ab_ij2ab_ji]

    def P_ij_ab(self, M):
        return M[self.ab_ij2ba_ji]
    
    def testfunction(self,a,b,i,j):
        #test that every contribution is reproduced exactly
        Np = self.Np
        Nh = self.Nh
        self.L1_t = 0.0
        self.L2_t = 0.0
        self.L3_t = 0.0
        
        self.Q1_t = 0.0
        self.Q2_t = 0.0
        self.Q3_t = 0.0
        self.Q4_t = 0.0
        
        for c in range(Np):
            for d in range(Np):
                self.L1_t += self.Vpppp[a + b*Np, c+ d*Np]*self.Tpphh[c + d*Np, i + j*Nh]
        #self.L1_t*=.5
        for k in range(Nh):
            for l in range(Nh):
                self.L2_t += self.Vhhhh[k + l*Nh, i + j*Nh]*self.Tpphh[a + b*Np, k + l*Nh]
        #self.L2_t*=.5
        l1,l2,l3,l4 = 0.0,0.0,0.0,0.0
        for k in range(Nh):
            for c in range(Np):
                self.L3_t += self.Vhpph[k + b*Nh, c + j*Np]*self.Tpphh[a + c*Np, i + k*Nh]
                self.L3_t -= self.Vhpph[k + b*Nh, c + i*Np]*self.Tpphh[a + c*Np, j + k*Nh]
                self.L3_t -= self.Vhpph[k + a*Nh, c + j*Np]*self.Tpphh[b + c*Np, i + k*Nh]
                self.L3_t += self.Vhpph[k + a*Nh, c + i*Np]*self.Tpphh[b + c*Np, j + k*Nh]
                l1 += self.Vhpph[k + b*Nh, c + j*Np]*self.Tpphh[a + c*Np, i + k*Nh]
                l2 -= self.Vhpph[k + b*Nh, c + i*Np]*self.Tpphh[a + c*Np, j + k*Nh]
                l3 -= self.Vhpph[k + a*Nh, c + j*Np]*self.Tpphh[b + c*Np, i + k*Nh]
                l4 += self.Vhpph[k + a*Nh, c + i*Np]*self.Tpphh[b + c*Np, j + k*Nh]
        #print l1,l2,l3,l4
        for k in range(Nh):
            for l in range(Nh):
                for c in range(Np):
                    for d in range(Np):
                        self.Q1_t += self.Vhhpp[k + l*Nh, c + d*Np]*self.Tpphh[c + d*Np, i + j*Nh]*self.Tpphh[a + b*Np, k + l*Nh]
                        
                        self.Q2_t += self.Vhhpp[k + l*Nh, c + d*Np]*self.Tpphh[a + c*Np, i + k*Nh]*self.Tpphh[b + d*Np, j + l*Nh]
                        self.Q2_t -= self.Vhhpp[k + l*Nh, c + d*Np]*self.Tpphh[a + c*Np, j + k*Nh]*self.Tpphh[b + d*Np, i + l*Nh]
                        
                        self.Q3_t += self.Vhhpp[k + l*Nh, c + d*Np]*self.Tpphh[d + c*Np, i + k*Nh]*self.Tpphh[a + b*Np, l + j*Nh]
                        self.Q3_t -= self.Vhhpp[k + l*Nh, c + d*Np]*self.Tpphh[d + c*Np, j + k*Nh]*self.Tpphh[a + b*Np, l + i*Nh]
                        
                        self.Q4_t += self.Vhhpp[k + l*Nh, c + d*Np]*self.Tpphh[a + c*Np, l + k*Nh]*self.Tpphh[d + b*Np, i + j*Nh]
                        self.Q4_t -= self.Vhhpp[k + l*Nh, c + d*Np]*self.Tpphh[b + c*Np, l + k*Nh]*self.Tpphh[d + a*Np, i + j*Nh]
        
        L1 = self.L1_t + self.L2_t + self.L3_t
        
        Q1 = self.Q1_t + self.Q2_t + self.Q3_t + self.Q4_t
        return self.L1_t, self.L2_t, self.L3_t, self.Q1_t, self.Q2_t, self.Q3_t, self.Q4_t
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
#Q.test_matrices()

print "Time spent initializing solver:", t1-t0
t0 = clock()
Q.advance()
t1 = clock()
print "Time spent on first iteration:", t1-t0
#Q.energy()
#print "        Correlation energy(0):", Q.C_energy
#print "                             :", Q.e0_()  
#for i in range(30):
#    t0 = clock()
#    Q.advance()
#    print clock()-t0
#Q.testadvance()
#from matplotlib.pyplot import *
#imshow(Q.Tpphh, cmap = "RdGy")
#show()
