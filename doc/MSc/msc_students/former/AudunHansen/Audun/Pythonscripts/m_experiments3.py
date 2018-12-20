from numpy import *
from scipy.sparse import coo_matrix, csr_matrix
from matplotlib.pyplot import *
from time import *

class consistency_test():
    def __init__(self):
        pass
    def test_Vhhhh(self,M, bs):
        #Old linear approach
        Np = bs.nstates-bs.nparticles
        Nh = bs.nparticles
        ret = True
        n = 0
        for i in range(Nh):
            for j in range(Nh):
                for k in range(Nh):
                    for l in range(Nh):
                        if M[i + j*Nh, k+l*Nh] != bs.v(i , j ,k, l):
                            print M[i + j*Nh, k+l*Nh], bs.v(i , j ,k, l)
                            n += 1
                            ret = False
        if ret == False:
            print "Found %i errors in hhhh" % n
        else:
            print "*************************"
            print "* Vhhhh is consistent!  *"
            print "*************************"
    
    def test_Vpppp(self,M, bs):
        Np = bs.nstates-bs.nparticles
        Nh = bs.nparticles
        n = 0
        ret = True
        for a in range(Np):
            for b in range(Np):
                for c in range(Np):
                    for d in range(Np):
                        if M[a + b*Np, c+d*Np] != bs.v(a + Nh, b + Nh,c + Nh, d + Nh):
                            #print M[a + b*Np, i+j*Nh], bs.v(a + Nh, b + Nh,i, j)
                            n += 1
                            ret = False
        if ret == False:
            print "Found %i errors in pppp" % n
        else:
            print "*************************"
            print "* Vpppp is consistent!  *"
            print "*************************"                    
                        
    def test_Vhpph(self,M, bs):
        #Old linear approach
        Np = bs.nstates-bs.nparticles
        Nh = bs.nparticles
        ret = True
        n = 0
        for a in range(Np):
            for b in range(Np):
                for i in range(Nh):
                    for j in range(Nh):
                        if M[i + a*Nh, b+j*Np] != bs.v(i,a+Nh,b+Nh, j):
                            #print M[a + b*Np, i+j*Nh], bs.v(a + Nh, b + Nh,i, j)
                            n += 1
                            ret = False
        if ret == False:
            print "Found %i errors in hpph" % n
        else:
            print "*************************"
            print "* Vhpph is consistent!  *"
            print "*************************"
    
    def test_Vpphh(self,M, bs):
        #Old linear approach
        Np = bs.nstates-bs.nparticles
        Nh = bs.nparticles
        ret = True
        n = 0
        for a in range(Np):
            for b in range(Np):
                for i in range(Nh):
                    for j in range(Nh):
                        if M[a + b*Np, i+j*Nh] != bs.v(a + Nh, b + Nh,i, j):
                            #print M[a + b*Np, i+j*Nh], bs.v(a + Nh, b + Nh,i, j)
                            n += 1
                            ret = False
        if ret == False:
            print "Found %i errors in pphh" % n
        else:
            print "*************************"
            print "* Vpphh is consistent!  *"
            print "*************************"
    
    def test_Vhhpp(self,M, bs):
        #Old linear approach
        Np = bs.nstates-bs.nparticles
        Nh = bs.nparticles
        ret = True
        n = 0
        for a in range(Np):
            for b in range(Np):
                for i in range(Nh):
                    for j in range(Nh):
                        if M[i+j*Nh,a + b*Np] != bs.v(i, j,a + Nh, b + Nh):
                            #print M[a + b*Np, i+j*Nh], bs.v(a + Nh, b + Nh,i, j)
                            n += 1
                            ret = False
        if ret == False:
            print "Found %i errors in hhpp" % n
        else:
            print "*************************"
            print "* Vhhpp is consistent!  *"
            print "*************************"

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

class flexitude():
    #T2 amplitudes, a flexible array for t2 storage (fast reorganization, still very high memory consumption)
    def __init__(self, Np, Nh):
        self.Np = Np
        self.Nh = Nh
        self.M = zeros((Np**2, Nh**2))
        self.init_index_reorg()
    def ab_ij(self):
        return self.M
    def ij_ab(self):
        return self.M.T
                
    def ba_ij(self):
        try:
            return self.Mp2p1_h1h2
        except:
            self.Mp2p1_h1h2 = self.M[self.ab_ij2ba_ij]
            return self.Mp2p1_h1h2
    def ab_ji(self):
        try:
            return self.Mp1p2_h2h1
        except:
            self.Mp1p2_h1h2 = self.M[self.ab_ij2ab_ji]
            return self.Mp1p2_h2h1
    def ba_ji(self):
        try:
            return self.Mp2p1_h2h1
        except:
            self.Mp2p1_h2h1 = self.M[self.ab_ij2ba_ji]
            return self.Mp2p1_h2h1
    def aj_bi(self):
        try:
            return self.Mp1h2_p2h1
        except:
            self.Mp1h2_p2h1 = self.M[self.ab_ij2aj_bi]
            return self.Mp1h2_p2h1
    
    def ai_bj(self):
        try:
            return self.Mp1h1_p2h2
        except:
            self.Mp1h1_p2h2 = self.M[self.ab_ij2ai_bj]
            return self.Mp1h1_p2h2

    def ia_bj(self):
        try:
            return self.Mh1p1_p2h2
        except:
            self.Mh1p1_p2h2 = self.M[self.ab_ij2ia_bj]
            return self.Mh1p1_p2h2
            
    def ia_jb(self):
        try:
            return self.Mh1p1_h2p2
        except:
            self.Mh1p1_h2p2 = self.M[self.ab_ij2ia_jb]
            return self.Mh1p1_h2p2
            
    def ai_jb(self):
        try:
            return self.Mp1h1_h2p2
        except:
            self.Mp1h1_h2p2 = self.M[self.ab_ij2ai_jb]
            return self.Mp1h1_h2p2
            
    def bj_ai(self):
        try:
            return self.Mp2h2_p1h1
        except:
            self.Mp2h2_p1h1 = self.M[self.ab_ij2bj_ai]
            return self.Mp2h2_p1h1
                       
    def i_jab(self):
        try:
            return self.Mh1_h2p1p2
        except:
            self.Mh1_h2p1p2 = self.M[self.ab_ij2i_jab]
            return self.Mh1_h2p1p2
                 
    def a_bji(self):
        try:
            return self.Mp1_p2h2h1
        except:
            self.Mp1_p2h2h1 = self.M[self.ab_ij2a_bji]
            return self.Mp1_p2h2h1     
 
    def b_aij(self):
        try:
            return self.Mp2_p1h1h2
        except:
            self.Mp2_p1h1h2 = self.M[self.ab_ij2b_aij]
            return self.Mp2_p1h1h2 
             
    #def b_aij(self):
    #    try:
    #        return self.Mp2_p1h1h2
    #    except:
    #        self.Mp2_p1h1h2 = self.M[self.ab_ij2b_aij]
    #        return self.Mp2_p1h1h2  
    """

        self.i_jab2ab_ij = [e1,e2]        

        self.a_bij2ab_ij = [e1,e2]

        self.ai_bj2ab_ij = [e1,e2]     

        self.ab_ji2ab_ij = [e1,e2]    
    """ 
    def init_index_reorg(self):
        #Reorganization madness using fancy indexing. NB! Do not freak out! Or even better: do not read.
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
        self.ab_ij2ba_ji = [e1,e2]   #
        
        ijji = arange(0,Nh**3).reshape(Nh**2,Nh)[:,0]%(Nh**2-1)
        ijji[Nh**2 - 1] = Nh**2 - 1
        self.ab_ij2ab_ji = [self.ab2ab, outer(ones(Np**2, dtype = int), ijji)]# 
        
        abba = arange(0,Np**3).reshape(Np**2,Np)[:,0]%(Np**2-1)
        abba[Np**2 - 1] = Np**2 - 1
        self.ab_ij2ba_ij = [outer(abba, ones(Nh**2, dtype = int)),self.ij2ij]#
        
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
        """
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
        """
        
        
        
class ashape():
    #class for quick reorganizing of array data in a+b*Np, i+j*Nh squences
    #For use with compressed row sparse matrices
    def __init__(self, M, Np, Nh, t0, t1):
        self.t0 = t0 #these are the coordinates of the elements the complete 4d tensor
        self.t1 = t1
        
        self.a = kron(ones(Nh**2), arange(Np**2)%Np) #preorganized 
        self.b = kron(ones(Nh**2), arange(Np**2)//Np)
        self.i = kron(arange(Nh**2)%Nh,ones(Np**2))
        self.j = kron(arange(Nh**2)//Nh,ones(Np**2))
        self.M = ravel(M.T)
    #The following functions is different representations of the same tensor, for use in matrix multiplications
    
    #########
    ##
    ##  PP-HH
    ##  
    #########
    def ab_ij(self):
        try:
            return self.Vab_ij
        except:
            self.Vab_ij = coo_matrix((data, (a + b*self.Np, i + j*self.Nh)), shape=(self.Np*self.Np, self.Nh*self.Nh)).tocsr()
            return self.Vaa_ij
    #########
    ##
    ##  PH-PH
    ##
    #########
         
    def ai_bj(self):
        try:
            return self.Vai_bj
        except:
            self.Vai_bj = coo_matrix((data, (a + i*self.Np, b + j*self.Np)), shape=(self.Np*self.Nh, self.Np*self.Nh)).tocsr()
            return self.Vai_bj
    def bi_aj(self):
        try:
            return self.Vbi_aj
        except:
            self.Vai_bj = coo_matrix((data, (b + i*self.Np, a + j*self.Np)), shape=(self.Np*self.Nh, self.Np*self.Nh)).tocsr()
            return self.Vbi_aj
    def aj_bi(self):
        try:
            return self.Vaj_bi
        except:
            self.Vaj_bi = coo_matrix((data, (a + j*self.Np, b + i*self.Np)), shape=(self.Np*self.Nh, self.Np*self.Nh)).tocsr()
            return self.Vaj_bi
    def a_bij(self):
        try:
            return self.Va_bij
        except:
            self.Va_bij = coo_matrix((data, (a, b + i*self.Np + j*self.Np*self.Nh)), shape=(self.Np,self.Nh*self.Np*self.Nh)).tocsr()
            return self.Va_bij 

class optiV():
    #quick initialization of interactions and amplitudes
    def __init__(self, bs):
        self.bs = bs
        self.Nh = bs.nparticles
        self.Np = bs.nstates-self.Nh
        self.Nm = bs.Nm
        self.Nm2 = self.Nm**2
        #self.sVpppp()
        #self.sVhhhh()
        #self.sVpphh()
        #self.sVphhp()
        #self.setup_hpph()
    def kdplt(self,oneN, states, P, Q, N):
        #kdplus, elementwise
        #oneNp = ones((Np**2, 1), dtype = int)
        return kron(oneN, states[:,N][P] + states[:,N][Q])
    
    def kdplq(self,states, P, Q, N):
        #kdplus, elementwise
        return states[:,N][P] + states[:,N][Q]
        
    def sVpppp(self):
        #Setting up V^ab_cd
        Np = self.Np
        Nh = self.Nh
        Nm = self.Nm
        Nm2= self.Nm2
        bs = self.bs
        
        #Setting up p,q indices
        AB = arange(Np**2)
        P = AB%Np + Nh
        Q = AB//Np+ Nh
        
        #identify different outcomes of summation of two states K-quantum numbers
        idents_pp = self.kdplq(bs.states, P,Q,1)+self.kdplq(bs.states, P,Q,2)*Nm+self.kdplq(bs.states, P,Q,3)*Nm2
        self.idents_pp = idents_pp
        
        #map unique outcomes
        uniques = unique(idents_pp)
        self.lVpppp = []
        for e in uniques:
            T = where(idents_pp==e)[0]
            O = ones(len(T), dtype = int)
            
            t0 = kron(O, T)
            t1 = kron(T, O)

            a = t0%Np  + Nh
            b = t0//Np + Nh
            c = t1%Np  + Nh
            d = t1//Np + Nh
            a_ = bs.states[a,:]
            b_ = bs.states[b,:]
            c_ = bs.states[c,:]
            d_ = bs.states[d,:]
            data = bs.V(a_.T,b_.T,c_.T,d_.T)
            if len(T) != 0:
                data.shape = (len(T), len(T))
                self.lVpppp.append([data, T, T]) #store array and corresponding rows in dot product (dot(Vpppp, M))
        return self.lVpppp
    def sVpppp2(self):
        #Setting up V^ab_cd
        Np = self.Np
        Nh = self.Nh
        Nm = self.Nm
        Nm2= self.Nm2
        bs = self.bs
        
        #Setting up p,q indices
        AB = arange(Np**2)
        P = AB%Np + Nh
        Q = AB//Np+ Nh
        idents_pp = bs.states[P,1:]+bs.states[Q,1:]
        idents_pp = sum(array([1,Nm, Nm2, Nm*Nm2]) * idents_pp, 1)
                
        #identify different outcomes of summation of two states K-quantum numbers
        #idents_pp = self.kdplq(bs.states, P,Q,1)+self.kdplq(bs.states, P,Q,2)*Nm+self.kdplq(bs.states, P,Q,3)*Nm2
        #self.idents_pp = idents_pp
        
        #map unique outcomes
        uniques = unique(idents_pp)
        self.lVpppp = []
        for e in uniques:
            T = where(idents_pp==e)[0]
            O = ones(len(T), dtype = int)
            
            t0 = kron(O, T)
            t1 = kron(T, O)

            a = t0%Np  + Nh
            b = t0//Np + Nh
            c = t1%Np  + Nh
            d = t1//Np + Nh
            a_ = bs.states[a,:]
            b_ = bs.states[b,:]
            c_ = bs.states[c,:]
            d_ = bs.states[d,:]
            data = bs.V(a_.T,b_.T,c_.T,d_.T)
            if len(T) != 0:
                data.shape = (len(T), len(T))
                self.lVpppp.append([data, T, T]) #store array and corresponding rows in dot product (dot(Vpppp, M))
        return self.lVpppp
                
    def sVhhhh(self):
        #Setting up V^kl_ij
        #Important for the L2 calculation; this matrix is symmetric (V^T = V) so that only reorganization of the amplitudes is necessary
        Np = self.Np
        Nh = self.Nh
        Nm = self.Nm
        Nm2= self.Nm2
        bs = self.bs
        
        #Setting up p,q indices
        IJ = arange(Nh**2)
        P = IJ%Np 
        Q = IJ//Np
        
        #identify different outcomes of summation of two states K-quantum numbers
        idents_hh = self.kdplq(bs.states, P,Q,1)+self.kdplq(bs.states, P,Q,2)*Nm+self.kdplq(bs.states, P,Q,3)*Nm2
        self.idents_hh = idents_hh
        #map unique outcomes
        uniques = unique(idents_hh)
        T0 = array([], dtype = int)
        T1 = array([], dtype = int)
        #ret = zeros((Np**2, Np**2))
        self.lVhhhh = []
        for e in uniques:
            T = where(idents_hh==e)[0]
            O = ones(len(T), dtype = int)
            
            t0 = kron(O, T)
            t1 = kron(T, O)

            k = t0%Nh
            l = t0//Nh
            i = t1%Nh
            j = t1//Nh
            k_ = bs.states[k,:]
            l_ = bs.states[l,:]
            i_ = bs.states[i,:]
            j_ = bs.states[j,:]
            data = bs.V(k_.T,l_.T,i_.T,j_.T)
            if len(T) != 0:
                data.shape = (len(T), len(T))
                self.lVhhhh.append([data, T, T]) #store array and corresponding rows in dot product (dot(Vpppp, M))
        return self.lVhhhh

    def sVhhpp(self):        
        Np = self.Np
        Nh = self.Nh
        Nm = self.Nm
        Nm2= self.Nm2
        bs = self.bs

        T0 = array([], dtype = int)
        T1 = array([], dtype = int)

        P = arange(Nh**2)
        I = P%Nh 
        J = P//Nh
        
        Q = arange(Np**2)
        A = Q%Np + Nh
        B = Q//Np + Nh
        
        
        idents_hh = bs.states[I,1:]+bs.states[J,1:]
        idents_hh = sum(array([1,Nm, Nm2, Nm*Nm2]) * idents_hh, 1)
        
        idents_pp = bs.states[A,1:]+bs.states[B,1:]
        idents_pp = sum(array([1,Nm, Nm2, Nm*Nm2]) * idents_pp, 1)
        
        uniques = unique(append(unique(idents_pp),unique(idents_hh)))
        #print uniques

        self.lVhhpp = []
        for e in uniques:
            p = where(idents_hh == e)[0]
            q = where(idents_pp == e)[0]
            if len(p) != 0 and len(q) != 0:
                t0 = kron(ones(len(q), dtype = int), p)
                t1 = kron(q, ones(len(p), dtype = int))
                T0 = append(T0, t0)
                T1 = append(T1, t1) 
                           
                i = t0%Nh
                j = t0//Nh
                a = t1%Np +Nh
                b = t1//Np +Nh
                
                i_ = bs.states[i,:]
                j_ = bs.states[j,:]
                a_ = bs.states[a,:]
                b_ = bs.states[b,:]  
                data = bs.V(i_.T,j_.T,a_.T,b_.T)  
                data.shape = (len(p), len(q))
                self.lVhhpp.append([data, p, q])

        i = T0%Nh
        j = T0//Nh 
        a = T1%Np + Nh
        b = T1//Np + Nh
        
        a_ = bs.states[a,:]
        b_ = bs.states[b,:]
        i_ = bs.states[i,:]
        j_ = bs.states[j,:]
        data = bs.V(i_.T,j_.T,a_.T,b_.T)

        return coo_matrix((data, (T0, T1)), shape=(Nh**2, Np**2)).toarray()
        #return self.lVhhpp


    def sVpphh(self):        
        Np = self.Np
        Nh = self.Nh
        Nm = self.Nm
        Nm2= self.Nm2
        bs = self.bs

        T0 = array([], dtype = int)
        T1 = array([], dtype = int)

        P = arange(Np**2)
        A = P%Np + Nh
        B = P//Np + Nh
        
        Q = arange(Nh**2)
        I = Q%Nh
        J = Q//Nh
        
        
        idents_hh = bs.states[I,1:]+bs.states[J,1:]
        idents_hh = sum(array([1,Nm, Nm2, Nm*Nm2]) * idents_hh, 1)
        idents_pp = bs.states[A,1:]+bs.states[B,1:]
        idents_pp = sum(array([1,Nm, Nm2, Nm*Nm2]) * idents_pp, 1)
        
        uniques = unique(append(unique(idents_pp),unique(idents_hh)))
        #print uniques

        self.lVpphh = []
        for e in uniques:
            p = where(idents_pp == e)[0]
            q = where(idents_hh == e)[0]
            if len(p) != 0 and len(q) != 0:
                t0 = kron(ones(len(q), dtype = int), p)
                t1 = kron(q, ones(len(p), dtype = int))
                T0 = append(T0, t0)
                T1 = append(T1, t1) 
                           
                a = t0%Np + Nh
                b = t0//Np + Nh
                i = t1%Nh 
                j = t1//Nh
                
                a_ = bs.states[a,:]
                b_ = bs.states[b,:]
                i_ = bs.states[i,:]
                j_ = bs.states[j,:]  
                data = bs.V(a_.T,b_.T,i_.T,j_.T)  
                data.shape = (len(p), len(q))
                self.lVpphh.append([data, p, q])

        a = T0%Np + Nh
        b = T0//Np + Nh
        i = T1%Nh 
        j = T1//Nh
        
        a_ = bs.states[a,:]
        b_ = bs.states[b,:]
        i_ = bs.states[i,:]
        j_ = bs.states[j,:]
        data = bs.V(a_.T,b_.T,i_.T,j_.T)

        return coo_matrix((data, (T0, T1)), shape=(Np**2, Nh**2)).toarray()
        #return self.lVpphh
        
    def sVhpph(self):
        Np = self.Np
        Nh = self.Nh
        Nm = self.Nm
        Nm2= self.Nm2
        bs = self.bs

        T0 = array([], dtype = int)
        T1 = array([], dtype = int)

        P = arange(Np*Nh)
        I = P%Nh 
        A = P//Nh + Nh
        
        Q = arange(Np*Nh)
        B = Q%Np + Nh 
        J = Q//Np
        
        
        idents_hp = bs.states[I,1:]+bs.states[A,1:]
        idents_hp = sum(array([1,Nm, Nm2, Nm*Nm2]) * idents_hp, 1)
        idents_ph = bs.states[B,1:]+bs.states[J,1:]
        idents_ph = sum(array([1,Nm, Nm2, Nm*Nm2]) * idents_ph, 1)
        
        uniques = unique(append(unique(idents_hp),unique(idents_ph)))
        #print uniques

        self.lVhpph = []
        for e in uniques:
            p = where(idents_hp == e)[0]
            q = where(idents_ph == e)[0]
            if len(p) != 0 and len(q) != 0:
                t0 = kron(ones(len(q), dtype = int), p)
                t1 = kron(q, ones(len(p), dtype = int))
                #T0 = append(T0, t0)
                #T1 = append(T1, t1)  
                
                i = t0%Nh
                a = t0//Nh + Nh
                b = t1%Np + Nh
                j = t1//Np 
                
                i_ = bs.states[i,:]
                a_ = bs.states[a,:]
                b_ = bs.states[b,:]
                j_ = bs.states[j,:]
                data = bs.V(i_.T,a_.T,b_.T,j_.T)
                data.shape = (len(p), len(q))
                self.lVhpph.append([data, p,q])
                         
        """
        i = T0%Nh
        a = T0//Nh + Nh
        b = T1%Np + Nh
        j = T1//Np 
        
        i_ = bs.states[i,:]
        a_ = bs.states[a,:]
        b_ = bs.states[b,:]
        j_ = bs.states[j,:]
        data = bs.V(i_.T,a_.T,b_.T,j_.T)
        return coo_matrix((data, (T0, T1)), shape=(Np*Nh, Nh*Np)).toarray()
        """
        return self.lVhpph

    def ident(self,v):
        #A unique identifying integer for the momentum combinations
        return v[0] + v[1]*self.bs.Nm + v[2]*self.bs.Nm**2 + v[3]*self.bs.Nm**3    
    


def setup_Vpppp(bs):
    #Old linear approach
    Np = bs.nstates-bs.nparticles
    Nh = bs.nparticles
    Z = zeros((Np**2, Np**2))
    for a in range(Np):
        for b in range(Np):
            for c in range(Np):
                for d in range(Np):
                    Z[a + b*Np, c+d*Np] = bs.v(a + Nh, b + Nh,c + Nh, d + Nh)
    return Z

def setup_Vpphh(bs):
    #Old linear approach
    Np = bs.nstates-bs.nparticles
    Nh = bs.nparticles
    Z = zeros((Np**2, Nh**2))
    for a in range(Np):
        for b in range(Np):
            for i in range(Nh):
                for j in range(Nh):
                    Z[a + b*Np, i+j*Nh] = bs.v(a + Nh, b + Nh,i, j)
    return Z

    
def composite_dot(M1,M2):
    Z = zeros((len(M2), len(M2[0])))
    for i in range(len(M1)):
        Z[M1[i][1]] = dot(M1[i][0], M2[M1[i][1]])
    return Z
            

bs = electronbasis(3,1.0,14)
op = optiV(bs)
C = consistency_test()

Np = bs.nstates-bs.nparticles
Nh = bs.nparticles
print "************************************************"
print "**"
print "**    Number of states:", Np+Nh
print "** Number of particles:", Np
print "**"
print "************************************************"
T2 = flexitude(Np, Nh)
T2.M = random.randint(0,10,(Np**2, Nh**2))




vhhpp = op.sVhhpp() #works
#C.test_Vhhpp(vhhpp, bs)
imshow(vhhpp)
show()
"""
#These are the only ones needed
t0 = clock()
vpppp = op.sVpppp() #works
vhpph = op.sVhpph() #works
vhhpp = op.sVhhpp() #works
vhhhh = op.sVhhhh() #untested, but assumed to work
t1 = clock()
r = composite_dot(vpppp, T2.ab_ij())
r = composite_dot(vhhhh, T2.ij_ab())
r = composite_dot(vhhpp, T2.ab_ij())
r = composite_dot(vhpph, T2.ai_bj())
t2 = clock()
print "Time spent initializing V            :", t1-t0
print "Time spent doing some amplitude calcs:", t2-t1
"""