from numpy import *
from time import *
from matplotlib.pyplot import *

class electronbasis():
    def __init__(self, N, rs, Nparticles):
        self.rs = rs
        self.states = []
        self.nstates = 0
        self.nparticles = Nparticles
        
        Nm = int(sqrt(N) + 1)
        self.Nm = N
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

class blockaccount():
    #Indexing and keeping tabs for block matrix
    #Class to set up V and T matrices
    def __init__(self, N, bs):
        self.N = N #size
        self.Ns = bs.nstates
        self.Np = bs.nparticles
        self.bs =bs
        self.N0 = 2*N
        self.N2 = N**2
        self.N3 = N**3
        self.Z = zeros((4*N)**3, dtype = int)
        self.Zmap = []
        for i in range((4*N)**3):
            self.Zmap.append([])
        self.colindex = zeros((self.Ns**2), dtype = int)
        self.V = zeros((self.Ns**2, self.Ns**2), dtype = float)
        self.T = zeros((self.Ns**2, self.Ns**2), dtype = float)
        #self.T = random.uniform(0,1,(self.Ns**2, self.Ns**2))
        
    def k_at(self, kx,ky,kz):
        return self.Z[(self.N0+kx) + (self.N0+ky)*self.N2 + (self.N0+kz)*self.N3]
    def k_set(self,kx,ky,kz, val):
        self.Zmap[(self.N0+kx) + (self.N0+ky)*self.N2 + (self.N0+kz)*self.N3].append(val)
    def V_at(self, p,q):
        return self.V[p + q*self.Ns, p + q*self.Ns]
    def V_set(self, p,q,val):
        self.V[p,q] = val  
    def T_at(self,a,b,i,j):
        #return the amplitude t^ab_ij
        tx = self.Zmap[a + self.Ns*b]
        ty = self.Zmap[i + self.Ns*j]
        return self.T[tx,ty]
    
    def T_set(self,a,b,i,j,val):
        tx = self.Zmap[a + self.Ns*b]
        ty = self.Zmap[i + self.Ns*j]
        self.T[tx,ty] = val

        
    def inverse_lookup(self, n):
        p = n%self.Ns
        q = n//self.Ns
        return p,q
    def sort_ks(self):
        for p in range(self.Ns):
            for q in range(self.Ns):
                val = self.bs.states[p][1:5] + self.bs.states[q][1:5]
                self.k_set(int(val[0]),int(val[1]),int(val[2]), p+q*self.Ns) #Add p,q to the corresponding K,M
                
        Ncount = 0
        self.Vblocks = []
        self.Tblocks = []
        self.blockmap = zeros((self.Ns,self.Ns), dtype = int) #a mapping from index to block
        self.indexmap = zeros((self.Ns,self.Ns), dtype = int) #a mapping form index to element in block
        for i in range(len(self.Zmap)):
            N = len(self.Zmap[i]) #number of elements in block
            self.Vblocks.append(zeros((N,N)))
            self.Tblocks.append(zeros((N,N)))
            for e in range(N):
                for u in range(e,N):
                    #initialize interactions
                    
                    p,q= self.inverse_lookup(self.Zmap[i][e])
                    r,s= self.inverse_lookup(self.Zmap[i][u])
                    val = self.bs.v(p,q,r,s)
                    self.Vblocks[i][e,u] = val
                    self.Vblocks[i][u,e] = val
                    self.V_set(Ncount + e, Ncount +u, val)
                    self.V_set(Ncount + u, Ncount +e, val)
                    
                    #save mapping to separate blocks and element
                    self.blockmap[p,q]= i #element p,q is in block i
                    self.indexmap[p,q]= e #element p,q is at column e in i
                    
                    
                    #initialize amplitudes
                    if r >= self.Np and s >= self.Np and p < self.Np and q < self.Np:
                        ampl = bs.v(r,s,p,q)/float(-bs.h(p,p) - bs.h(q,q) + bs.h(r,r) + bs.h(s,s))
                        #self.T_set(r,s,p,q,ampl)
                        self.Tblocks[i][e,u] = ampl
                        self.Tblocks[i][u,e] = ampl
                        self.T[Ncount + e, Ncount + u] = ampl
                        self.T[Ncount + u, Ncount + e] = ampl
                        
                        
                        #self.T_set(p,q,s,r,-ampl)
                        #self.T_set(q,p,r,s,-ampl)
                        #self.T_set(q,p,s,r,ampl)
                        
                        #self.T_set(r,s,p,q,ampl)
                        #self.T_set(s,r,p,q,-ampl)
                        #self.T_set(r,s,q,p,-ampl)
                        #self.T_set(s,r,q,p,ampl)
                        
            Ncount += N
        print self.colindex, Ncount == self.Ns**2
        print self.T.max()
        #print self.inverse_lookup(self.colindex[0])
        #print self.V
    def multiply(self):
        N = len(self.Zmap)
        self.TV = []
        for i in range(self.N):
            self.TV.append(self.Tblocks[i]*self.Vblocks[i])
    
    #Access block elements
    def blockT_at(self,p,q,r,s):
        #P,Q = p + q*self.Ns, r+s*self.Ns
        val = 0
        Ip = self.blockmap[p,q]
        Iq = self.blockmap[r,s]
        if Ip==Iq:
            P = self.indexmap[p,q]
            Q = self.indexmap[r,s]
            val = self.Tblocks[Ip][P,Q]
        #Eq = self.indexmap[P,Q]
        return val
    
    def blockV_at(self,p,q,r,s):
        #P,Q = p + q*self.Ns, r+s*self.Ns
        val = 0
        Ip = self.blockmap[p,q]
        Iq = self.blockmap[r,s]
        if Ip==Iq:
            P = self.indexmap[p,q]
            Q = self.indexmap[r,s]
            val = self.Vblocks[Ip][P,Q]
        #Eq = self.indexmap[P,Q]
        return val
    
    
                
        

        
                

t0 = clock()
bs = electronbasis(3,1,14)
print "Time spent setting up basis:", clock()-t0

#session = CCDsolver(bs)
#session.sort_states()

H = 2
B = blockaccount(H, bs)
t0 = clock()
B.sort_ks()
print "Time spent setting up matrix:", clock()-t0

t0 = clock()
B.multiply()
print "Time spent on listwise multiplication:", clock()-t0

t0 = clock()
VT = B.V*B.T
print "Time spent on full multiplication:", clock()-t0

t0 = clock()
for i in range(14):
    for j in range(14):
        for a in range(14,bs.nstates):
            for b in range(14,bs.nstates):
                if B.blockV_at(i,j,a,b) != bs.v(i,j,a,b):
                    print i,j,a,b, "Not true"
print "Scanned all interactions in:", clock()-t0, "seconds."



#print 
#print B.Zmap
#imshow(B.T, cmap = "RdGy")
#show()
