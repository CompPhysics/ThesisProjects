from numpy import *
from matplotlib.pyplot import *
from time import *


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
        
        self.setup_matrices_optimized()
        self.Tpphh = random.uniform(0,1,(self.Np**2, self.Nh**2))
    def ident(self,v):
        #A unique identifying integer for the momentum combinations
        return v[0] + v[1]*bs.Nm + v[2]*bs.Nm**2 + v[3]*bs.Nm**3
        
    def setup_pppp(self):
        Np = self.Np
        ppa = kron(arange(Np**2)//Np, ones((Np**2,1), dtype = int))
        ppb = kron(arange(Np**2)%Np, ones((Np**2,1), dtype = int))
        
        combs_pp = 20000*ones((Np**2,Np**2), dtype = int) #arbitrary large number since identifier will include zeros
        idents = zeros((Np**2))
        for p in range(Np):
            for q in range(Np):
                v = bs.states[p][1:5]+bs.states[q][1:5]
                iv =  self.ident(v)
                combs_pp[p + q*Np, :] = iv  #this one should not be zero, as most elements in array is already zero, or?
                idents[p+q*Np] = iv
        spectrum = unique(idents)
        print spectrum
        combs_pp[combs_pp!=combs_pp.T]=20000 #identify each pair of quantum numbers sharing the same added momentum
        
        #create all submatrices 
        submat = []
        submat_dim = zeros((len(spectrum),2))
        submat_castback = []
        
        for m in range(len(spectrum)):
            t = where(combs_pp==spectrum[m])
            t = array(t)
            
            #calculate all contribs from t
            
            
            Ntx = len(unique(t[0]))
            Nty = len(unique(t[1]))
            submat_dim[m,:] = array([Ntx,Nty])
            t0 = t[0].reshape(Ntx,Nty)
            t1 = t[1].reshape(Ntx,Nty)
            
            V = self.Vpppp[t0,t1]
            
            #  Extract corresponding subselection in T
            #T = self.Tpphh[t1[0],:] #select only rows from T corresponding to elements in V
            #  Broadcast values back to 
            #self.Tpphh[t1[0],:] += dot(V,T) 
            
            #set up interactions here? (optimization possibility)
            
            submat.append(V)  #replace combs_pp -> V_pp
            submat_castback.append([t0,t1]) #in case of back-transform requirements, broadcasts values back to their original indices
        #print submat
        #print submat_castback
        self.sub_pppp = submat
        self.sub_pppp_castback = submat_castback
        self.sub_ppp_dim = submat_dim

            
        
        
    def maprange(self):
        pq_occ = zeros((self.Ns, self.Ns)) #combinations already accounted for
        Np = self.Np
        Ns = self.Ns
        Nh = self.Nh
        Nm = self.bs.Nm #Max possible momentum
        #combs_pq = unique(combs_pq) #retain only unique resuts
        
        hhi = kron(arange(Nh**2)//Nh, ones((Nh**2, 1), dtype = int))
        hhj = kron(arange(Nh**2)%Nh, ones((Nh**2, 1), dtype = int))
        ppa = kron(arange(Np**2)//Np, ones((Np**2,1), dtype = int))
        ppb = kron(arange(Np**2)%Np, ones((Np**2,1), dtype = int))
        #print i
        #print j

        
        combs_pp = zeros((Np**2,Np**2), dtype = int)      
        for p in range(Np):
            for q in range(Np):
                v = bs.states[p][1:5]+bs.states[q][1:5]
                combs_pp[p + q*Np, :] = self.ident(v)
        
        combs_pp[combs_pp!=combs_pp.T]=0 #identify each pair of quantum numbers sharing the same added momentum
        
        #decompose into smaller matrices
        
        
        figure(1)
        title("V_pppp")
        imshow(combs_pp)
        
        
        combs_hh = zeros((Nh**2,Nh**2), dtype = int)      
        for p in range(Nh):
            for q in range(Nh):
                v = bs.states[p][1:5]+bs.states[q][1:5]
                combs_hh[p + q*Nh, :] = self.ident(v)
        #combs_pq = unique(combs_pq) #retain only unique resuts
        combs_hh[combs_hh!=combs_hh.T]=0 #identify each pair of quantum numbers sharing the same added momentum
        figure(2)
        title("V_hhhh")
        t = where(combs_hh==1)
        t= array(t)
        Ntx = len(unique(t[0]))
        Nty = len(unique(t[1]))
        print Ntx, Nty
        print t[0]
        t0 = t[0].reshape(Ntx,Nty)
        t1 = t[1].reshape(Ntx,Nty)
        print t0, t1
        print combs_hh[t0,t1]
        #try to use "choose" (numpy)
        
        imshow(combs_hh)  
        show()     

        """
        combs_rs = zeros((Nrs), dtype = int)
        for r in rr:
            for s in rs:
                v = bs.states[r][1:5]+bs.states[s][1:5]
                combs_rs[:,s + r*Nrs] = self.ident(v)
        #combs_rs = unique(combs_rs) #retain only unique resuts
        
        imshow(combs_rs==combs_pq)
        
        
        self.blocks = []
        n = 0

        for e in range(len(combs)):
            indices = []
            for q in rq:
                for p in rp:
                    v = bs.states[p][1:5]+bs.states[q][1:5]
                    if self.ident(v)==combs[e]:
                        indices.append([p,q])
            self.blocks.append([zeros((len(indices), len(indices)), dtype = int),zeros((len(indices), len(indices)), dtype = int)])
            for pq in range(len(indices)):
                for rs in range(len(indices)):
                    self.blocks[e][0][pq,rs] = indices[pq][0] + Ns*indices[pq][1]
                    self.blocks[e][1][pq,rs] = indices[rs][0] + Ns*indices[rs][1]
            #Identify indices below and above FL

        self.blockdiag = zeros((self.Ns**2, self.Ns**2))
        self.dims = zeros(len(self.blocks), dtype = int)
        n = 0
        for e in range(len(self.blocks)):
            N = len(self.blocks[e][0][0])
            self.dims[e] = N
            self.blockdiag[n:N,n:N] += self.V[self.blocks[e]]
            #n += N
        print self.dims, sum(self.dims)
        
        imshow(self.blockdiag)
        show()
        """
    def map2(self):
        pq_occ = zeros((self.Ns, self.Ns)) #combinations already accounted for
        Np = self.Np
        Ns = self.Ns
        Nh = self.Nh
        Nm = self.bs.Nm #Max possible momentum
        
        combs = zeros((Ns**2), dtype = int)
        for p in range(Ns):
            for q in range(Ns):
                v = bs.states[p][1:5]+bs.states[q][1:5]
                combs[p + q*Ns] = self.ident(v)
        combs = unique(combs) #retain only unique resuts
        self.blocks = []
        n = 0

        for e in range(len(combs)):
            indices = []
            for q in range(Ns):
                for p in range(Ns):
                    v = bs.states[p][1:5]+bs.states[q][1:5]
                    if self.ident(v)==combs[e]:
                        indices.append([p,q])
            self.blocks.append([zeros((len(indices), len(indices)), dtype = int),zeros((len(indices), len(indices)), dtype = int)])
            for pq in range(len(indices)):
                for rs in range(len(indices)):
                    self.blocks[e][0][pq,rs] = indices[pq][0] + Ns*indices[pq][1]
                    self.blocks[e][1][pq,rs] = indices[rs][0] + Ns*indices[rs][1]
            #Identify indices below and above FL

        self.blockdiag = zeros((self.Ns**2, self.Ns**2))
        self.dims = zeros(len(self.blocks), dtype = int)
        n = 0
        for e in range(len(self.blocks)):
            N = len(self.blocks[e][0][0])
            self.dims[e] = N
            self.blockdiag[n:N,n:N] += self.V[self.blocks[e]]
            #n += N
        print self.dims, sum(self.dims)
        
        imshow(self.blockdiag)
        show()
                
        

        
    def map(self):
        Np = self.Np
        Ns = self.Ns
        Nh = self.Nh
        Nm = self.bs.Nm #Max possible momentum
        
        alpha = zeros((Nm*2 + 1, Nm*2 + 1, Nm*2 + 1))
        beta  = zeros((Nm*2 + 1, Nm*2 + 1, Nm*2 + 1))
        aNb = zeros(((Nm*2 + 1)**3 ))
        n = 0
        combs = []
        for p in range(Ns):
            for q in range(Ns):
                v = bs.states[p][1:5]+bs.states[q][1:5]
                i = ident(v)
                combs.append([ident,p,q])
        combs.sort()
        print combs
        print len(combs)
    def generate_all(self):
        Ns = self.Ns
        self.V = zeros((Ns**2, Ns**2))
        for p in range(Ns):
            for q in range(Ns):
                for r in range(Ns):
                    for s in range(Ns):
                        self.V[p + q*Ns, r + s*Ns] = q+r+s #self.bs.v(p,q,r,s)
    def setup_matrices_optimized(self):
        t0 = clock()

        Nh = self.Nh
        Np = self.Np
        for i in range(Nh):
            for j in range(i,Nh):
                for a in range(Np):
                    for b in range(a,Np):
                        val = self.bs.v(i,j,a+Nh,b+Nh)
                        self.Vhhpp[i + j*Nh, a + b*Np] = val
                        self.Vhhpp[j + i*Nh, b + a*Np] = val
                        self.Vhhpp[j + i*Nh, a + b*Np] = -val
                        self.Vhhpp[i + j*Nh, b + a*Np] = -val
                        #self.Vhhpp[i + j*Nh, a + b*Np] = self.bs.v(i,j,a+Nh,b+Nh)
                        
                        #self.Vpphh[a + b*Np, i + j*Nh] = self.bs.v(a+Nh,b+Nh,i,j)
                        self.Vpphh[a + b*Np, i + j*Nh] = val
                        self.Vpphh[a + b*Np, j + i*Nh] = -val
                        self.Vpphh[b + a*Np, i + j*Nh] = -val
                        self.Vpphh[b + a*Np, j + i*Nh] = val
                        
                        #print self.bs.v(a+Nh,i,j,b+Nh) == self.bs.v(j,a+Nh,b+Nh,i)
                        
                        val = self.bs.v(a+Nh,i,j,b+Nh)
                        self.Vphhp[a + i*Np, j + b*Nh] = val
                        self.Vhpph[i + a*Nh, b + j*Np] = val
                        self.Vphhp[b + j*Np, i + a*Nh] = val             
                        self.Vhpph[j + b*Nh, a + i*Np] = val
                        
                        val = self.bs.v(j,a+Nh,b+Nh,i)
                        self.Vphhp[a + j*Np, i + b*Nh] = val
                        self.Vhpph[j + a*Nh, b + i*Np] = val
                        self.Vphhp[b + i*Np, j + a*Nh] = val
                        self.Vhpph[i + b*Nh, a + j*Np] = val
                        
                        
                        
                        
                        #val = self.bs.v(i,a+Nh,b+Nh,j)
                        #self.Vhpph[i + a*Nh, b + j*Np] = self.bs.v(i,a+Nh,b+Nh,j)
                        

                        
                        val = self.bs.v(a+Nh,b+Nh,i,j)
                        eps = self.bs.h(i,i) + self.bs.h(j,j) -self.bs.h(a+Nh,a+Nh) - self.bs.h(b+Nh,b+Nh)
                        
                        self.Tpphh[a + b*Np, i + j*Nh] = val/eps
                        self.Epphh[a + b*Np, i + j*Nh] = eps
                        
                        self.Tpphh[a + b*Np, j + i*Nh] = -val/eps
                        self.Epphh[a + b*Np, j + i*Nh] = eps
                        
                        self.Tpphh[b + a*Np, i + j*Nh] = -val/eps
                        self.Epphh[b + a*Np, i + j*Nh] = eps
                        
                        self.Tpphh[b + a*Np, j + i*Nh] = val/eps
                        self.Epphh[b + a*Np, j + i*Nh] = eps
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
                        
                        self.Vhhhh[k + l*Nh, i+ j*Nh] = val
                        self.Vhhhh[j + i*Nh, l+ k*Nh] = val
                        self.Vhhhh[j + i*Nh, k+ l*Nh] = -val
                        self.Vhhhh[i + j*Nh, l+ k*Nh] = -val
                        
                        
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

class transformer():
    def __init__(self, Np, Nh):
        self.Nh = Nh
        self.Np = Np
        self.init_index_reorg()
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
            
bs = electronbasis(2,1.0,14)
b = blocks(bs)
#b.generate_all()
b.setup_pppp()