from numpy import *
from time import *
from matplotlib.pyplot import *
from scipy.sparse import csr_matrix, coo_matrix

#Main goal for this implementation: avoid poor design choices

class electronbasis():
    def __init__(self, N, rs, Nparticles):
        self.rs = rs
        self.states = []
        self.nstates = 0
        self.nparticles = Nparticles
        self.nshells = N - 1
        self.Nm = N + 1
        
        self.k_step = 2*(self.Nm + 1)
        Nm = N
        n = 0 #current shell
        ene_integer = 0
        while n <= self.nshells:
            is_shell = False
            for x in range(-Nm, Nm+1):
                for y in range(-Nm, Nm+1):
                    for z in range(-Nm,Nm+1):
                        e = x*x + y*y + z*z
                        if e   == ene_integer:
                            is_shell = True
                            self.nstates += 2
                            self.states.append([e, x,y,z,1])
                            self.states.append([e, x,y,z, -1])
                            
            if is_shell:
                n += 1
            ene_integer += 1
        self.L3 = (4*pi*self.nparticles*self.rs**3)/3.0
        self.L2 = self.L3**(2/3.0)
        self.L = pow(self.L3, 1/3.0)
        
        #for i in range(self.nstates):
        #    self.states[i][0] *= 2*(pi**2)/self.L**2 #Multiplying in the missing factors in the single particle energy
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
    def V2(self, kp,kq,kr,ks):
        #k = (energy, kx, ky, kz, ms)
        # Vectorized interaction
        # This function assumes the the first criterion (comment line below) has been asserted to be true

        kdplus = 4*pi/self.L3 #(kp[:,1]+kq[:,1]==kr[:,1]+ks[:,1])*(kp[:,2]+kq[:,2]==kr[:,2]+ks[:,2])*(kp[:,3]+kq[:,3]==kr[:,3]+ks[:,3])*4*pi/self.L3#d_k+k k+k

        kdspin1 = (kp[:,4]==kr[:,4])*(kq[:,4]==ks[:,4])*1
        kdwave1 = abs((kp[:,1]==kr[:,1])*(kp[:,2]==kr[:,2])*(kp[:,3]==kr[:,3])-1)
        absdiff2_1 = ((kr[:,1]-kp[:,1])**2+(kr[:,2]-kp[:,2])**2+(kr[:,3]-kp[:,3])**2) #absdiff2
        term1=(4.0*absdiff2_1*pi**2)/self.L2
        term1[term1==0] = 1
        
        kdspin2 = (kp[:,4]==ks[:,4])*(kq[:,4]==kr[:,4])*1
        kdwave2 = abs((kp[:,1]==ks[:,1])*(kp[:,2]==ks[:,2])*(kp[:,3]==ks[:,3])-1)
        absdiff2_2 = ((ks[:,1]-kp[:,1])**2+(ks[:,2]-kp[:,2])**2+(ks[:,3]-kp[:,3])**2) #absdiff2
        term2=(4.0*absdiff2_2*pi**2)/self.L2
        term2[term2==0] = 1
        
        return kdplus*(kdspin1*kdwave1/term1 - kdspin2*kdwave2/term2)
        
    def V(self, kp,kq,kr,ks):
        #k = (energy, kx, ky, kz, ms)
        # Vectorized interaction
        # This function assumes the the first criterion (comment line below) has been asserted to be true

        kdplus = (kp[1,:]+kq[1,:]==kr[1,:]+ks[1,:])*(kp[2,:]+kq[2,:]==kr[2,:]+ks[2,:])*(kp[3,:]+kq[3,:]==kr[3,:]+ks[3,:])*4*pi/self.L3#d_k+k k+k //FACTOR 2 originally 4

        kdspin1 = (kp[4,:]==kr[4,:])*(kq[4,:]==ks[4,:])*1
        kdwave1 = abs((kp[1,:]==kr[1,:])*(kp[2,:]==kr[2,:])*(kp[3,:]==kr[3,:])-1)
        absdiff2_1 = ((kr[1,:]-kp[1,:])**2+(kr[2,:]-kp[2,:])**2+(kr[3,:]-kp[3,:])**2) #absdiff2
        term1=(4.0*absdiff2_1*pi**2)/self.L2
        term1[term1==0] = 1
        
        kdspin2 = (kp[4,:]==ks[4,:])*(kq[4,:]==kr[4,:])*1
        kdwave2 = abs((kp[1,:]==ks[1,:])*(kp[2,:]==ks[2,:])*(kp[3,:]==ks[3,:])-1)
        absdiff2_2 = ((ks[1,:]-kp[1,:])**2+(ks[2,:]-kp[2,:])**2+(ks[3,:]-kp[3,:])**2) #absdiff2
        term2=(4.0*absdiff2_2*pi**2)/self.L2
        term2[term2==0] = 1
        
        return kdplus*(kdspin1*kdwave1/term1 - kdspin2*kdwave2/term2)
    def v_preapproved(self, p,q,r,s):
        terms = 0.0
        term1 = 0.0
        term2 = 0.0
        val = 1.0/self.L3

        if self.kdspin(p,r)*self.kdspin(q,s)==1:
            if self.kdwave(p,r) != 1.0:
                #term1=(4*self.absdiff2(r,p)*pi**2)/self.L2
                #terms += 1.0/term1
                term1 = self.L2/(pi*self.absdiff2(r,p))
        if self.kdspin(p,s)*self.kdspin(q,r)==1:
            if self.kdwave(p,s) != 1.0:
                #term2=(4*self.absdiff2(s,p)*pi**2)/self.L2
                #terms -= 1.0/term2
                term2 = self.L2/(pi*self.absdiff2(s,p))
        return val*(term1-term2)

    def v(self,p,q,r,s):
        #Two body interaction
        #To optimize bottleneck: vectorize this function ! (remove if-tests)
        val = 0
        terms = 0.0
        term1 = 0.0
        term2 = 0.0
        kdpl = self.kdplus(p,q,r,s)
        if kdpl != 0:
            val = 1.0/self.L3
            
            if self.kdspin(p,r)*self.kdspin(q,s)==1:
                if self.kdwave(p,r) != 1.0:
                    #term1=(4*self.absdiff2(r,p)*pi**2)/self.L2
                    #terms += 1.0/term1
                    term1 = self.L2/(pi*self.absdiff2(r,p))
            if self.kdspin(p,s)*self.kdspin(q,r)==1:
                if self.kdwave(p,s) != 1.0:
                    #term2=(4*self.absdiff2(s,p)*pi**2)/self.L2
                    #terms -= 1.0/term2
                    term2 = self.L2/(pi*self.absdiff2(s,p))
        return val*(term1-term2)

    
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
    def unique(self, p):
        return self.states[p,1] + self.states[p,2]*self.k_step + self.states[p,3]*self.k_step**2 + self.states[p,4]*self.k_step**3
        #return self.states[p,4] + self.states[p,3]*self.k_step + self.states[p,1]*self.k_step**2 + self.states[p,2]*self.k_step**3 + self.k_step**4




class blockmap():
    """
    This tensorclass keeps track of and stores
    (1) The dense blocks of amplitudes or interactions
    (2) The assosciated unique identifiers resulting from conservation of quantum numbers
    (3) The reorganization patterns needed to align matrices when performing contractions
    """
    def __init__(self, basis, Np, Nq, Nr, Ns):
        self.bs = basis
        #print "Number of states:", self.bs.nstates
        self.k_step = 2*(self.bs.Nm + 1) #steplength of momentum vector, used to uniquely identify combinations of quantum numbers
        self.elements = [] #1D array where each element is stored, common to all blocks
        self.elements = array([], dtype = float)
        self.blockmap = [] #1D array where for each element, the mapped blocks and indices are stored, so that 
        self.blocks = [] #nested list containing blocks of pointers to elements
        self.configs = [] #nested list containing the incoming and outgoing quantum numbers of blocks
        self.blocklengths = []
        self.ordering = []
        
        self.Np = Np
        self.Nq = Nq
        self.Nr = Nr
        self.Ns = Ns
        
        self.a = 0 #timers for profiling
        self.b = 0
        self.c = 0
        self.d = 0
        
        #self.Ns = [self.Np, self.Nq, self.Nr, self.Ns]
        
        self.configurations = []
        self.amplitude = False
    
    ###########################
    # Transform between 1D to /from 4D
    ############################
    def to(self, p,q,r,s):
        #translate to transformed index 
        return p + q*self.Np + r*self.Np*self.Nq + s*self.Np*self.Nq*self.Nr
    def of(self, i):
        #translate from transformed index
        s = i//(self.Np*self.Nq*self.Nr)
        r = (i-s*self.Np*self.Nq*self.Nr)//(self.Np*self.Nq)
        q = (i-s*self.Np*self.Nq*self.Nr-r*self.Np*self.Nq)//self.Np
        p = (i-s*self.Np*self.Nq*self.Nr-r*self.Np*self.Nq-q*self.Np)       
        return [p,q,r,s]
        
    ###########################
    # Control elements
    ############################

    def zeros(self):
        #set all elements to zero
        self.elements *= 0
    def init(self):
        for i in range(len(self.elements)):
            p,q,r,s = self.of(self.elements[i]) #a,b,i,j
            self.elements[i] = self.bs.v(p,q,r+self.Np,s+self.Nq)
    
    def init_as_t2amp(self):
        #initialize elements as t2 amplitude
        self.energy_denom = zeros(len(self.elements), dtype = float)
        self.reset_elements = zeros(len(self.elements), dtype = float)
        
        for i in range(len(self.elements)):
            p,q,r,s = self.of(self.elements[i]) #a,b,i,j
            v = self.bs.v(p+self.Nr,q+self.Nr,r,s)
            self.elements[i] = v #/(bs.states[r,0] + bs.states[s,0] - bs.states[self.Nr + p,0] - bs.states[self.Nr + q,0])
            self.reset_elements[i] = v
            self.energy_denom[i] = bs.states[r,0] + bs.states[s,0] - bs.states[self.Nr + p,0] - bs.states[self.Nr + q,0]
        self.elements /= self.energy_denom
    def reset_amplitude(self):
        self.elements = self.reset_elements
    def init_as_pppp(self):
        #initialize elements as t2 amplitude
        self.energy_denom = zeros(len(self.elements), dtype = float)
        for i in range(len(self.elements)):
            p,q,r,s = self.of(self.elements[i]) #a,b,i,j
            self.elements[i] = self.bs.v(p+14,q+14,r+14,s+14)
    def energy_division(self):
        for i in range(len(self.elements)):
            self.elements[i] /= self.energy_denom[i]

    
    #######################################
    ##
    ##   The different subindex initialization routines
    ##
    #######################################    
    
    def unpack(self, row = arange(10*12), L = [[0,1], [1,1]], MN=[10, 12,14,5]):  
        #Unpack a compressed index of any size
        M = [1]
        mn = 1
        for i in range(len(L)):
            mn*=MN[L[i][0]]
            M.append(mn)
        M.reverse()
        indices = []
        for i in range(len(L)):
            P = row.copy()
            for e in range(i):
                P -= indices[e]*M[e+1]
            indices.append(P//M[i+1])
        indices.reverse()
        return indices
    
    def map_regions(self, L = [[0,1, 14],[1,1, 14]], R=[[2,1,0],[3,1,0]]):
        
        ###########################################################
        ##
        ##  Set up configuration specified by L, R, without any storage (to later generate elements on the fly)
        ##
        ###########################################################
        
        #(1) label rows and columns
        Ns = [self.Np, self.Nq, self.Nr, self.Ns]
        #(2) determine dimensions and indices of transformed matrix
        Nrows = 1
        for i in range(len(L)):
            Nrows*=Ns[L[i][0]]
        rows = arange(Nrows)
        Ncols = 1
        for i in range(len(R)):
            Ncols*=Ns[R[i][0]]
        cols = arange(Ncols)
        
        #(3) Transform matrix indices to element indices in tensor
        left = self.unpack(rows, L, Ns)  ##indices in the transformed matrix, do I know in which order these occurs? (p-q-r-s) yes, from L
        right = self.unpack(cols, R, Ns)
        pqrs = []
        PQRS = [0,0,0,0]
        for i in range(len(left)):
            pqrs.append(left[i])
            PQRS[L[i][0]] = left[i]
        for i in range(len(right)):
            pqrs.append(right[i])
            PQRS[R[i][0]] = right[i]
        
        #(3) Identify blocks
        #basically ident(row) == ident(col) produces the blocks
        """
        k = self.bs.k_step 
        LHS = zeros(Nrows, dtype = int)
        #kx =zeros(Nrows, dtype = int)
        ky =zeros(Nrows, dtype = int)
        kz = zeros(Nrows, dtype = int)
        ms = zeros(Nrows, dtype = int)
        for i in range(len(L)):
            pn = PQRS[L[i][0]]+L[i][2]
            kx += self.bs.states[pn,1]
            ky += self.bs.states[pn,2]
            kx += self.bs.states[pn,3]
            ms += self.bs.states[pn,4]
            #LHS += L[i][1]*(self.bs.states[pn,4] + self.bs.states[pn,2]*k + self.bs.states[pn,3]*k**2 + self.bs.states[pn,1]*k**3)
            LHS += self.bs.unique(PQRS[L[i][0]]+L[i][2])*L[i][1]
        
        #LHS = kx + ky*k + kz*k**2 + ms * k**3
        RHS = zeros(Ncols, dtype = int)
        kx =zeros(Ncols, dtype = int)
        ky =zeros(Ncols, dtype = int)
        kz = zeros(Ncols, dtype = int)
        ms = zeros(Ncols, dtype = int)
        for i in range(len(R)):
            pn = PQRS[R[i][0]]+R[i][2]
            kx += self.bs.states[pn,1]
            ky += self.bs.states[pn,2]
            kx += self.bs.states[pn,3]
            ms += self.bs.states[pn,4]
            #RHS += R[i][1]*(self.bs.states[pn,4] + self.bs.states[pn,2]*k + self.bs.states[pn,3]*k**2 + self.bs.states[pn,1]*k**3)
            RHS += self.bs.unique((PQRS[R[i][0]]+R[i][2]))*R[i][1]
        #RHS = kx + ky*k + kz*k**2 + ms * k**3
        
        
        
        
        """
        LHS = zeros(Nrows, dtype = int)
        for i in range(len(L)):
            LHS += self.bs.unique(PQRS[L[i][0]]+L[i][2])*L[i][1]

        RHS = zeros(Ncols, dtype = int)
        for i in range(len(R)):
            RHS += self.bs.unique((PQRS[R[i][0]]+R[i][2]))*R[i][1]
        
        #now, we will need to find blocks where RHS==LHS, and sort them according to L and R
        #We begin by finding the intersection of RHS and LHS
        uniques = intersect1d(LHS, RHS)

        nblocks = len(uniques) #number of blocks in configuration
        self.blocklengths.append(nblocks)
        
        self.configurations.append(array(uniques, dtype = int)) #used for inter-tensor contraction mapping
        
        self.blocks.append([]) #list for block storage
        
        U = len(self.blocks)-1 #current block configuration
        

        self.ordering.append([L, R])
        #traverse uniques, consolidate blocks
        for ni in range(len(uniques)):
            u = uniques[ni]
            row_indices = LHS==u
            col_indices = RHS==u
            row = rows[row_indices] #rows and columns in the transformed matrix
            col = cols[col_indices]  
            self.blocks[U].append([row, col])

    def map_regions_pppp(self, L = [[0,1, 14],[1,1, 14]], R=[[2,1,14],[3,1,14]]):
        
        ###########################################################
        ##
        ##  Set up configuration specified by L, R, without any storage (to later generate elements on the fly)
        ##
        ###########################################################
        
        #(1) label rows and columns
        Ns = [self.Np, self.Nq, self.Nr, self.Ns]
        #(2) determine dimensions and indices of transformed matrix
        
        #Nrows = self.Np**2
        pq = arange(self.Np**2)
        a = pq%self.Np
        b = pq//self.Np
        
        LHS = self.bs.states[a,1] + self.bs.states[b,1]
        LHS += (self.bs.states[a,2] + self.bs.states[b,2])*self.bs.k_step
        LHS += (self.bs.states[a,3] + self.bs.states[b,3])*self.bs.k_step**2
        LHS += (self.bs.states[a,4] + self.bs.states[b,4])*self.bs.k_step**3
        
        uniques = unique(LHS)
        
        #now, we will need to find blocks where RHS==LHS, and sort them according to L and R
        #We begin by finding the intersection of RHS and LHS
        #uniques = intersect1d(LHS, RHS)

        nblocks = len(uniques) #number of blocks in configuration
        self.blocklengths.append(nblocks)
        
        self.configurations.append(array(uniques, dtype = int)) #used for inter-tensor contraction mapping
        
        self.blocks.append([]) #list for block storage
        
        U = len(self.blocks)-1 #current block configuration
        

        self.ordering.append([L, R])
        #traverse uniques, consolidate blocks
        for ni in range(len(uniques)):
            u = uniques[ni]
            row_indices = LHS==u
            #col_indices = RHS==u
            row = pq[row_indices] #rows and columns in the transformed matrix
            #col = cols[col_indices]  
            self.blocks[U].append([row, row])        

    def genblock(self, u, i):
        #Generate block as interaction (not stored)
        a,b,c,d = 0,0,0,0
        t0 = clock()
        Ns = [self.Np, self.Nq, self.Nr, self.Ns]
        row, col = self.blocks[u][i]
        Nx, Ny = len(row), len(col)
        L, R = self.ordering[u]
        block = zeros((Nx, Ny), dtype = float)
        PQRS = [0,0,0,0]
        #print "::::::::::::::::::", clock()-t0
        a += clock()-t0
        for nx in range(Nx):
            for ny in range(Ny):
                t0 = clock()
                lhs = self.unpack(row[nx], L, Ns)
                rhs = self.unpack(col[ny], R, Ns)
                b += clock()-t0 
                t0 = clock()
                for i in range(len(lhs)):
                    PQRS[L[i][0]] = lhs[i] + L[i][2]
                for i in range(len(rhs)):
                    PQRS[R[i][0]] = rhs[i] + R[i][2]
                c += clock()-t0
                t0 = clock()
                p,q,r,s = PQRS

                block[nx,ny] = self.bs.v_preapproved(p,q,r,s)
                
                d += clock()-t0 #this process steals most time
        #print "::::::::::::::::::a", a
        #print "::::::::::::::::::b", b
        #print "::::::::::::::::::c", c
        self.a += a
        self.b += b
        self.c += c
        self.d += d
        return block

    def genblock_symmetric(self, u, i):
        #Generate block as interaction (not stored)
        a,b,c,d = 0,0,0,0
        t0 = clock()
        Ns = [self.Np, self.Nq, self.Nr, self.Ns]
        row, col = self.blocks[u][i]
        Nx, Ny = len(row), len(col)
        L, R = self.ordering[u]
        block = zeros((Nx, Ny), dtype = float)
        PQRS = [0,0,0,0]
        #print "::::::::::::::::::", clock()-t0
        a += clock()-t0
        for nx in range(Nx):
            for ny in range(nx, Ny):
                t0 = clock()
                lhs = self.unpack(row[nx], L, Ns)
                rhs = self.unpack(col[ny], R, Ns)
                b += clock()-t0 
                t0 = clock()
                for i in range(len(lhs)):
                    PQRS[L[i][0]] = lhs[i] + L[i][2]
                for i in range(len(rhs)):
                    PQRS[R[i][0]] = rhs[i] + R[i][2]
                c += clock()-t0
                t0 = clock()
                p,q,r,s = PQRS
                v = self.bs.v_preapproved(p,q,r,s)
                block[nx,ny] = v
                block[ny,nx] = v
                
                d += clock()-t0 #this process steals most time
        #print "::::::::::::::::::a", a
        #print "::::::::::::::::::b", b
        #print "::::::::::::::::::c", c
        self.a += a
        self.b += b
        self.c += c
        self.d += d
        return block

    def profile(self):
        print "Time spent at genblock processes:"
        print "a:", self.a
        print "b:", self.b
        print "c:", self.c
        print "d:", self.d
    
    def map(self, L = [[0,1, 14],[1,1, 14]], R=[[2,1,0],[3,1,0]]):
        
        ###########################################################
        ##
        ##  Set up and consolidate configuration specified by L, R
        ##
        ###########################################################
        
        #(1) label rows and columns
        Ns = [self.Np, self.Nq, self.Nr, self.Ns]
        
        #(2) determine dimensions and indices of transformed matrix
        Nrows = 1
        for i in range(len(L)):
            Nrows*=Ns[L[i][0]]
        rows = arange(Nrows)
        Ncols = 1
        for i in range(len(R)):
            Ncols*=Ns[R[i][0]]
        cols = arange(Ncols)
        
        
        
        #(3) Transform matrix indices to element indices in tensor
        left = self.unpack(rows, L, Ns)  ##indices in the transformed matrix, do I know in which order these occurs? (p-q-r-s) yes, from L
        right = self.unpack(cols, R, Ns)
        pqrs = []
        PQRS = [0,0,0,0]
        for i in range(len(left)):
            pqrs.append(left[i])
            PQRS[L[i][0]] = left[i]
        for i in range(len(right)):
            pqrs.append(right[i])
            PQRS[R[i][0]] = right[i]
        
        #(3) Identify blocks
        #basically ident(row) == ident(col) produces the blocks
        LHS = zeros(Nrows, dtype = int)
        for i in range(len(L)):
            LHS += self.bs.unique(PQRS[L[i][0]]+L[i][2])*L[i][1]

        RHS = zeros(Ncols, dtype = int)
        for i in range(len(R)):
            RHS += self.bs.unique((PQRS[R[i][0]]+R[i][2]))*R[i][1]
        """
        k = self.bs.k_step 
        LHS = zeros(Nrows, dtype = int)
        kx =zeros(Nrows, dtype = int)
        ky =zeros(Nrows, dtype = int)
        kz = zeros(Nrows, dtype = int)
        ms = zeros(Nrows, dtype = int)
        for i in range(len(L)):
            pn = PQRS[L[i][0]]+L[i][2]
            kx += self.bs.states[pn,1]*L[i][1]
            ky += self.bs.states[pn,2]*L[i][1]
            kx += self.bs.states[pn,3]*L[i][1]
            ms += self.bs.states[pn,4]*L[i][1]
            #LHS += L[i][1]*(self.bs.states[pn,4] + self.bs.states[pn,2]*k + self.bs.states[pn,3]*k**2 + self.bs.states[pn,1]*k**3)
            LHS += self.bs.unique(PQRS[L[i][0]]+L[i][2])*L[i][1]
        
        #LHS = kx + ky*k + kz*k**2 + ms * k**3
        RHS = zeros(Ncols, dtype = int)
        kx =zeros(Ncols, dtype = int)
        ky =zeros(Ncols, dtype = int)
        kz = zeros(Ncols, dtype = int)
        ms = zeros(Ncols, dtype = int)
        for i in range(len(R)):
            pn = PQRS[R[i][0]]+R[i][2]
            kx += self.bs.states[pn,1]*R[i][1]
            ky += self.bs.states[pn,2]*R[i][1]
            kx += self.bs.states[pn,3]*R[i][1]
            ms += self.bs.states[pn,4]*R[i][1]
            #RHS += R[i][1]*(self.bs.states[pn,4] + self.bs.states[pn,2]*k + self.bs.states[pn,3]*k**2 + self.bs.states[pn,1]*k**3)
            RHS += self.bs.unique((PQRS[R[i][0]]+R[i][2]))*R[i][1]
        """
        #RHS = kx + ky*k + kz*k**2 + ms * k**3
        #print LHS, RHS
        """  
        ivec KABx = bs.vKx.elem(A+iNh)+bs.vKx.elem(B+iNh);
        ivec KABy = iNmax*(bs.vKy.elem(A+iNh)+bs.vKy.elem(B+iNh));
        ivec KABz = iNmax2*(bs.vKz.elem(A+iNh)+bs.vKz.elem(B+iNh));
        ivec KABms = iNmax*iNmax2*(bs.vMs(A+iNh) + bs.vMs(B + iNh));
    
        ivec KAB = KABx+KABy+KABz + KABms;
        ivec KAB_unique = unique(KAB);
        """

        #now, we will need to find blocks where RHS==LHS, and sort them according to L and R
        #We begin by finding the intersection of RHS and LHS
        uniques = intersect1d(LHS, RHS)

        nblocks = len(uniques) #number of blocks in configuration
        self.blocklengths.append(nblocks)
        
        self.configurations.append(array(uniques, dtype = int)) #used for inter-tensor contraction mapping
        
        tempElements = []
        tempBlockmap = []
        self.blocks.append([]) #list for block storage

        
        U = len(self.blocks)-1 #current block configuration
        
        
        #traverse uniques, consolidate blocks
        for ni in range(len(uniques)):
            u = uniques[ni]
            row_indices = LHS==u
            col_indices = RHS==u
            row = rows[row_indices] #rows and columns in the transformed matrix
            col = cols[col_indices]  
            PQRS = [0,0,0,0]
            Nx = len(row)
            Ny = len(col)
            block = zeros((Nx,Ny), dtype = int)
            for nx in range(Nx):
                for ny in range(Ny):
                    lhs = self.unpack(row[nx], L, Ns)
                    rhs = self.unpack(col[ny], R, Ns)
                    for i in range(len(lhs)):
                        PQRS[L[i][0]] = lhs[i]
                    for i in range(len(rhs)):
                        PQRS[R[i][0]] = rhs[i]
                    p,q,r,s = PQRS
                    index = self.to(p,q,r,s)
                    tempElements.append(index)
                    tempBlockmap.append([ni, nx, ny]) #the corresponding pointer back to the element in the current block
                    block[nx,ny] = index
            self.blocks[U].append(block)

        #######
        ##
        ##   CONSOLIDATE ELEMENTS (link doubly occuring indices)
        ##
        #########
        #(1) sort temp elements and temp blockmap
        n = argsort(tempElements)
        tempElements = array(tempElements, dtype = float)[n]
        tempBlockmap = array(tempBlockmap, dtype = int)[n]
        
        #traverse sorted arrays simultaneously and compare elements, for doubly occuring elements; change pointers in block so they point to the same adress
        tempN = 0
        trueN = 0
        tempL = len(tempElements)
        trueL = len(self.elements)
        all_resolved = False
        while trueN < trueL:
            if self.elements[trueN] == tempElements[tempN]:
                #identical indices found, resolve doubly occuring indices by map to one unique index, update pointer in block
                block_n, nx, ny = tempBlockmap[tempN]
                self.blocks[U][block_n][nx,ny] = trueN #self.elements[trueN]
                trueN += 1
                tempN += 1
            else:
                if self.elements[trueN] < tempElements[tempN]:
                    trueN += 1
                else:
                    tempN += 1
                    if tempN >= tempL:
                        all_resolved = True #all elements accounted for in already existing self.elements
                        break

                    
        #append the remaining elements to self.elements
        if not all_resolved:
            #resolve remaining elements (unique in current config)
            tempRemaining = zeros(tempL-tempN, dtype = int)
            tN = 0
            while tempN<tempL:
                block_n, nx, ny = tempBlockmap[tempN]
                tempRemaining[tN] = tempElements[tempN]
                self.blocks[U][block_n][nx,ny] = trueN + tN #self.elements[trueN]
                tempN += 1
                tN += 1
            self.elements = append(self.elements, tempRemaining)
        
    #################################
    ##
    ##  Functions that communicate with other tensors
    ##
    #################################

    def matchblock(self, u, identifier):
        return where(self.configurations[u]==identifier)[0][0]
    def matchconfig(self, u, config):       
        unique_c = intersect1d(config, self.configurations[u]) #unique quantum numbers 
        pattern1 = zeros(len(unique_c), dtype = int) #pattern stores the index of the block where ident(quantum numbers) = config[i]
        pattern2 = zeros(len(unique_c), dtype = int) #pattern stores the index of the block where ident(quantum numbers) = config[i]
        for i in range(len(unique_c)):
            pattern1[i] = self.matchblock(u, unique_c[i])
            pattern2[i] = where(config==unique_c[i])[0][0]
        return pattern1, pattern2

    
    def getblock(self, u, i):
        #Get prestored block as amplitude
        #u = config
        #i = block number
        return self.elements[self.blocks[u][i]]
    
    def getraw(self, u, i):
        return self.blocks[u][i]
        

        
    def setblock(self, u, i, block):
        self.elements[self.blocks[u][i]] = block #this acutally works beautifully!!! :D
        
    def addblock(self, u, i, block):
        self.elements[self.blocks[u][i]] += block #this acutally works beautifully!!! :D

 

class CCD_block():
    def __init__(self, bs, Nh):
        #set up all diagrams needed as channelmap objects
        self.bs = bs
        self.Np = bs.nstates-Nh
        self.channels = []
        
        #initializing t2
        self.t2 = blockmap(bs, Np, Np, Nh, Nh)
        self.t2.map([[0, 1, Nh],[1,1,Nh]], [[2,1,0],[3,1,0]])
        #self.t2.map([[1, 1, Nh],[0,1,Nh]], [[2,1,0],[3,1,0]])
        self.t2.init_as_t2amp()
        print self.t2.elements
        
        self.t2new = blockmap(bs, Np, Np, Nh, Nh) #next iteration of amplitude
        self.t2new.map([[0, 1, Nh],[1,1,Nh]], [[2,1,0],[3,1,0]])
        self.t2new.init_as_t2amp()
        #print self.t2new.elements
        
        #initializing interactions
        self.vhhpp = blockmap(bs, Nh, Nh, Np, Np)
        self.vhhpp.map([[0, 1, 0],[1,1,0]], [[2,1,Nh],[3,1,Nh]])
        self.vhhpp.init()
        #self.vhhpp.map_regions([[0, 1, 0],[1,1,0]], [[2,1,Nh],[3,1,Nh]])
        
        self.vpphh = blockmap(bs, Np, Np, Nh, Nh)
        self.vpphh.map_regions([[0, 1, Nh],[1,1,Nh]], [[2,1,0],[3,1,0]])
        
        #vhhpp = blockmap(bs, Nh, Nh, Np, Np)
        
        self.vhhhh = blockmap(bs, Nh, Nh, Nh, Nh)
        self.vhhhh.map_regions([[0, 1, 0],[1,1,0]], [[2,1,0],[3,1,0]])
        
        vhpph = blockmap(bs, Nh, Np, Np, Nh)
        
        
        #self.vpppp = blockmap(bs, Np, Np, Np, Np)
        #self.vpppp.map_regions([[0, 1, Nh],[1,1,Nh]], [[2,1,Nh],[3,1,Nh]])
        #print "Number of elements in vpppp:", len(self.vpppp.elements)
        #print "Number of nonzero elements in vpppp:", len(self.vpppp.elements[self.vpppp.elements != 0])
        
        self.vpppp = laddermap(bs, Np)
        self.vpppp.map()
        
        self.align_diagrams()
        print "initialization energy:", self.energy()
        
        #self.compare()
        #c1,c2 = M1.matchconfig(u1, M2.configurations[u2])
        #print self.t2new.elements.max()
    def compare(self):
        print "1:", self.vpppp.blocklengths
        print "2:", self.vpppp2.blocklengths
        n = self.vpppp.blocklengths[0]
        for i in range(n):
            print self.vpppp.genblock(0,i)
            print self.vpppp2.genblock(0,i)
            print "...."
            
        
    def align_diagrams(self):
        #generate all channels 
        self.c1,self.c2 = self.t2.matchconfig(0, self.vhhpp.configurations[0])
        self.c3,self.c4 = self.t2.matchconfig(0, self.vpppp.configurations[0])
    def advance(self):
        #perform one iteration
        t0 = clock()
        #self.t2new.zeros()
        #for i in range(self.vpphh.blocklengths[0]):
        #    self.t2new.setblock(0, i, self.vpphh.genblock(0,i))
        self.t2new.reset_amplitude()
        print "(1):", clock()-t0
        t0 = clock()
        self.contract(self.t2, 0, self.vpppp, 0, 0, self.t2new, 0, .5, self.c3, self.c4) #Adding in L1
        print "(2):", clock()-t0
        t0 = clock()
        
        #######################
        # Updating amplitude  #
        #######################
        self.t2.zeros()
        for i in range(self.t2new.blocklengths[0]):
            self.t2.setblock(0, i, self.t2new.getblock(0,i))     
        self.t2.energy_division()   
        print "(3):", clock()-t0
        t0 = clock()
        print "energy:", self.energy()
        print "(4):", clock()-t0
        t0 = clock()
        

    def energy(self):
        e_corr = 0
        c1,c2 = self.c1, self.c2
        for i in range(len(c1)):
            block = dot(self.vhhpp.getblock(0,c2[i]),self.t2.getblock(0,c1[i]))
            e_corr += .25*sum(block.diagonal())
        return e_corr
    
    def contract_symmetric(self, M1, u1, M2, u2, channel, MR, ur, factor = 1, c1 = None, c2 = None):
        #contract M1*M2 = MR guided by channel
        if c1 == None or c2 == None:
            c1,c2 = M1.matchconfig(u1, M2.configurations[u2])
        #print MR.elements
        a,b,c,d = 0,0,0,0
        for i in range(len(c1)):
            t0 = clock()
            b1 = M2.genblock_symmetric(u2,c2[i])
            a+= clock()-t0
            t0 = clock()
            b2 = M1.getblock(u1,c1[i])
            b+= clock()-t0
            t0 = clock()
            
            block = factor * dot(b1,b2)
            c+= clock()-t0
            t0 = clock()
            
            MR.addblock(ur, c1[i], block)
            d+= clock()-t0        
        
    def contract(self, M1, u1, M2, u2, channel, MR, ur, factor = 1, c1 = None, c2 = None):
        #contract M1*M2 = MR guided by channel
        if c1 == None or c2 == None:
            c1,c2 = M1.matchconfig(u1, M2.configurations[u2])
        #print MR.elements
        a,b,c,d = 0,0,0,0
        for i in range(len(c1)):
            t0 = clock()
            b1 = M2.genblock(u2,c2[i])
            #print b1.shape
            a+= clock()-t0
            t0 = clock()
            b2 = M1.getblock(u1,c1[i])
            b+= clock()-t0
            t0 = clock()
            
            #print b1.shape, b2.shape
            
            
            block = factor * dot(b1,b2)
            c+= clock()-t0
            t0 = clock()
            
            MR.addblock(ur, c1[i], block)
            d+= clock()-t0
        #print ":::::a", a
        #print ":::::b", b
        #print ":::::c", c
        #print ":::::d", d
            #MR.setblock(ur, c1[i], block) #where to broadcast? ur, i
            
            
class laddermap():
    def __init__(self, basis, Np, Nh = 14):
        self.bs = basis
        #print "Number of states:", self.bs.nstates
        self.k_step = 2*(self.bs.Nm + 1) #steplength of momentum vector, used to uniquely identify combinations of quantum numbers
        self.elements = [] #1D array where each element is stored, common to all blocks
        self.elements = array([], dtype = float)
        self.blockmap = [] #1D array where for each element, the mapped blocks and indices are stored, so that 
        self.blocks = [] #nested list containing blocks of pointers to elements
        self.configurations = [] #nested list containing the incoming and outgoing quantum numbers of blocks
        self.blocklengths = []
        self.ordering = []
        self.Np = Np
        self.Nh = Nh
        self.stresstest = [0,0,0,0]
    def map(self):
        ###########################################################
        ##
        ##  Set up configuration specified by L, R, without any storage (to later generate elements on the fly)
        ##
        ###########################################################
        
        #(1) label rows and columns
        #Ns = [self.Np, self.Nq, self.Nr, self.Ns]
        #(2) determine dimensions and indices of transformed matrix
        Nh = self.Nh
        #Nrows = self.Np**2
        pq = arange(self.Np**2)
        a = pq%self.Np
        b = pq//self.Np
        
        #LHS = self.bs.unique(a) + self.bs.unique(b)
        
        
        LHS = self.bs.states[a+Nh,1] + self.bs.states[b+Nh,1]
        LHS += (self.bs.states[a+Nh,2] + self.bs.states[b+Nh,2])*self.bs.k_step
        LHS += (self.bs.states[a+Nh,3] + self.bs.states[b+Nh,3])*self.bs.k_step**2
        LHS += (self.bs.states[a+Nh,4] + self.bs.states[b+Nh,4])*self.bs.k_step**3
        
        uniques = unique(LHS)
        
        #now, we will need to find blocks where RHS==LHS, and sort them according to L and R
        #We begin by finding the intersection of RHS and LHS
        #uniques = intersect1d(LHS, RHS)

        nblocks = len(uniques) #number of blocks in configuration
        self.blocklengths.append(nblocks)
        
        self.configurations.append(array(uniques, dtype = int)) #used for inter-tensor contraction mapping
        
        self.blocks.append([]) #list for block storage
        
        U = len(self.blocks)-1 #current block configuration
        

        #self.ordering.append([L, R])
        #traverse uniques, consolidate blocks
        for ni in range(len(uniques)):
            u = uniques[ni]
            #row_indices = LHS==u
            row_a = a[LHS==u] #rows and columns in the transformed matrix
            row_b = b[LHS==u] #rows and columns in the transformed matrix
            print len(row_a)
            #col = cols[col_indices]  
            self.blocks[U].append([row_a, row_b]) 
    
    def genblock(self, u,i):
        t0 = clock()
        row_a = self.blocks[u][i][0] + self.Nh
        row_b = self.blocks[u][i][1] + self.Nh
        self.stresstest[0] += clock()-t0
        t0 = clock()
        n = len(row_a)
        o = ones(n, dtype = int)
        #expand
        self.stresstest[1] += clock()-t0
        t0 = clock()
        block = zeros((n,n) ,dtype = float)
        
        row_a_ = kron(o, row_a)
        row_b_ = kron(o, row_b)
        row_c = kron(row_a,o)
        row_d = kron(row_b,o)
        self.stresstest[2] += clock()-t0
        t0 = clock()
        block = self.bs.V2(self.bs.states[row_a_,:],self.bs.states[row_b_],self.bs.states[row_c],self.bs.states[row_d])
        self.stresstest[3] += clock()-t0
        t0 = clock()
        block.shape = (n,n)
        """
        for nx in range(n):
            a,b = row_a[nx], row_b[nx]
            block[nx,nx] = self.bs.v(a,b,a,b)
            for ny in range(nx+1, n):
                c,d = row_a[ny], row_b[ny]
                v = self.bs.v(a,b,c,d)
                block[nx,ny] = v
                block[ny,nx] = v
        
        """
        return block
    def profile(self):
        print self.stresstest

def ppp(eBs):
    #triple particle row
    Np = eBs.nstates-eBs.nparticles
    Nh = eBs.nparticles
    ndim = Np*(Np+1)*(Np+2)/6
    #print ndim
    #pq = zeros(ndim, dtype = int)
    p = zeros(ndim, dtype = int)
    q = zeros(ndim, dtype = int)
    r = zeros(ndim, dtype = int)
    
    count = 0
    for np in range(Np):
        for nq in range(np+1):
            for nr in range(nq+1):
                p[count] = np
                q[count] = nq
                r[count] = nr
                count +=1
    U = bs.unique(p+Nh) + bs.unique(q+Nh) + bs.unique(r+Nh)
    #print count, ndim
    return U,p,q,r

def hhh(eBs):
    #triple particle row
    Np = eBs.nstates-eBs.nparticles
    Nh = eBs.nparticles
    ndim = Nh*(Nh+1)*(Nh+2)/6
    #print ndim
    #pq = zeros(ndim, dtype = int)
    p = zeros(ndim, dtype = int)
    q = zeros(ndim, dtype = int)
    r = zeros(ndim, dtype = int)
    
    count = 0
    for np in range(Nh):
        for nq in range(np+1):
            for nr in range(nq+1):
                p[count] = np
                q[count] = nq
                r[count] = nr
                count +=1
    U = bs.unique(p) + bs.unique(q) + bs.unique(r)
    #print count, ndim
    return U,p,q,r

def pp(eBs):
    Np = eBs.nstates-eBs.nparticles
    Nh = eBs.nparticles
    ndim = Np*(Np+1)/2
    #print ndim
    #pq = zeros(ndim, dtype = int)
    p = zeros(ndim, dtype = int)
    q = zeros(ndim, dtype = int)
    
    count = 0
    for np in range(Np):
        for nq in range(np+1):
            p[count] = np
            q[count] = nq
            count +=1
    U = bs.unique(p+Nh) + bs.unique(q+Nh)
    
    return U,p,q
    
def ph(eBs):
    Nh = eBs.nstates-eBs.nparticles
    Np = eBs.nparticles
    p = array(kron(arange(Np),ones(Nh)), dtype = int)
    q = array(kron(ones(Np), arange(Nh)), dtype = int)
    U = bs.unique(p+Nh) + bs.unique(q)
    return U, p, q

def hp(eBs):
    Nh = eBs.nstates-eBs.nparticles
    Np = eBs.nparticles
    p = array(kron(arange(Nh),ones(Np)), dtype = int)
    q = array(kron(ones(Nh), arange(Np)), dtype = int)
    U = bs.unique(p) + bs.unique(q+Nh)
    return U, p, q
    
def hh(eBs):
    Nh = eBs.nstates-eBs.nparticles
    Np = eBs.nparticles
    p = array(kron(arange(Nh),ones(Nh)), dtype = int)
    q = array(kron(ones(Nh), arange(Nh)), dtype = int)
    U = bs.unique(p) + bs.unique(q)
    return U, p, q

def collect_symmetric_ppp(Np, pn,qn,rn,Kpq, K):
    Np2 = Np*Np
    k_sorted = argsort(Kpq)
    #print Kpq
    #collect indices in Kpq where Kpq == K(i)
    nj = 0
    i = 0
    C = K[i]
    
    #align k_sorted(i ) with K(j)
    while Kpq[k_sorted[nj]] < C:
        nj +=1
    
    collect = False
    blocks = []
    row = zeros(100000, dtype = int)
    nx = 0
    
    while nj < len(k_sorted):
        n = k_sorted[nj]
        #lhs_i = Kpq[n]
        #print nj, i, len(K)
        if Kpq[n] == C:
            p = pn[n]
            q = qn[n]
            r = rn[n]
            #collect
            collect = True
            row[nx] = p + q*Np + r*Np2
            nx += 1
            if r!=q:
                row[nx] = p + r*Np + q*Np2
                nx += 1
                if p!= q:
                    #none three equal - all possibilities
                    #pqr prq qpr qrp rpq rqp                     
                    row[nx] = q + p*Np + r*Np2
                    nx += 1
                    row[nx] = q + r*Np + p*Np2
                    nx += 1
                    #it follows that p!=r, so
                    row[nx] = r + p*Np + q*Np2
                    nx += 1
                    row[nx] = r + q*Np + p*Np2
                    nx += 1
            else:
                #it follows that r==q
                if p!=q:
                    row[nx] = q + p*Np + r*Np2
                    nx += 1
                if p!=r:
                    row[nx] = r + q*Np + p*Np2
                    nx += 1        
        else:
            if(collect):
                #assemble block
                collect = False
                
                #print row[0:nx]
                blocks.append(sort(row[0:nx]))
                nx = 0
                i += 1
                nj -= 1
                try:
                    C = K[i]
                except:
                    pass
        nj += 1
    blocks.append(sort(row[0:nx])) 
    return blocks

def collect_symmetric(Np, pn,qn,Kpq, K):
    k_sorted = argsort(Kpq)
    #print Kpq
    #collect indices in Kpq where Kpq == K(i)
    nj = 0
    i = 0
    C = K[i]
    
    #align k_sorted(i ) with K(j)
    while Kpq[k_sorted[nj]] < C:
        nj +=1
    
    collect = False
    blocks = []
    row = zeros(10000, dtype = int)
    nx = 0
    
    while nj < len(k_sorted):
        n = k_sorted[nj]
        #lhs_i = Kpq[n]
        #print nj, i, len(K)
        if Kpq[n] == C:
            p = pn[n]
            q = qn[n]
            
            #collect
            collect = True
            row[nx] = p + q*Np
            nx += 1
            if p!=q:
                row[nx] = q + p*Np
                nx += 1
        
        else:
            if(collect):
                #assemble block
                collect = False
                
                #print row[0:nx]
                blocks.append(sort(row[0:nx]))
                nx = 0
                i += 1
                nj -= 1
                try:
                    C = K[i]
                except:
                    pass
        nj += 1
    blocks.append(sort(row[0:nx])) 
    return blocks
                
def collect(Np, pn,qn,Kpq, K):
    k_sorted = argsort(Kpq)
    #print Kpq
    #collect indices in Kpq where Kpq == K(i)
    nj = 0
    i = 0
    C = K[i]
    
    #align k_sorted(i ) with K(j)
    while Kpq[k_sorted[nj]] < C:
        nj +=1
    
    collect = False
    blocks = []
    row = zeros(10000, dtype = int)
    nx = 0
    
    while nj < len(k_sorted):
        n = k_sorted[nj]
        #lhs_i = Kpq[n]
        #print nj, i, len(K)
        if Kpq[n] == C:
            p = pn[n]
            q = qn[n]
            
            #collect
            collect = True
            row[nx] = p + q*Np
            nx += 1
            #if p!=q:
            #    row[nx] = q + p*Np
            #    nx += 1
        
        else:
            if(collect):
                #assemble block
                collect = False
                
                #print row[0:nx]
                blocks.append(sort(row[0:nx]))
                nx = 0
                i += 1
                nj -= 1
                try:
                    C = K[i]
                except:
                    pass
        nj += 1
    blocks.append(sort(row[0:nx])) 
    return blocks

def php(eBs):
    Np = eBs.nstates-eBs.nparticles
    Nh = eBs.nparticles
    ndim = Np*Nh*(Np+1)/2
    #print ndim
    #pq = zeros(ndim, dtype = int)
    p = zeros(ndim, dtype = int)
    q = zeros(ndim, dtype = int)
    r = zeros(ndim, dtype = int)
    count = 0
    for np in range(Np):
        for nq in range(Nh):
            for nr in range(np+1):
                p[count] = np
                q[count] = nq
                r[count] = nr
                count +=1
    U = bs.unique(p+Nh) + bs.unique(q) + bs.unique(r+Nh)
    return U,p,q,r
    
def collect_php(Np, Nh, pn,qn,rn, Kpqr, K):
    k_sorted = argsort(Kpqr)
    #print Kpq
    #collect indices in Kpq where Kpq == K(i)
    nj = 0
    i = 0
    C = K[i]
    
    #align k_sorted(i ) with K(j)
    while Kpqr[k_sorted[nj]] < C:
        nj +=1
    
    collect = False
    blocks = []
    row = zeros(10000, dtype = int)
    nx = 0
    Nph = Np*Nh
    
    while nj < len(k_sorted):
        n = k_sorted[nj]
        #lhs_i = Kpq[n]
        #print nj, i, len(K)
        if Kpqr[n] == C:
            p = pn[n]
            q = qn[n]
            r = rn[n]
            #collect
            collect = True
            row[nx] = p + q*Np + r*Nph
            nx += 1
            if p!=r:
                row[nx] = r + q*Np + p*Nph
                nx += 1
        
        else:
            if(collect):
                #assemble block
                collect = False
                
                #print row[0:nx]
                blocks.append(sort(row[0:nx]))
                nx = 0
                i += 1
                nj -= 1
                try:
                    C = K[i]
                except:
                    pass
        nj += 1
    blocks.append(sort(row[0:nx])) 
    return blocks
    

    
    
def build_blocks(eBs):
    Np = eBs.nstates-eBs.nparticles
    Nh = eBs.nparticles
    
    Kabc, a, b, c = php(eBs)
    Kijk, i, j, k = hhh(eBs)
    K_unique = intersect1d(Kabc, Kijk)
    #print K_unique
    
    blocks_p = collect_symmetric_ppp(Np, a,b,c, Kabc, K_unique)
    blocks_h = collect_symmetric_ppp(Nh, i,j,k, Kijk, K_unique)
    ns = 0
    for i in range(len(K_unique)):
        ns += len(blocks_p[i])*len(blocks_h[i])
    print 8*64*ns*1.25e-10, "Gb required for t3 amplitude"
    """
    
    Kab,a,b = pp(eBs)
    Kij,i,j = hh(eBs)
    K = intersect1d(Kab, Kij)
    #K = unique(Kab)
    print len(K)
    Iij = collect_symmetric(Nh, i,j,Kij, K)
    Iab = collect_symmetric(Np, a,b,Kab, K)
    print len(Iab)
    print len(Iij)
    #print Iij
    
    #testing ppp
    Kab,a,b,c = ppp(eBs)
    
    K = unique(Kab)
    
    Iab = collect_symmetric_ppp(Np, a,b,c,Kab, K) #this really works!!!! :) :D
    print Iab[-2:]
    """
    
    """
    ABC = array(linspace(0,Np**3-1, Np**3), dtype = int)
    C = array(floor(ABC/(Np**2)), dtype =int)
    B = array(floor((ABC-C*Np**2)/Np), dtype =int)
    A = array(floor(ABC-C*Np**2-B*Np), dtype =int)
    print A,B, C
    
    
    
    #A = array(AB%Np, dtype = int)
    #B = array(AB//Np, dtype = int)
    LHS = eBs.unique(A+Nh) + eBs.unique(B +Nh) + eBs.unique(C+Nh)
    K_unique = unique(LHS)
    print len(K_unique)
    blocks = []
    for i in K_unique:
        blocks.append(ABC[LHS==i])
    print blocks[-2:]
    """
    
    
    #for comparison
    """
    AB = array(linspace(0,Np**2-1, Np**2), dtype = int)
    A = array(AB%Np, dtype = int)
    B = array(AB//Np, dtype = int)
    LHS = eBs.unique(A+Nh) + eBs.unique(B +Nh)
    K_unique = unique(LHS)
    print len(K_unique)
    blocks = []
    for i in K_unique:
        blocks.append(AB[LHS==i])
    print blocks[-2:]
    """
    
    
    
                         
        
def boxdraw(ny,nx, rx,ry):
    plot([nx-rx,nx + rx], [ny -ry, ny - ry], color = (0,0,0))
    plot([nx-rx,nx + rx], [ny +ry, ny + ry], color = (0,0,0))
    
    plot([nx+rx,nx + rx], [ny -ry, ny + ry], color = (0,0,0))
    plot([nx-rx,nx - rx], [ny -ry, ny + ry], color = (0,0,0))
    
"""
Nh = 14
Nshells =7
bs = electronbasis(Nshells, 1.0, Nh)
Np = bs.nstates-Nh
Ns = bs.nstates
print "Number of states:", bs.nstates
"""

#ccd = CCD_block(bs, Nh)
#print pp(bs)
#print ph(bs)
#print hp(bs)
#print ppp(bs)

#build_blocks(bs)


eBs = electronbasis(5, 1.0, 14)
print eBs.nstates
for i in range(eBs.nstates):
    for e in eBs.states[i]:
        print e, "&",
    print "\\\\"


if False:
    #calculate density
    NN = 17
    Nh = 14
    T = zeros(NN-3)
    NS = zeros(NN-3)
    
    for i in range(3,NN):
        Nshells =i
        bs = electronbasis(Nshells, 1.0, Nh)
        Np = bs.nstates-Nh
        Ns = bs.nstates
        print "Number of states:", bs.nstates
        
        pp = arange(Np**2)
        a = pp%Np
        b = pp//Np
        #print a, b 
        K = bs.unique(a + Nh) + bs.unique(b+Nh)
        K_unique = unique(K)
        for e in range(len(K_unique)):
            T[i-3] += len(pp[K==K_unique[e]])**2
        T[i-3] = T[i-3]/float(Np**4)
        NS[i-3] = Np
    print T
    
    plot(NS,T*100)
    title("Density for the pp-pp interaction matrix")
    xlabel("Number of particle states [$N_s$]")
    ylabel("Density [% of nonzero elements]")
    show()



"""
ab = linspace(0,Np**2-1, Np**2)

a = array(ab%Np, dtype = int)
b = array(ab//Np, dtype = int)
LHS = bs.unique(a) + bs.unique(b)
K = unique(LHS)
rows = []
for i in K:
    rows.append(ab[LHS==i])
    print rows[-1]
    
#print rows

c = []
dim = Np*(Np+1)*(Np+2)/6 #*((Np+2)/3)
AB = zeros(dim, dtype =int)
A = zeros(dim, dtype =int)
B = zeros(dim, dtype =int)
C = zeros(dim, dtype =int)

cn = 0
print Np, (Np+1)/2
for aa in range(Np):
    for bb in range(aa+1):
        #cn += 1
        for cc in range(bb+1):
            
            val = bs.unique(aa)+bs.unique(bb)
                
                
            AB[cn] = val
            A[cn] = aa
            #B[cn] = bb
            #C[cn] = cc
            #c.append(bs.unique(aa)+bs.unique(bb))
                
            cn += 1
        #print cn, dim

print cn, dim
"""


"""
K = unique(AB)

iAB = argsort(AB)
l_c = AB[iAB[0]]
"""
    
    
    
#print A+B
#print B
#scan for blocks
i = 0
nx = 0
ni = 0
    


#for i in range(bs.nstates):
#    print bs.states[i]



"""
g = []
for p in range(Ns):
    for q in range(Ns):
        for r in range(Ns):
            for s in range(Ns):
                if bs.unique(p)+bs.unique(q) == bs.unique(r) + bs.unique(s):
                    a =  [bs.states[p,0], bs.states[q,0], bs.states[r,0], bs.states[s,0]]
                    if a not in g:
                        g.append(a)
                        print a
p = 0
q = 0
"""


#while p < Ns and q < Ns:
#    if p <= q:
#        #print bs.states[p,1:4] + bs.states[q,1:4]
#        print bs.unique(p) + bs.unique(q)
#        p += 1
#    else:
#        #print bs.states[p,1:4] + bs.states[q,1:4]
#        print bs.unique(p) + bs.unique(q)
#        q += 1
        
        


"""
figure(1)
#plot blocks indices
hold("on")
N = ccd.t2.blocklengths[0]
dx = 0
dy = 0
for i in range(N):
    block = ccd.t2.getraw(0,i)
    Nx = len(block)
    Ny = len(block[0])
    for nx in range(Nx):
        for ny in range(Ny):
            text(dx + nx-.3,dy + ny, str(block[nx,ny]), size = 7)
            boxdraw(dy + ny, dx + nx , .5, .5)
    dx += Nx
    dy += Ny
hold("off")
axis("off")
show()
"""


"""
#Plot connections in block setup
Nx = 0
figure(1)
hold("on")
dx = 0
dy = 10
s = 2	
N = ccd.t2.blocklengths[0]
for i in range(ccd.t2.blocklengths[0]):
    block = ccd.t2.getraw(0,i)
    Nx = len(block)
    Ny = len(block[0])
    
    
    for nx in range(Nx):
        for ny in range(Ny):
            plot([dx + .01*nx, block[nx,ny]], [dy + .01*ny, 0], color = (i/float(N),i/float(N), .5))
    
    boxdraw(dy, dx, s*.1*Nx, s*.01*Ny)        
    dx += 3.7*Nx
boxdraw(0, dx/2, dx/2, .05)
#dy += .1*Ny
text(550, 10.5, "Blocks")
text(550, -.5, "Elements")
#text(0,0,"test")
axis("off")
show()
"""




"""
#plot blockmap dependencies
print len(ccd.t2.elements)
Nx = 0
for i in range(ccd.vpppp.blocklengths[0]):
    Nx += len(ccd.vpppp.genblock(0,i))

Z = zeros((Nx,Nx))
nx = 0
for i in range(ccd.vpppp.blocklengths[0]):
    block = ccd.vpppp.genblock(0,i)
    dx = len(block)
    #Z[nx:nx+dx, nx:nx+dx] = block
    Z[nx:nx+dx, nx:nx+dx] = block
    nx += dx

#imshow(Z, cmap = "gist_ncar")

imshow(Z, cmap = "RdGy")
axis('off')
show()
"""


"""
for i in range(10):
    ccd.advance()
print "Total time:", clock()-t
ccd.vpppp.profile()
"""




"""
t2 = blockmap(bs, Np, Np, Nh, Nh)
t2.map()
t2.init_as_t2amp()
vpppp = blockmap(bs, Np, Np, Np, Np)
vpppp.map_regions_pppp([[0, 1, Nh],[1,1,Nh]], [[2,1,Nh],[3,1,Nh]])
#vpppp.map_regions([[0, 1, Nh],[1,1,Nh]], [[2,1,Nh],[3,1,Nh]])
print "intersections:", len(intersect1d(vpppp.configurations[0], t2.configurations[0]))
print len(vpppp.configurations[0])
print vpppp.blocklengths
t0 = clock()
for i in range(vpppp.blocklengths[0]):
    a = vpppp.genblock_symmetric(0,i)
print clock()-t0
print vpppp.profile()



vpppp = laddermap(bs, Np)
vpppp.map()
print vpppp.blocklengths
t0 = clock()
for i in range(vpppp.blocklengths[0]):
    a = vpppp.getblock(0,i)
print clock()-t0
print vpppp.stresstest
#print vpppp.profile()

#vpppp.map_regions_pppp()
#vpppp.init_as_pppp()
#print len(vpppp.elements)
#print len(vpppp.elements[vpppp.elements == 0])
"""


if False:
    v1 = zeros((Nh**2, Np**2))
    v2 = zeros((Np**2, Nh**2))
    #manual MBPT(2) energy (-0.525588309385 for 66 states)
    psum2 = 0
    for i in range(Nh):
        for j in range(Nh):
            for a in range(Np):
                for b in range(Np):
                    v1[i + j*Nh, a+b*Np] = bs.v(i,j,a+Nh,b+Nh)
                    v2[a+b*Np,i + j*Nh] = bs.v(a+Nh,b+Nh,i,j)/(bs.states[i,0] + bs.states[j,0] - bs.states[a + Nh, 0] - bs.states[b+Nh,0])
    psum = .25*sum(dot(v1,v2).diagonal())
    print "Standard calc:", psum
"""
a: 0.00547195664302
b: 16.6360407147
c: 3.49098722157
d: 47.7523773696
"""


