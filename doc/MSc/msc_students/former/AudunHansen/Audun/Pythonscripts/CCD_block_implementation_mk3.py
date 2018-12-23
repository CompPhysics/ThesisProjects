from numpy import *
from time import *
from matplotlib.pyplot import *
from scipy.sparse import csr_matrix, coo_matrix

#Main goal for this implementation: avoid poor design choices
class electronbasis_pol():
    def __init__(self, N, rs, Nparticles):
        self.rs = rs
        self.states = []
        self.nstates = 0
        self.nparticles = Nparticles #particle pairs
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
                            self.nstates += 1
                            self.states.append([e, x,y,z])
                            #self.states.append([e, x,y,z,-1])
                            
            if is_shell:
                n += 1
            ene_integer += 1
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
                    term1 = self.L2/(pi*self.absdiff2(r,p))
            if self.kdspin(p,s)*self.kdspin(q,r)==1:
                if self.kdwave(p,s) != 1.0:
                    term2 = self.L2/(pi*self.absdiff2(s,p))
        return val*(term1-term2)

    def direct(self,p,q,r,s):
        #Two body interaction
        #To optimize bottleneck: vectorize this function ! (remove if-tests)
        val = 0
        terms = 0.0
        term1 = 0.0
        term2 = 0.0
        kdpl = self.kdplus(p,q,r,s)
        if kdpl != 0:
            val = 1.0/self.L3           
            #if self.kdspin(p,r)*self.kdspin(q,s)==1:
            if self.kdwave(p,r) != 1.0:
                term1 = self.L2/(pi*self.absdiff2(r,p))
        
        return val*term1
    def exchange(self,p,q,r,s):
        #Two body interaction
        #To optimize bottleneck: vectorize this function ! (remove if-tests)
        val = 0
        terms = 0.0
        term1 = 0.0
        term2 = 0.0
        kdpl = self.kdplus(p,q,r,s)
        if kdpl != 0:
            val = 1.0/self.L3
            #if self.kdspin(p,s)*self.kdspin(q,r)==1:
            if self.kdwave(p,s) != 1.0:
                #term2=(4*self.absdiff2(s,p)*pi**2)/self.L2
                #terms -= 1.0/term2
                term2 = self.L2/(pi*self.absdiff2(s,p))
        return val*term2
        
    
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
                            self.states.append([e, x,y,z, 1])
                            self.states.append([e, x,y,z,-1])
                            
            if is_shell:
                n += 1
            ene_integer += 1
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

    def direct(self,p,q,r,s):
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
            
        return val*(term1-term2)
    def exchange(self,p,q,r,s):
        #Two body interaction
        #To optimize bottleneck: vectorize this function ! (remove if-tests)
        val = 0
        terms = 0.0
        term1 = 0.0
        term2 = 0.0
        kdpl = self.kdplus(p,q,r,s)
        if kdpl != 0:
            val = 1.0/self.L3
            
            
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

class channelmap():
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
        self.Np = Np
        self.Nq = Nq
        self.Nr = Nr
        self.Ns = Ns
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
    
    def map2(self, L, R):
        #map as 2-2 diagram
        #Map the conserved channels in the tensor where k_p + k_q == k_r + k_s & m_p + m_q == m_r + m_s
        #L and R is p , q , ... = L; r , s , ... = R, giving the range of the rows an columns to be evaluates 
        #The resulting rows and cols are lists containing rows and columns of dense (nonzero) blocks in the full tensor
        
        N = len(L)
        N_rows = len(L[0][0]) #Number of indices in L
        N_cols = len(R[0][0]) #Number of indices in R
        
        ident_cols = zeros(N_cols, dtype = int)
        ident_rows = zeros(N_rows, dtype = int)
        bs = self.bs
        for i in range(len(L)):
            ident_rows += bs.unique(L[i][0]+L[i][2])*L[i][1]
        for i in range(len(R)):
            ident_cols += bs.unique(R[i][0]+R[i][2])*R[i][1]    
        #ident_rows += bs.unique(L[0]) +  bs.unique(L[1]) #bs.states[L[i],1] + bs.states[L[i],2]*self.k_step + bs.states[L[i],3]*self.k_step**2 + bs.states[L[i],4]*self.k_step**3
        #ident_cols += bs.unique(R[0]) +  bs.unique(R[1]) #bs.states[R[i],1] + bs.states[R[i],2]*self.k_step + bs.states[R[i],3]*self.k_step**2 + bs.states[R[i],4]*self.k_step**3                
        
        uniques = intersect1d(ident_rows, ident_cols)
        blocks = []
        rows = []
        cols = []
        configurations = [] #store unique values for different configs
        for u in uniques:
            row_indices = ident_rows==u
            col_indices = ident_cols==u
            
            #The following is not general
            p = L[0][0][row_indices]
            q = L[1][0][row_indices]
            r = R[0][0][col_indices]
            s = R[1][0][col_indices]
            if self.amplitude == False:
                p += L[0][2]
                q += L[1][2]
                r += R[0][2]
                s += R[1][2]
            rows.append([p,q])
            cols.append([r,s])
            configurations.append(u) #LHS = RHS = u
            """
            #store indices in a temporary 1D array
            Nx = len(row_indices)
            Ny = len(col_indices)
            for nx in range(Nx):
                for ny in range(Ny):
            """     

        nblocks = len(rows)
        self.blocklengths.append(nblocks)
        if self.amplitude:
    
            self.consolidate(rows, cols,  configurations, nblocks) #this process should only be used for amplitudes, for interactions: calculate elements on the fly
        else:
            #store only blocks, not elements
            self.blocks.append([])
            for i in range(len(rows)):
                self.blocks[-1].append([rows[i], cols[i]])
            self.configurations.append(configurations)
            
            
        #print unique(self.ident_rows)
    def broadcast(self):
        #distribute all elements to the configurations
        pass
    
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
        
        
    def consolidate(self, rows, cols, configurations, nblocks):
        #consolidate the block mapping to the self.elements indexing and return a list of properly mapped blocks
        """
        This function needs to set up each block, where the block elements contain a pointer to the location in the self.elements array where the elements are kept
        The challenge here is that the block elements are likely to be already present in the element array (unless first configuration), so that we need only append new elements
        Lookup routines may be unbearably slow unless some preconcieved ordering is utilized
        #config ID = ?
        #does it work? remains to bee seen ;)
        """
        self.configurations.append(array(configurations, dtype = int))
        
        tempElements = []
        tempBlockmap = []
        self.blocks.append([])
        
        u = len(self.blocks)-1 #current block configuration
        for i in range(nblocks):
            p,q = rows[i]
            r,s = cols[i]
            Nx = len(rows[i][0])
            Ny = len(cols[i][0])
            #print p,q,r,s, Nx, Ny
            currentblock = zeros((Nx, Ny), dtype = int)
            for nx in range(Nx):
                for ny in range(Ny):
                    index = self.to(p[nx], q[nx], r[ny], s[ny])
                    #print ":", index
                    tempElements.append(index)
                    tempBlockmap.append([i, nx, ny]) #the corresponding pointer back to the element in the current block
                    currentblock[nx,ny] = index
                    #if index == 0:
                    #    print i
            self.blocks[u].append(currentblock)
        #print rows[62]
        #print cols[62]

        #(1) sort temp elements and temp blockmap
        n = argsort(tempElements)
        tempElements = array(tempElements, dtype = float)[n]
        tempBlockmap = array(tempBlockmap, dtype = int)[n]
        #print "temp", tempElements
        #print "te0", tempElements[0]
        
        #traverse sorted arrays simultaneously and compare elements, for doubly occuring elements; change pointers in block so they point to the same adress
        tempN = 0
        trueN = 0
        tempL = len(tempElements)
        trueL = len(self.elements)
        all_resolved = False
  
        while trueN < trueL:
            #print "mapping all elements"
            if self.elements[trueN] == tempElements[tempN]:
                #identical indices found, resolve doubly occuring indices by map to one unique index, update pointer in block
                block_n, nx, ny = tempBlockmap[tempN]
                self.blocks[u][block_n][nx,ny] = trueN #self.elements[trueN]
                trueN += 1
                tempN += 1
            else:
                if self.elements[trueN] < tempElements[tempN]:
                    trueN += 1
                #if self.elements[trueN] > tempElements[tempN]:
                else:
                    tempN += 1
                    if tempN >= tempL:
                        all_resolved = True #all elements accounted for in already existing self.elements
                        break

                    
        #append the remaining elements to self.elements
        if not all_resolved:
            #
            #print "resolve remaining"
            tempRemaining = zeros(tempL-tempN, dtype = int)
            tN = 0
            while tempN<tempL:
                block_n, nx, ny = tempBlockmap[tempN]
                tempRemaining[tN] = tempElements[tempN]
                self.blocks[u][block_n][nx,ny] = trueN + tN #self.elements[trueN]
                tempN += 1
                tN += 1
            self.elements = append(self.elements, tempRemaining)
        #print "elements:", self.elements.min()


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
    def getvblock(self, u, i):
        #Get block as interaction (not stored)
        p,q = self.blocks[u][i][0][0], self.blocks[u][i][0][1]
        r,s = self.blocks[u][i][1][0], self.blocks[u][i][1][1]
        Nx = len(p)
        Ny = len(r)
        block = zeros((Nx, Ny), dtype = float)
        for nx in range(Nx):
            for ny in range(Ny):
                block[nx,ny] = bs.v(p[nx], q[nx], r[ny], s[ny])
        return block
        
    
    def getblock(self, u, i):
        #Get block as amplitude
        #u = config
        #i = block number
        return self.elements[self.blocks[u][i]]
    
    def setblock(self, u, i, block):
        self.elements[self.blocks[u][i]] = block #this acutally works beautifully!!! :D
        
    def getblock2(self, i):
        #set up and return a dense matrix for block i
        nx = len(self.rows[i][0])
        ny = len(self.cols[i][0])
        block = zeros((nx,ny))
        for x in range(nx):
            for y in range(ny):
                block[x,y] = self.bs.v(self.rows[i][0][x],self.rows[i][1][x],self.cols[i][0][y],self.cols[i][1][y])
        #print "unique:", self.bs.unique(self.rows[i][0][0]) + self.bs.unique(self.rows[i][1][0])
        #print "identifier:", self.configurations[i]
        return block
    def zeros(self):
        #set all elements to zero
        self.elements *= 0
        
    def init(self):
        #as vhhpp
        for i in range(len(self.elements)):
            p,q,r,s = self.of(self.elements[i]) #i,j,a,b
            self.elements[i] = self.bs.v(p,q,r+self.Np,s+self.Np)     
    
    def init_as_t2amp(self):
        
        self.energy_denom = zeros(len(self.elements), dtype = float)
        for i in range(len(self.elements)):
            p,q,r,s = self.of(self.elements[i]) #a,b,i,j
            self.elements[i] = self.bs.v(p+self.Nr,q+self.Nr,r,s)/(bs.states[r,0] + bs.states[s,0] - bs.states[self.Nr + p,0] - bs.states[self.Nr + q,0])
            self.energy_denom[i] = bs.states[r,0] + bs.states[s,0] - bs.states[self.Nr + p,0] - bs.states[self.Nr + q,0]
        #print "energy:", self.energy_denom
        #self.elements/=self.energy_denom
    
    
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
    
    def setup_as_pq_rs(self, L = [[0,1],[1,1]], R=[[2,1],[3,1]]):
        #(1) label rows and columns
        #this is general to all setups
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
        
        
        left = self.unpack(rows, L, Ns)  ##indices in the transformed matrix, do I know in which order these occurs? (p-q-r-s) yes, from L
        right = self.unpack(cols, R, Ns)
        #print left
        pqrs = []
        PQRS = [0,0,0,0]
        for i in range(len(left)):
            pqrs.append(left[i])
            PQRS[L[i][0]] = left[i]
        for i in range(len(right)):
            pqrs.append(right[i])
            PQRS[R[i][0]] = right[i]
        
        #to sort these as they occur in the unaligned matrix, use L, R
        
        
        #p,q,r,s = PQRS
        #(3) Identify blocks
        #basically ident(row) == ident(col) produces the blocks
        LHS = pqrs[L[0][0]]*L[0][1]
        for i in range(len(L)-1):
            LHS += self.bs.unique(pqrs[L[i+1][0]]*L[i+1][1])
            
        
        RHS = pqrs[R[0][0]]*R[0][1]
        for i in range(len(R)-1):
            RHS += self.bs.unique(pqrs[R[i+1][0]]*R[i+1][1])
        
        #now, we will need to find blocks where RHS==LHS, and sort them according to L and R
        
        #We begin by finding the intersection of RHS and LHS
        uniques = intersect1d(LHS, RHS)
        
        nblocks = len(uniques)
        
        #traverse uniques
        for u in uniques:
            row_indices = LHS==u
            col_indices = RHS==u
            row = rows[row_indices] #rows and columns in the transformed matrix
            col = cols[col_indices]
            
            #left = self.unpack(rows, L, Ns)  ##indices in the transformed matrix, do I know in which order these occurs? (p-q-r-s) yes, from L
            #right = self.unpack(cols, R, Ns)
            
            
            
            
            PQRS = [0,0,0,0]
            Nx = len(row)
            Ny = len(col)
            block = zeros((Nx,Ny), dtype = int)
            for nx in range(Nx):
                for ny in range(Ny):
                    lhs = self.unpack(row[nx], L, Ns)
                    rhs = self.unpack(col[ny], R, Ns)
                    #print lhs, rhs
                    for i in range(len(lhs)):
                        PQRS[L[i][0]] = lhs[i]
                    for i in range(len(rhs)):
                        PQRS[R[i][0]] = rhs[i]
                    #print PQRS
                    p,q,r,s = PQRS
                    block[nx,ny] = self.to(p,q,r,s)
                    
            
            
            
            """
            #now, find p,q,r,s and broadcast to aligned matrix
            pqrs_ = []
            for i in range(len(L)):
                print pq[L[i][0]]
                pqrs_.append(L[i][0][row_indices])
            for i in range(len(R)):
                pqrs_.append(R[i][0][col_indices])
            p,q,r,s = pqrs_
            """
            
            
            #print pq[row_indices]
            #print rs[col_indices]
            #print " "
            #The following is not general
            #p = L[0][0][row_indices]
            #q = L[1][0][row_indices]
            #r = R[0][0][col_indices]
            #s = R[1][0][col_indices]
        
    
    def t2(self):
        bs = self.bs
        self.amplitude = True
        Np = bs.nstates-bs.nparticles
        Nh = bs.nparticles #conflicting naming should be resolved in final implementation
        
        ##############
        # find basic arrangement
        ###############
        
        row = arange(Nh**2)
        i = row % Nh
        j = row // Nh
        
        col = arange(Np**2)
        a = col % Np
        b = col // Np
        
        L = [[a,1, Nh],[b,1, Nh]]
        R = [[i,1,0],[j,1, 0]]
        self.map2(L,R)
        
        
        ##############
        # find ck-ai, as used in Q, and L3
        ###############
        row = arange(Nh*Np)
        c = row // Np
        k = row % Np
        
        col = arange(Nh*Np)
        i = col // Nh
        a = col % Nh
        
        L = [[c,1,Nh],[k,-1, 0]]
        R = [[i,1, Nh],[a,-1, Nh]]
        self.map2(L,R)
        
        
        self.init_as_t2amp()
        
    def vpphh(self):
        bs = self.bs
        Np = bs.nstates-bs.nparticles
        Nh = bs.nparticles #conflicting naming should be resolved in final implementation
        
        ##############
        # find basic arrangement
        ###############
        
        row = arange(Nh**2)
        i = row % Nh
        j = row // Nh
        
        col = arange(Np**2)
        a = col % Np
        b = col // Np
        
        L = [[a,1, Nh],[b,1, Nh]]
        R = [[i,1,0],[j,1, 0]]
        self.map2(L,R)
                
    
    def vpppp(self):
        bs = self.bs
        Np = bs.nstates-bs.nparticles
        Nh = bs.nparticles #conflicting naming should be resolved in final implementation

        col = arange(Np**2)
        a = col % Np
        b = col // Np
        
        L = [[a,1,Nh],[b,1,Nh]]
        R = [[a,1,Nh],[b,1,Nh]]
        self.map2(L,R)

    def vhhhh(self):
        bs = self.bs
        Np = bs.nstates-bs.nparticles
        Nh = bs.nparticles #conflicting naming should be resolved in final implementation

        col = arange(Nh**2)
        i = col % Nh
        j = col // Nh
        
        L = [[i,1,0],[j,1,0]]
        R = [[i,1,0],[j,1,0]]
        self.map2(L,R)

    def vhpph(self):
        bs = self.bs
        Np = bs.nstates-bs.nparticles
        Nh = bs.nparticles #conflicting naming should be resolved in final implementation

        row = arange(Nh*Np)
        a = row // Nh
        i = row % Nh
        
        col = arange(Nh*Np)
        j = col // Nh
        b = col % Nh
        
        L = [[i,1,0],[a,-1, Nh]]
        R = [[b,1, Nh],[j,-1, 0]]
        self.map2(L,R)

        
    
    def vhhpp(self):
        bs = self.bs
        Np = bs.nstates-bs.nparticles
        Nh = bs.nparticles #conflicting naming should be resolved in final implementation
        
        ##############
        # find basic arrangement (used in Q1)
        ###############
        
        row = arange(Nh**2)
        i = row % Nh
        j = row // Nh
        
        col = arange(Np**2)
        a = col % Np
        b = col // Np
        
        L = [[i,1, 0],[j,1,0]]
        R = [[a,1,Nh],[b,1, Nh]]
        self.map2(L,R)
        

        #####################
        # find vhpph (used for L3)
        ####################

        row = arange(Nh*Np)
        a = row // Nh
        i = row % Nh
        
        col = arange(Nh*Np)
        j = col // Nh
        b = col % Nh
        
        L = [[i,1,0],[a,-1, Nh]]
        R = [[b,1, Nh],[j,-1, 0]]
        self.map2(L,R)


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
    
    def init_as_t2amp(self):
        #initialize elements as t2 amplitude
        self.energy_denom = zeros(len(self.elements), dtype = float)
        for i in range(len(self.elements)):
            p,q,r,s = self.of(self.elements[i]) #a,b,i,j
            self.elements[i] = self.bs.v(p+self.Nr,q+self.Nr,r,s)/(bs.states[r,0] + bs.states[s,0] - bs.states[self.Nr + p,0] - bs.states[self.Nr + q,0])
            self.energy_denom[i] = bs.states[r,0] + bs.states[s,0] - bs.states[self.Nr + p,0] - bs.states[self.Nr + q,0]

    
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
            

    def genblock(self, u, i):
        #Generate block as interaction (not stored)
        Ns = [self.Np, self.Nq, self.Nr, self.Ns]
        row, col = self.blocks[u][i]
        Nx, Ny = len(row), len(col)
        L, R = self.ordering[u]
        block = zeros((Nx, Ny), dtype = float)
        PQRS = [0,0,0,0]
        for nx in range(Nx):
            for ny in range(Ny):
                lhs = self.unpack(row[nx], L, Ns)
                rhs = self.unpack(col[ny], R, Ns)
                for i in range(len(lhs)):
                    PQRS[L[i][0]] = lhs[i] + L[i][2]
                for i in range(len(rhs)):
                    PQRS[R[i][0]] = rhs[i] + R[i][2]
                p,q,r,s = PQRS
                block[nx,ny] = self.bs.v(p,q,r,s)
        return block
    
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
        

        
    def setblock(self, u, i, block):
        self.elements[self.blocks[u][i]] = block #this acutally works beautifully!!! :D
        


 

class CCD_block():
    def __init__(self, bs, Nh):
        #set up all diagrams needed as channelmap objects
        pass
    def align_diagrams(self):
        #generate all channels 
        pass
    def advance(self):
        #perform one iteration
        pass
    def energy(self):
        pass
    def contract(self, M1, u1, M2, u2, channel, MR, ur):
        #contract M1*M2 = MR guided by channel
        c1,c2 = channel
        for i in range(len(c2)):
            #MR.setblock( dot(M1.getblock(u, c1[i]), M2.getblock(u2, c2)) #where to broadcast? ur, i
            pass
            
             
        



def show_regions(bs):
    Z = zeros((Nh**2, Np**2), dtype = int)   
    Ns = bs.nstates
    for p in range(Nh):
        for q in range(Nh):
            for r in range(Np):
                for s in range(Np):
                    if bs.kdplus(p,q,r,s):
                        Z[p+q*Nh,r+s*Np] += 1

    imshow(Z)
    show()

Nh = 14
Nshells =5
bs = electronbasis(Nshells, 1.0, Nh)
Np = bs.nstates-Nh
print "Number of states:", bs.nstates
vpphh = channelmap(bs, Np, Np, Nh, Nh)
vhhpp = blockmap(bs, Nh, Nh, Np, Np)
vhhhh = channelmap(bs, Nh, Nh, Nh, Nh)
vhpph = channelmap(bs, Nh, Nh, Nh, Nh)
vpppp = channelmap(bs, Np, Np, Np, Np)
t2amp = channelmap(bs, Np, Np, Nh, Nh)
t2prev = channelmap(bs,Np, Np, Nh, Nh)





#vpphh.vpphh()
##vpppp.map()
#vhhhh.map()
#vhpph.map()
vhhpp.map_regions([[0, 1, 0],[1,1,0]], [[2,1,Nh],[3,1,Nh]])

#v = blockmap(bs, Nh, Nh, Np , Np)
#v.map_regions([[0, 1, 0],[1,1,0]], [[2,1,Nh],[3,1,Nh]])
#for i in range(v.blocklengths[0]):
#    print v.genblock(0,i)


t2amp.map([[0, 1, Nh],[1,1,Nh]], [[2,1,0],[3,1,0]])
t2amp.map([[1, 1, Nh],[0,1,Nh]], [[2,1,0],[3,1,0]])
#t2amp.t2()
t2amp.init_as_t2amp()
print "Elements:", t2amp.elements, len(t2amp.elements)

print "number of nonzero elements:", len(nonzero(t2amp.elements)[0])
c1,c2 = t2amp.matchconfig(0, vhhpp.configurations[0])
#print c1,c2
psum = 0

for i in range(len(c1)):
    #print " "
    #print vpppp.getvblock(0,c2[i]).shape
    #print t2amp.getblock(0,c1[i]).shape
    block = dot(vhhpp.genblock(0,c2[i]),t2amp.getblock(0,c1[i]))
    #t2prev.setblock(0,c1[i], block)
    psum += .25*sum(block.diagonal())
#print t2prev.elements
print "Block calc:", psum


B = electronbasis_pol(Nshells, 1.0, Nh)
Ns = B.nstates - 7
ec = 0
for i in range(Nh/2):
    for j in range(Nh/2):
        for a in range(Ns):
            for b in range(Ns):
                Di = B.direct(i,j,a+Nh/2,b+Nh/2)
                Ex = B.exchange(i,j,a+Nh/2,b+Nh/2)
                
                TDi = B.direct(a+Nh/2,b+Nh/2,i,j)
                TEx = B.exchange(a+Nh/2,b+Nh/2,i,j)
                
                Eabij = (B.h(a+Nh/2,a+Nh/2) + B.h(b+Nh/2,b+Nh/2) - B.h(i,i) - B.h(j,j))
                
                
                ec += (2*Di - Ex)*TDi/Eabij
                #ec += Di**2/Eabij
                
print -ec

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

if False:
    v1 = zeros(((Nh/2)**2, (Np/2)**2))
    v2 = zeros(((Np/2)**2, (Nh/2)**2))
    #manual MBPT(2) energy (-0.525588309385 for 66 states)
    psum2 = 0
    e_c = 0
    for i in range(Nh/2):
        for j in range(Nh/2):
            for a in range(Np/2):
                for b in range(Np/2):
                    v1[i + j*Nh/2, a+b*Np/2] = bs.v(2*i,2*j,2*a+Nh,2*b+Nh)
                    v2[a+b*Np/2,i + j*Nh/2] = bs.v(2*a+Nh,2*b+Nh,2*i,2*j)/(bs.states[2*i,0] + bs.states[2*j,0] - bs.states[2*a + Nh, 0] - bs.states[2*b+Nh,0])
                    e_c += (2*bs.direct(2*i,2*j,2*a+Nh,2*b+Nh)+bs.exchange(2*i,2*j,2*a+Nh,2*b+Nh))*(2*bs.direct(2*a+Nh,2*b+Nh,2*i,2*j)-bs.exchange(2*a+Nh,2*b+Nh,2*i,2*j))/(bs.states[2*i,0] + bs.states[2*j,0] - bs.states[2*a + Nh, 0] - bs.states[2*b+Nh,0])
    
    
    psum = .25*sum(dot(v1,v2).diagonal())
    print "Standard calc:", psum, .25*e_c