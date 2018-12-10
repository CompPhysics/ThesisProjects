from numpy import *
from time import clock

class amplitude():
    def __init__(self, Ns, Np):
        self.elements = []
        self.N = 0 #length of list
        self.Ns = Ns
        self.Np = Np
        self.p1sort = []
        self.p2sort = []
        self.h1sort = []
        self.h2sort = []

        self.p1V = []
        self.p2V = []
        self.h1V = []
        self.h2V = []
        
        self.p1p2 = []
        self.p1h1 = []
        self.p1h2 = []
        
        self.p2h1 = []
        self.p2h2 = []
        
        self.h1h2 = []
        for p1 in range(self.Ns-self.Np):
            self.p1p2.append([])
            self.p1h1.append([])
            self.p1h2.append([])
            self.p2h1.append([])
            self.p2h2.append([])
            
            for p2 in range(self.Ns-self.Np):
                self.p1p2[p1].append([])
            for h1 in range(self.Np):
                self.p1h1[p1].append([])
                self.p1h2[p1].append([])
                self.p2h2[p1].append([])
                self.p2h1[p1].append([])
        for h1 in range(self.Np):
            self.h1h2.append([])
            for h2 in range(self.Np):
                self.h1h2[h1].append([])
                
                
            
                
    def reset(self):
        del (self.elements)
        del (self.p1sort)
        del (self.p2sort)
        del (self.h1sort)
        del (self.h2sort)
        del (self.p1V)
        del (self.p2V)
        del (self.h1V)
        del (self.h2V)
        
        self.elements = []
        self.N = 0 #length of list
        self.p1sort = []
        self.p2sort = []
        self.h1sort = []
        self.h2sort = []

        self.p1V = []
        self.p2V = []
        self.h1V = []
        self.h2V = []
        
        self.p1V = []
        self.p2V = []
        self.h1V = []
        self.h2V = []                        
        #temporarily implementation     
    def at(self, i):
        #return the ith element
        return self.elements[i]
    def add(self, element):
        #append element to end of list
        self.elements.append(element)
        #if element[0] == self.p1N:
        #    self.p1sort.append(self.N)

        self.N += 1
    def map_tensor(self, maxindex):
        #Map out sequences of amplitudes corresponding to fixed indices
        self.p1N = 0
        self.p2N = 0
        self.h1N = 0
        self.h2N = 0
        for m in range(maxindex):
            self.p1V.append(self.p1N)
            self.p2V.append(self.p2N)
            self.h1V.append(self.h1N)
            self.h2V.append(self.h2N)
            
            for n in range(self.N):
                p1,p2,h1,h2, val = self.at(n)
                if p1 == m:
                    self.p1sort.append(n)
                    self.p1N += 1
                if p2 == m:
                    self.p2sort.append(n)
                    self.p2N += 1
                if h1 == m:
                    self.h1sort.append(n)
                    self.h1N += 1
                if h2 == m:
                    self.h2sort.append(n)
                    self.h2N += 1
                
                for m2 in range(maxindex):
                    if p1==m and p2==m2:
                        self.p1p2[m-self.Np][m2-self.Np].append(n)
                    
                    if p1==m and h1==m2:
                        self.p1h1[m-self.Np][m2].append(n)
                    
                    if p1==m and h2==m2:
                        self.p1h2[m-self.Np][m2].append(n)
                    
                    if p2==m and h1==m2:
                        self.p2h1[m-self.Np][m2].append(n)
                    
                    if p2==m and h2==m2:
                        self.p2h2[m-self.Np][m2].append(n)
                    
                    if h1==m and h2==m2:
                        self.h1h2[m][m2].append(n)
                    
                        
                
        self.p1V.append(self.p1N)
        self.p2V.append(self.p2N)
        self.h1V.append(self.h1N)
        self.h2V.append(self.h2N)
        #print len(self.h1V)
        #print len(self.h1sort)
    def getlist_p1p2(self, n1,n2):
        return self.p1p2[n1-self.Np][n2-self.Np]
    def getlist_p1h1(self, n1,n2):
        return self.p1h1[n1-self.Np][n2]
    def getlist_p1h2(self, n1,n2):
        return self.p1h2[n1-self.Np][n2]
    def getlist_p2h1(self, n1,n2):
        return self.p2h1[n1-self.Np][n2]
    def getlist_p2h2(self, n1,n2):
        return self.p2h2[n1-self.Np][n2]
    def getlist_h1h2(self, n1,n2):
        return self.h1h2[n1][n2]
    
    def getlist_p1(self,n):
        #return self.elements[self.p1sort[self.p1V[n]:self.p1V[n+1]]]
        return self.p1sort[self.p1V[n]:self.p1V[n+1]]
    def getlist_p2(self,n):
        #return self.elements[self.p1sort[self.p1V[n]:self.p1V[n+1]]]
        return self.p2sort[self.p2V[n]:self.p2V[n+1]]
    def getlist_h1(self,n):
        #return self.elements[self.p1sort[self.p1V[n]:self.p1V[n+1]]]
        #print self.h1sort
        return self.h1sort[self.h1V[n]:self.h1V[n+1]]
    def getlist_h2(self,n):
        #return self.elements[self.p1sort[self.p1V[n]:self.p1V[n+1]]]
        return self.h2sort[self.h2V[n]:self.h2V[n+1]]
        
def initialize_t2amps(Np, t2amps, basis):
    Ns = basis.nstates
    #Initialize amplitudes for Np particles and Ns states
    for i in range(Np):
        for j in range(i,Np):
            for a in range(Np,Ns):
                for b in range(a, Ns):
                    if basis.kdfullplus(a,b,i,j)==1:
                        val = basis.v(i,j,a,b)/(basis.h(i,i) + basis.h(j,j) - basis.h(a,a) - basis.h(b,b))
                        #t2amps2[a,b,i,j] = val
                        if abs(val)>.00000001:
                            t2amps.add([a,b,i,j, val])
                            #t2amps.add([a,b,j,i, -val])
                            #t2amps.add([b,a,i,j, -val])
                            #t2amps.add([b,a,j,i, val])


class electronbasis():
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

def CCD_linear_optimized(a,b,i,j, t2amps, bs):
    #optimized routine II
    l2a = 0.0
    l2b = 0.0
    l2c = 0.0
    ops = 0
    for n in t2amps.getlist_h1(i):
        p1,p2,h1,h2, t2val = t2amps.at(n)
        if h2 == j:
            l2a += bs.v(a,b,p1,p2)*t2val
            ops += 1
            
    for n in t2amps.getlist_p2(b):
        p1,p2,h1,h2, t2val = t2amps.at(n)

        if p1 ==a:
            l2b += bs.v(h1,h2,i,j)*t2val
            ops += 1


        if h1==i:
            l2c -= bs.v(a,h2,p1,j)*t2val
            ops += 1
        if h1==j:
            l2c += bs.v(a,h2,p1,i)*t2val
            ops += 1
    for n in t2amps.getlist_p2(a):
        p1,p2,h1,h2, t2val = t2amps.at(n)

        if h1==i:
            l2c += bs.v(b,h2,p1,j)*t2val
            ops += 1
    
        if h1==j:
            l2c -= bs.v(b,h2,p1,i)*t2val
            ops += 1
    #print "Number of operations in list-implementation:", ops
    return .5*l2a + .5*l2b  + l2c

def CCD_quadratic_optimized_mk2(a,b,i,j,t2amps,bs):
    #//Quadratic contributions to the CCD amplitude equation
    #//Call the elements in the following order: vmin(i,j)(a,b)
    Qa = 0
    Qb = 0
    Qc = 0
    Qd = 0
    
    for na in t2amps.getlist_h1h2(i,j):
        for nb in t2amps.getlist_p1p2(a,b):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            Qa += bs.v(h1b,h2b,p1a,p2a)*t2aval*t2bval
    
    
    #for na in t2amps.getlist_h1(i):
    #    for nb in t2amps.getlist_p1(a):
    #        p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
    #        p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
    #        if h2a == j and p2b == b:
    #            Qa += bs.v(h1b,h2b,p1a,p2a)*t2aval*t2bval
    #Substitutes          
    #for k in range(Np):
    #    for l in range(Np):
    #        for c in range(Np, Ns):
    #            for d in range(Np, Ns):
    #                Qa += v(k,l,c,d)*t2amps[c,d,i,j]*t2amps[a,b,k,l];
    
    for na in t2amps.getlist_p1h1(a,i):
        for nb in t2amps.getlist_p1h1(b,j):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            Qb += bs.v(h2a,h2b,p2a,p2b)*t2aval*t2bval
    for na in t2amps.getlist_p1h1(a,j):
        for nb in t2amps.getlist_p1h1(b,i):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            Qb -= bs.v(h2a,h2b,p2a,p2b)*t2aval*t2bval            
    
    #for na in t2amps.getlist_p1(a): #iterates over all amplitudes where <p1,p2||h1,h2> = <a,p2||h1,h2>
    #    for nb in t2amps.getlist_p1(b):
    #        p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
    #        p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
    #        if h1a == i and h1b == j:
    #            Qb += bs.v(h2a,h2b,p2a,p2b)*t2aval*t2bval
    #        if h1a == j and h1b == i:
    #            Qb -= bs.v(h2a,h2b,p2a,p2b)*t2aval*t2bval
    #Substitutes
    #for k in range(Np):
    #    for l in range(Np):
    #        for c in range(Np, Ns):
    #            for d in range(Np, Ns):
    #                Qb += v(k,l,c,d)*t2amps[a,c,i,k]*t2amps[b,d,j,l];
    #                Qb -= v(k,l,c,d)*t2amps[a,c,j,k]*t2amps[b,d,i,l];

    for na in t2amps.getlist_h1(i):
        for nb in t2amps.getlist_p1p2(a,b):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            if h2b == j:
                Qc -= bs.v(h2a,h1b,p2a,p1a)*t2aval*t2bval
    for na in t2amps.getlist_h1(j):
        for nb in t2amps.getlist_p1p2(a,b):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            if h2b == i:
                Qc -= bs.v(h2a,h1b,p2a,p1a)*t2aval*t2bval
            
                               
    #Q_c contribution
    #for na in t2amps.getlist_h1(i): #iterates over all amplitudes where <p1,p2||h1,h2> = <p1,p2||i,h2>
    #    for nb in t2amps.getlist_p1(a):
    #        p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
    #        p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
    #        if p2b == b and h2b == j:
    #            Qc -= bs.v(h2a,h1b,p2a,p1a)*t2aval*t2bval
    #for na in t2amps.getlist_h1(j): #iterates over all amplitudes where <p1,p2||h1,h2> = <p1,p2||j,h2>
    #    for nb in t2amps.getlist_p1(a):
    #        p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
    #        p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
    #        if p2b == b and h2b == i:
    #            Qc += bs.v(h2a,h1b,p2a,p1a)*t2aval*t2bval   
    #Substitutes                     
    #for k in range(Np):
    #    for l in range(Np):
    #        for c in range(Np,Ns):
    #            for d in range(Np,Ns):
    #                Qc -= v(k,l,c,d)*t2amps[d,c,i,k]*t2amps[a,b,l,j];
    #                Qc += v(k,l,c,d)*t2amps[d,c,j,k]*t2amps[a,b,l,i];

    for na in t2amps.getlist_p1(a):
        for nb in t2amps.getlist_p2h1(b,i):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            if h2b == j:
                Qd -= bs.v(h2a,h1a,p2a,p1b)*t2aval*t2bval
                
    for na in t2amps.getlist_p1(b):
        for nb in t2amps.getlist_p2h1(a,i):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            if h2b == j:
                Qd += bs.v(h2a,h1a,p2a,p1b)*t2aval*t2bval            
    
    #for na in t2amps.getlist_p1(a): #iterates over all amplitudes where <p1,p2||h1,h2> = <a,p2||h1,h2>
    #    for nb in t2amps.getlist_p2(b):
    #        p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
    #        p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
    #        if h1b == i and h2b == j:
    #            Qd -= bs.v(h2a,h1a,p2a,p1b)*t2aval*t2bval
    #for na in t2amps.getlist_p1(b): #iterates over all amplitudes where <p1,p2||h1,h2> = <b,p2||h1,h2>
    #    for nb in t2amps.getlist_p2(a):
    #        p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
    #        p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
    #        if h1b == i and h2b == j:
    #            Qd += bs.v(h2a,h1a,p2a,p1b)*t2aval*t2bval  
    #Substitutes          
    #for k in range(Np):
    #    for l in range(Np):
    #        for c in range(Np,Ns):
    #            for d in range(Np, Ns):
    #                Qd -= v(k,l,c,d)*t2amps[a,c,l,k]*t2amps[d,b,i,j];
    #                Qd += v(k,l,c,d)*t2amps[b,c,l,k]*t2amps[d,a,i,j];

    return 0.25*Qa + Qb + 0.5*(Qc + Qd)

def CCD_quadratic_optimized(a,b,i,j,t2amps,bs):
    #//Quadratic contributions to the CCD amplitude equation
    #//Call the elements in the following order: vmin(i,j)(a,b)
    Qa = 0
    Qb = 0
    Qc = 0
    Qd = 0
    
    for na in t2amps.getlist_h1(i):
        for nb in t2amps.getlist_p1(a):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            if h2a == j and p2b == b:
                Qa += bs.v(h1b,h2b,p1a,p2a)*t2aval*t2bval
    #Substitutes          
    #for k in range(Np):
    #    for l in range(Np):
    #        for c in range(Np, Ns):
    #            for d in range(Np, Ns):
    #                Qa += v(k,l,c,d)*t2amps[c,d,i,j]*t2amps[a,b,k,l];
    
    
    for na in t2amps.getlist_p1(a): #iterates over all amplitudes where <p1,p2||h1,h2> = <a,p2||h1,h2>
        for nb in t2amps.getlist_p1(b):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            if h1a == i and h1b == j:
                Qb += bs.v(h2a,h2b,p2a,p2b)*t2aval*t2bval
            if h1a == j and h1b == i:
                Qb -= bs.v(h2a,h2b,p2a,p2b)*t2aval*t2bval
    #Substitutes
    #for k in range(Np):
    #    for l in range(Np):
    #        for c in range(Np, Ns):
    #            for d in range(Np, Ns):
    #                Qb += v(k,l,c,d)*t2amps[a,c,i,k]*t2amps[b,d,j,l];
    #                Qb -= v(k,l,c,d)*t2amps[a,c,j,k]*t2amps[b,d,i,l];


    #Q_c contribution
    for na in t2amps.getlist_h1(i): #iterates over all amplitudes where <p1,p2||h1,h2> = <p1,p2||i,h2>
        for nb in t2amps.getlist_p1(a):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            if p2b == b and h2b == j:
                Qc -= bs.v(h2a,h1b,p2a,p1a)*t2aval*t2bval
    for na in t2amps.getlist_h1(j): #iterates over all amplitudes where <p1,p2||h1,h2> = <p1,p2||j,h2>
        for nb in t2amps.getlist_p1(a):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            if p2b == b and h2b == i:
                Qc += bs.v(h2a,h1b,p2a,p1a)*t2aval*t2bval   
    #Substitutes                     
    #for k in range(Np):
    #    for l in range(Np):
    #        for c in range(Np,Ns):
    #            for d in range(Np,Ns):
    #                Qc -= v(k,l,c,d)*t2amps[d,c,i,k]*t2amps[a,b,l,j];
    #                Qc += v(k,l,c,d)*t2amps[d,c,j,k]*t2amps[a,b,l,i];


    for na in t2amps.getlist_p1(a): #iterates over all amplitudes where <p1,p2||h1,h2> = <a,p2||h1,h2>
        for nb in t2amps.getlist_p2(b):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            if h1b == i and h2b == j:
                Qd -= bs.v(h2a,h1a,p2a,p1b)*t2aval*t2bval
    for na in t2amps.getlist_p1(b): #iterates over all amplitudes where <p1,p2||h1,h2> = <b,p2||h1,h2>
        for nb in t2amps.getlist_p2(a):
            p1a,p2a,h1a,h2a,t2aval = t2amps.at(na)
            p1b,p2b,h1b,h2b,t2bval = t2amps.at(nb)
            if h1b == i and h2b == j:
                Qd += bs.v(h2a,h1a,p2a,p1b)*t2aval*t2bval  
    #Substitutes          
    #for k in range(Np):
    #    for l in range(Np):
    #        for c in range(Np,Ns):
    #            for d in range(Np, Ns):
    #                Qd -= v(k,l,c,d)*t2amps[a,c,l,k]*t2amps[d,b,i,j];
    #                Qd += v(k,l,c,d)*t2amps[b,c,l,k]*t2amps[d,a,i,j];

    return 0.25*Qa + Qb + 0.5*(Qc + Qd)

def CCD_quadratic_standard(a,b,i,j,t2amps,bs):
    Np = bs.nparticles
    Ns = bs.nstates
    #//Quadratic contributions to the CCD amplitude equation
    #//Call the elements in the following order: vmin(i,j)(a,b)
    Qa = 0
    Qb = 0
    Qc = 0
    Qd = 0

    for k in range(Np):
        for l in range(Np):
            for c in range(Np, Ns):
                for d in range(Np, Ns):
                    Qa += bs.v(k,l,c,d)*t2amps[c,d,i,j]*t2amps[a,b,k,l];


    for k in range(Np):
        for l in range(Np):
            for c in range(Np, Ns):
                for d in range(Np, Ns):
                    Qb += bs.v(k,l,c,d)*t2amps[a,c,i,k]*t2amps[b,d,j,l];
                    Qb -= bs.v(k,l,c,d)*t2amps[a,c,j,k]*t2amps[b,d,i,l];


    for k in range(Np):
        for l in range(Np):
            for c in range(Np,Ns):
                for d in range(Np,Ns):
                    Qc -= bs.v(k,l,c,d)*t2amps[d,c,i,k]*t2amps[a,b,l,j];
                    Qc += bs.v(k,l,c,d)*t2amps[d,c,j,k]*t2amps[a,b,l,i];



    for k in range(Np):
        for l in range(Np):
            for c in range(Np,Ns):
                for d in range(Np, Ns):
                    Qd -= bs.v(k,l,c,d)*t2amps[a,c,l,k]*t2amps[d,b,i,j];
                    Qd += bs.v(k,l,c,d)*t2amps[b,c,l,k]*t2amps[d,a,i,j];

    return 0.25*Qa + Qb + 0.5*(Qc + Qd)









def CCD_linear_standard(a,b,i,j,t2amp, bs):
    Np = bs.nparticles
    Ns = bs.nstates
    #Standard implementation
    l2a = 0.0
    l2b = 0.0
    l2c = 0.0
    ops = 0
    for c in range(Np, Ns):
        for d in range(Np, Ns):
            ops += 1
            l2a += bs.v(a,b,c,d)*t2amp[c,d,i,j]

    for k in range(Np):
        for l in range(Np):
            ops += 1
            l2b += bs.v(k,l,i,j)*t2amp[a,b,k,l]

    for k in range(Np):
        for c in range(Np, Ns):
            ops += 4
            l2c -= bs.v(a,k,c,j)*t2amp[c,b,i,k];
            l2c += bs.v(a,k,c,i)*t2amp[c,b,j,k];
            l2c += bs.v(b,k,c,j)*t2amp[c,a,i,k];
            l2c -= bs.v(b,k,c,i)*t2amp[c,a,j,k];
            #EQ. (9.120 in S-B)
    #print "Number of operations in standard implementation:", ops
    return .5*l2a + .5*l2b  + l2c

def corr_energy(t2amps,bs):
    e0 = 0
    for n in range(t2amps.N):
        p1,p2,h1,h2,t2val = t2amps.at(n)
        #e0 += .25*bs.v(p1,p2,h1,h2)*t2val
        e0 += .25*bs.v(h1,h2,p1,p2)*t2val
    return e0
    
def advance(Np, Ns, t2amps, t2amp_new, bs):
    print "Initialize new amplitudes"
    #t2amp_new = amplitude()
    for i in range(Np):
        for j in range(i, Np):
            for a in range(Np,Ns):
                for b in range(a, Ns):
                    ampval = CCD_linear_optimized(a,b,i,j,t2amps,bs) + CCD_quadratic_optimized(a,b,i,j,t2amps,bs)
                    
                    #if ampval > 0.0001:
                    t2amp_new.add([a,b,i,j, ampval])
                    #t2amp_new.add([b,a,i,j,-ampval])
                    #t2amp_new.add([a,b,j,i,-ampval])
                    #t2amp_new.add([b,a,j,i, ampval])
            print i,j
    print t2amp_new.N
    return t2amp_new

def stresstest_standard(Nshells, Np):
    bs = electronbasis(Nshells,1.0, Np) 
    bs.liststates()
    print "=============================================="
    print "Performing stresstest for array-implementation"
    print "Number of states   :", bs.nstates
    print "Number of shells   :", Nshells
    print "Number of particles:", Np
    print " "
    print "Results:"
    Ns = bs.nstates
    t2amps = zeros((Ns,Ns,Ns,Ns))
    #Initialization
    t0 = clock()
    for i in range(Np):
        for j in range(i,Np):
            for a in range(Np,Ns):
                for b in range(a,Ns):
                    val = bs.v(i,j,a,b)/(bs.h(i,i) + bs.h(j,j) - bs.h(a,a) - bs.h(b,b))
                    t2amps[a,b,i,j] = val
                    t2amps[b,a,i,j] = -val
                    t2amps[a,b,j,i] = -val
                    t2amps[b,a,j,i] = val
    print "         Initialization time:", clock()-t0, "s."

    Nt = 1
    NT = range(Nt)
    t0 = clock()
    for t in NT:
        val = CCD_quadratic_standard(18,19,4,5, t2amps, bs)
    tq = (clock()-t0)/float(Nt)
    print " Quadratic contribution time:", tq , "s."
    
    t0 = clock()
    for t in NT:
        val = CCD_linear_standard(18,19,4,5, t2amps, bs)
    tl = (clock()-t0)/float(Nt)
    print "    Linear contribution time:", tl, "s."
    print ""
    print "Estimated time per iteration:", (tq+tl)*(Np**2)*(Ns-Np)**2, "s."
    
def stresstest_optimized(Nshells,Np):
    bs = electronbasis(Nshells,1.0, Np)
    Ns = bs.nstates
    print "============================================="
    print "Performing stresstest for list-implementation"
    print "Number of states   :", bs.nstates
    print "Number of shells   :", Nshells
    print "Number of particles:", Np
    print " "
    print "Results:"
    t2amp = amplitude(Ns,Np)
    t0 = clock()
    initialize_t2amps(Np,t2amp,bs)
    init_time = clock()-t0
    print "         Initialization time:", init_time, "s."
    print " Number of non-zero elements:", t2amp.N
    
    t0 = clock()
    t2amp.map_tensor(Ns)
    map_time = clock()-t0
    #print t2amp.getlist_p1(15)
    print "Time spent on sorting tensor:", map_time, "s."
    
    for n in t2amp.getlist_p1(15):
        print t2amp.at(n)
    print " "
    for n in t2amp.getlist_p1p2(15,50):
        print t2amp.at(n)
    
    lin_time = stresstest_linear(100,t2amp,bs)
    quad_time = stresstest_quadratic(100,t2amp,bs)
    print " Quadratic contribution time:", quad_time , "s."
    print "    Linear contribution time:", lin_time, "s."
    print " "
    print "Estimated time per iteration:", map_time + (lin_time+quad_time)*2280 #86100 #(Np**2)*(Ns-Np)**2, "s."
    print " "
    print "Testing actual iteration."
    P = range(Np, Ns)
    H = range(Np)
    t0 = clock()
    for i in H:
        for j in H:
            for a in P:
                for b in P:
                    if bs.kdfullplus(a,b,i,j):
                        val = CCD_quadratic_optimized_mk2(a,b,i,j, t2amp, bs) +  CCD_linear_optimized(a,b,i,j, t2amp, bs)
    print "Actual time spent on 1 iteration (s) :", clock()-t0
    
    
    
    
    
    
    
                    

def CCD(Nt):
    Nshells = 2
    Np = 14
    bs = electronbasis(Nshells,1.0, Np)
    Ns = bs.nstates
    stresstest_interaction(10000,bs)
    
    print "Number of states:", Ns
    #t2amp2 = zeros((Ns,Ns,Ns,Ns), dtype = float)
    t2amp = amplitude(Ns,Np)
    t2amp_new = amplitude()
    initialize_t2amps(Np,t2amp, bs)
    #print "Initial correlation energy:", corr_energy(t2amp,bs), " hartrees."
    t2amp.map_tensor(bs.nstates)
    stresstest_quadratic(1000,t2amp,bs)
    stresstest_linear(1000,t2amp,bs)
    print t2amp.N
    #for t in range(Nt):
    e = corr_energy(t2amp,bs)
    print e
    #print "Correlation energy:", corr_energy(t2amp,bs), " hartrees. (%i)" %1
    
    print "Done mapping."
    advance(Np,Ns,t2amp,t2amp_new, bs)
    print t2amp_new.N
    #print "Final correlation energy:", corr_energy(t2amp,bs), " hartrees."

def stresstest_interaction(Nt,bs):
    Ntt = range(Nt)
    t0 = clock()
    for t in Ntt:
        val = 0*bs.v(4,4,6,12)
    print "Mean interaction time:", (clock()-t0)/float(Nt)
    
def stresstest_quadratic(Nt,t2amps,bs):
    Ntt = range(Nt)
    t0 = clock()
    for t in Ntt:
        val = CCD_quadratic_optimized_mk2(15+t%20,18,3+t%7,1, t2amps, bs)
    tm = (clock()-t0)/float(Nt)
    #print "Mean quadratic time:", (clock()-t0)/float(Nt)
    return tm

def stresstest_linear(Nt,t2amps,bs):
    Ntt = range(Nt)
    t0 = clock()
    for t in Ntt:
        val = CCD_linear_optimized(18,18,3,1, t2amps, bs)
    tm = (clock()-t0)/float(Nt)
    #print "Mean quadratic time:", (clock()-t0)/float(Nt)
    return tm

def setup_matrix():
    bs = electronbasis(4, 1, 14)
    Ns = bs.nstates
    Np = bs.nparticles
    print Ns
    Vab = zeros((Ns**2,Ns**2))
    Ntot = 0
    Ncon = 0
    Nmatrices = 0
    for a in range(Np,Ns):
        for b in range(a,Ns):
            #Mabc = 0
            for i in range(Np):
                for j in range(i,Np):
                    #Vab[a + b*Ns, i + j*Ns] = bs.v(a,b,i,j)
                    #if not bs.v(a,b,i,j) == -bs.v(b,a,i,j) == -bs.v(a,b,j,i) == bs.v(b,a,j,i):
                    #    print a,b,i,j, "unsymmetric"
                    #if not bs.v(a,b,i,j) == bs.v(i,j,a,b):
                    #    print "unsymmetir" 
                    val = bs.v(a,b,i,j)
                    Vab[AB(a,b,i,j,Ns)] = val
                    Vab[AB(b,a,i,j,Ns)] = -val
                    Vab[AB(a,b,j,i,Ns)] = -val
                    Vab[AB(b,a,j,i,Ns)] = val
                    
                    Vab[AB(i,j,a,b,Ns)] = val
                    Vab[AB(i,j,b,a,Ns)] = -val
                    Vab[AB(j,i,a,b,Ns)] = -val
                    Vab[AB(j,i,b,a,Ns)] = val
                    
                    
                    #Ntot += 1
                    #Mab = bs.kdfullplus(a,b,i,j)
                    #Mabc+= Mab
                    #if Mab != 0:
                        
                        #Ncon += 1
                        #print a,b,i,j, bs.kdplus(a,b,i,j)
            #if Mabc != 0:
            #Nmatrices += 1
    #print Ntot, Ncon, Ncon/float(Ntot), Nmatrices
    return Vab

class blockmatrix():
    #taking as input a electron basis and setting up a block-diagonal matrix
    #1. identify all
    def __init__(self, bs):
        self.Ns = bs.nstates
        self.Np = bs.nparticles
        
        

def setup_matrix2():
    bs = electronbasis(3, 1, 14)
    Ns = bs.nstates
    Np = bs.nparticles
    KMs = []
    zmap = zeros((6,6,6), dtype = int)
    #zmap = zeros((6**3), dtype = int)
    zmap_pq = zeros((Ns,Ns), dtype = int)
    for p in range(Ns):
        for q in range(p,Ns):
            kms = array(bs.states[p,1:5]) + array(bs.states[q,1:5])
            zmap[2+kms[0],2+kms[1], 2+kms[2]] += 2
            zmap_pq[p,q] = 0
            #if kms not in KMs:
            #KMs.append(kms)
    #print len(KMs),KMs
    print zmap
    return KMs

def AB(p,q,r,s, Ns):
    return p + q*Ns, r + s*Ns

def A(p,q, Ns):
    return p + q*Ns

def B(r,s, Ns):
    #i know, it only for human convenience
    return r + Ns*s

def count_ops():
    ops = 0
    bs = electronbasis(3, 1, 14)
    Ns = bs.nstates
    Np = bs.nparticles
    for a in range(Np,Ns):
        for b in range(Np,Ns):
            for i in range(Np):
                for j in range(Np):
                    if bs.kdfullplus(a,b,i,j):
                        ops += 1
    print "Number of operations needed:", ops




#CCD(4)
#stresstest_standard(3,14)
#count_ops()
#stresstest_optimized(3,14)
V = setup_matrix2()
#print V.max(), V.min()
#print V[A(6,4,54), :]
#from matplotlib.pyplot import *
#imshow(V)
#show()
