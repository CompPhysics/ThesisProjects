from IPython.display import display, Math, Latex 
#%matplotlib inline  

from numpy import *
from itertools import *
from matplotlib.pyplot import *
from copy import copy

class Combiner():
    #Normal ordered operator for cluster algebra (diagrammatic)
    def __init__(self, H, T):       
        self.q_a = H
        self.I = []
        HN = len(self.q_a)
        self.combine(H,T)
        
    def combine(self, H, ops, excitation = None):
        #Assosciate a list of T-operators (ops) to the current operator instance    
        #Find all possible ways to combine the operators using self.distinct_combinations() 
        T = []
        for i in ops:
            T.append(i)

        self.T_operator = T
        
        #Finding acceptable combinations of internal contractions between the operators
        self.I = self.distinct_combinations(self.q_a, T)
               
    def distinct_combinations(self,H,T):
        #Returns all possible combinations of H and T
        #I - all q-particle annihilation operators
        #T list of list with T operators. ex. [[-1,1],[-1,1]]  = T_1 T_1
       lenH  = len(H)
       lenT  = len(T)
       lenTi = []
       for i in range(lenT):
           lenTi.append(len(T[i]))
       #H+=[0 for i in range(lenT-1)] #adding zeros to denote separations in the cluster-operators

       H2 = []
       for i in H:
           H2.append(i)
       H2 += [0 for i in range(lenT-1)] #adding zeros to denote separations in the cluster-operators

       #Creating countlist for T-operator to keep track of q-operators in each cluster operator
       T_budget = self.itemcount(T)
       
       
       #Create all permutations
       H_permuted = permutations(H2)
       
       #Sort out indistinct diagrams and cancelling terms
       excluded = []
       excluded_budgets = []
       accepted = []
       for i in H_permuted:
           if self.acceptable(i, T, excluded, excluded_budgets):
               #print "Accepted"
               accepted.append(self.splitlist(i,0))       
       self.combined = accepted
               
           
       return accepted
       #def evaluate_q_c_total(self):
                
    def assess_excitation(self):
        #Assess current excitation level (also if combined operator)
        self.E = self.q_c.count(-1) + self.q_c.count(1) - (self.q_a.count(1) + self.q_a.count(-1))
        for i in range(len(self.T_operator)):
            self.E += self.T_operator[i].count(1) + self.T_operator[i].count(-1) 
        #fill in contracted elements
        self.E/= 2.0
    def scan_extract(L, e):
        #Exctract element e from list L
        ret = None
        for i in range(len(L)):
            if L[i] == e:
                L = delete(L, i)
                ret = e
    def nloops(self,x,y):
        #returns number of loops in budget
        return (x+y - abs(x-y))//2
    
  
    def nozeroedges(self,i):
        #Assert that there are no borders in the endpoints of the connection pattern i
        ret = True
        try:
            if i[0] == 0 or i[-1] == 0:
                ret = False
        except:
            ret = False
        return ret
    
    def nozerocontact(self,i):
        #Assert that there are no neighbouring borders in the connection pattern i
        ret = True
        for e in range(len(i) - 1):
            if i[e+1] == 0 and i[e] == 0:
                ret = False
        return ret
    
    def splitlist(self,L, d):
        #Returns a split list HL from a list L into constituents, d denotes barrier
        HL = [[]]
        n = 0
        for i in range(len(L)):
            if L[i] != d:
                HL[n].append(L[i])        
            if L[i] == d:
                HL.append([])
                n += 1        
        return HL
    
    def itemcount(self,T):
        #Count number of particle- and hole lines in each constituent part of T
        #Returns a list object of the form [[#particles, #holes], ...]
        itemnumber = []
        for i in range(len(T)):
            itemnumber.append([])
            itemnumber[i].append(T[i].count(1))  #number of q-particles
            itemnumber[i].append(T[i].count(-1)) #number of q-holes
        return itemnumber
    
    def contractable(self,L,T):
        #Asserts that the number of contractions in each p- and hline of L does not superseed # of p-h in T 
        #input two itemcount items lists, returns bool
        ret = True
        for i in range(len(T)):
            for e in range(len(T[i])):
                if T[i][e]<L[i][e]:
                    ret = False
        return ret
        
    def find_identical(self,T):
        #Find identical operators in the list of operators T
        #returns a list of pairs of indices that denotes permutations that does not alter the T operator
        identicals = []
        for i in range(len(T)):
            for e in range(i,len(T)):
                if T[i] == T[e] and i!=e:
                    identicals.append([i,e])
        return identicals
                    
    def permute_elements(self,e1,e2,L):
        #Returns a list where elements at indices e1,e2 in L is permuted
        L_ret = []        
        for i in L:
            L_ret.append(i)
        L_ret[e1] = L[e2]
        L_ret[e2] = L[e1]
        return L_ret
    
    def acceptable(self,i, T, excluded, excluded_budgets):
        #Test if a potential connection pattern is distinct
        #Returns bool
        ret = False
        identicals = self.find_identical(T)
        T_budget   = self.itemcount(T)
        
        if self.nozeroedges(i) and self.nozerocontact(i):
            I = self.splitlist(i, 0)
            I_budget = self.itemcount(I)
            if I_budget not in excluded_budgets:
                excluded_budgets.append(I_budget)
                for e in identicals:
                    excluded_budgets.append(self.permute_elements(e[0], e[1], I_budget))
                if self.contractable(I_budget,T_budget):
                    if I not in excluded:
                        excluded.append(I)
                        for e in identicals:
                            excluded.append(self.permute_elements(e[0], e[1], I))
                        ret = True
        return ret

class qnode():
    #Object for keeping track of vertices
    def __init__(self, qa):
        #self.pos = pos
        self.connected = 0
        self.qa = qa #creation or annihilation +1/-1
        self.group = []
    def connect(i):
        self.connected = i

class qop():
    #a collection of vertices -- in principle an operator
    def __init__(self,qnodes):
        self.qnodes = qnodes
        self.pos = [0,0]
        self.setpos(self.pos)
    def setpos(self, pos):
        self.pos = pos
        for i in range(len(self.qnodes)):
            self.qnodes[i].qpos = self.pos[0] + i
            

class qvert():
    #A vertex in the operator with input, output
    def __init__(self, qin, qout):
        self.qin = qin   #qnode in
        self.qout = qout #qnode out
        self.qpos = 0.0


def qlist(o_cop):
    O = o_cop
    excluded = []
    olist = []
    for i in range(len(O)):
        if O[i] == +1:
            #find -1
            for e in range(len(O)):
                if e!=i:
                    if O[e] == -1:
                        if e not in excluded:
                            if i not in excluded:
                                excluded.append(i)
                                excluded.append(e)
                                olist.append(qvert(qnode(1), qnode(-1)))
    return qop(olist)

    
def connect(H,T,I):
    #connect H and T using I as pattern
    #Return a mapping of the connections (a list of H-connections to T)
    connected_H = []
    connected_T = []    
    map_T = []
    map_H = []
    for i in range(len(H)):
        map_H.append(None)
    for e in range(len(T)):
        map_T.append([])
        for u in range(len(T[e])):
            map_T[e].append(None)
    
    for i in range(len(I)):
        for e in range(len(I[i])):
            t = I[i][e] #connect this quantity to H and T
            scan = True
            for e1 in range(len(T[i])):
                if T[i][e1] == t and map_T[i][e1] == None:
                    for i1 in range(len(H)):
                        if H[i1] == t and map_H[i1] == None:
                            if scan:
                                #Perform connection
                                map_H[i1] = [i,e1]
                                map_T[i][e1] = i1
                                scan = False
    return map_T


class Operator():
    def __init__(self, q_a, q_c, pos = [0,0]):
        self.q_a = q_a
        self.q_c = q_c
        self.internals = len(q_c)
        self.children = []
        self.E = self.excitation()
        self.pos = pos
        self.pair_up_lines()
    def excitation(self):
        #Evaluate operator excitation level
        Ex = (self.q_c.count(-1) + self.q_c.count(1) - (self.q_a.count(1) + self.q_a.count(-1)))/2.0
        for i in range(len(self.children)):
            Ex += self.children[i].excitation()
        return Ex
    def adopt(self, O, I):
        self.children = O
        self.relation = I
        for i in self.children:
            i.set_pos([self.pos[0], self.pos[1]-1])
        #for i in self.relation:
        #    if 
        #append new child operator to this operator, include more outgoing lines, use connection pattern I
        #C = Combiner(self.q_c, O.q_a) #This must be done externally
        return 0
    def get_schematic(self):
        schema = []
        for i in self.children:
            schema.append(i.get_schematic())
        #return a list of connected nodes and type of connections
        
    def set_pos(self, pos):
        self.pos = pos
        for i in self.children:
            i.pos =[pos[0], pos[1]-1]
    def get_c_pos(self,ci):
        ret = None
        if ci<self.internals:
            for i in range(len(self.pos_group)):
                if ci in self.pos_group[i][0]:
                    ret = [self.pos[0]+i, self.pos[1]]
                    break
        else:
            ret = 0
            
        return ret
        
    def pair_up_lines(self):
        #Pair up ***internal*** vertices
        p0 = 0 #Vertical position for iterations
        #Set up a mapping
        pos_groups = [] #Grouping 2 and 2 nodes together
        paired = [[],[]]
        #for i in range((len(self.q_a)+len(self.q_c))/2):
        #    pos_groups.append([])
        
        #Use same principle as below (schematic) to pair them up and map them out
        self.m_qa = []
        self.m_qc = []

        self.inout = [[],[]] #keep track of vertices in operator

        for i in range(len(self.q_c)):
            self.m_qc.append(None)
        for i in range(len(self.q_a)):
            self.m_qa.append(None)

        
        #identify pairs passing through the operator
        for i in range(self.internals):
            for e in range(len(self.q_a)):
                if self.q_a[e] == self.q_c[i] and self.m_qc[i] == None and self.m_qa[e] == None:
                    #paired[0].append(i)
                    #paired[1].append(e)
                    #pos_groups.append([[i],[e]])
                    self.m_qc[i] = e+1
                    self.m_qa[e] = i+1
                    break
                    
        #identify pairs beneath the interaction
        for i in range(len(self.q_a)):
            for e in range(len(self.q_a)):
                if i!=e and self.q_a[e] == -self.q_a[i] and self.m_qa[i] == None and self.m_qa[e] == None:
                    self.m_qa[i] = -e-1 #negative index implies index inside self
                    self.m_qa[e] = -i-1
                    break


        #identify pairs above the interaction
        for i in range(self.internals):
            for e in range(self.internals):
                if i!=e and self.q_c[e] == -self.q_c[i] and self.m_qc[i] == None and self.m_qc[e] == None:
                    self.m_qc[i] = -e-1
                    self.m_qc[e] = -i-1
                    break                                   
        
        

    def get_paired_c(self, p):
        #return the corresponding index of the line below the vertex
        ret = None
        for i in range(len(self.pos_group)):
            if p in self.pos_group[i][0]:
                ret = self.pos_group[i][1][0]
        return ret

class contracted():
    def __init__(self, H, T, pattern):
        self.H = H
        self.T = T
        self.pattern = pattern
        
        

def contract(H,T):
    T_q_c = []
    for i in T:
        i.set_pos([H.pos[0], H.pos[1]-1])
        T_q_c.append(i.q_c)
    C = Combiner(H.q_a,T_q_c) #Combines annihilators from H (Hqa) and creators from T (Tqc)
    diagrams = []
    for i in C.I:
        diagrams.append([copy(H), copy(T), connect(H.q_a, T_q_c, i), i])
    #Pending: optional combinations of operators (T^-1 H T) for property calculations
    return diagrams


class schematic():
    def __init__(self, diagram):
        print "Connection " , diagram[2]
        print "Number ", diagram[3]
        self.H = diagram[0]            
        self.T = diagram[1] 
        self.schema = diagram[2]
        self.i = diagram[3]
    
    def scan_extract(L, e):
        #Exctract element e from list L
        ret = None
        for i in range(len(L)):
            if L[i] == e:
                L = delete(L, i)
                ret = e
    def nloops(self,x,y):
        #returns number of loops in budget
        return (x+y - abs(x-y))//2
    def itemcount(self,T):
        #Count number of particle- and hole lines in each constituent part of T
        #Returns a list object of the form [[#particles, #holes], ...]
        itemnumber = []
        for i in range(len(T)):
            itemnumber.append([])
            itemnumber[i].append(T[i].count(1))  #number of q-particles
            itemnumber[i].append(T[i].count(-1)) #number of q-holes
        return itemnumber
    
    def contractable(self,L,T):
        #Asserts that the number of contractions in each p- and hline of L does not superseed # of p-h in T 
        #input two itemcount items lists, returns bool
        ret = True
        for i in range(len(T)):
            for e in range(len(T[i])):
                if T[i][e]<L[i][e]:
                    #print T[i][e],L[i][e]
                    ret = False
        return ret
        
    def find_identical(self,T):
        #Find identical operators in the list of operators T
        #returns a list of pairs of indices that denotes permutations that does not alter the T operator
        identicals = []
        for i in range(len(T)):
            for e in range(i,len(T)):
                if T[i] == T[e] and i!=e:
                    identicals.append([i,e])
        return identicals
                    
    def permute_elements(self,e1,e2,L):
        #Returns a list where elements at indices e1,e2 in L is permuted
        L_ret = []        
        for i in L:
            L_ret.append(i)
        L_ret[e1] = L[e2]
        L_ret[e2] = L[e1]
        return L_ret
    
    def acceptable(self,i, T, excluded, excluded_budgets):
        #Test if a potential connection pattern is distinct
        #Returns bool
        ret = False
        identicals = self.find_identical(T)
        T_budget   = self.itemcount(T)
        
        if self.nozeroedges(i) and self.nozerocontact(i):
            I = self.splitlist(i, 0)
            I_budget = self.itemcount(I)
            if I_budget not in excluded_budgets:
                excluded_budgets.append(I_budget)
                for e in identicals:
                    excluded_budgets.append(self.permute_elements(e[0], e[1], I_budget))
                if self.contractable(I_budget,T_budget):
                    if I not in excluded:
                        excluded.append(I)
                        for e in identicals:
                            excluded.append(self.permute_elements(e[0], e[1], I))
                        ret = True
        return ret
    def setup(self):
        plabels_t = ['r','s','x','z','b','a'] #Target labels
        hlabels_t = ['t','u','q', 'w', 'j','i']

        plabels = ['h','g','f','e','d','c']
        hlabels = ['p','o','n','m','l','k']

                        
        #plabels_ex =  ['h','g','f','e','d','c']
        #hlabels_ex =  ['p','o','n','m','l','k']

        self.T_labels = []
        for i in range(len(self.T)):
            self.T_labels.append([])
            for e in range(len(self.T[i].q_c)):
                self.T_labels[i].append(None)

        self.H_labels_c = []
        self.H_labels_a = []        
        for i in range(len(self.H.q_a)):
            self.H_labels_a.append(None)
        for i in range(len(self.H.q_c)):
            self.H_labels_c.append(None)

        #Label all external lines exiting the cluster operators (targets)
        for i in range(len(self.schema)):
            for e in range(len(self.schema[i])):
                if self.schema[i][e] == None:
                    if self.T[i].q_c[e] == 1:
                        self.T_labels[i][e] = plabels_t.pop()
                    if self.T[i].q_c[e] == -1:
                        self.T_labels[i][e] = hlabels_t.pop()
                
        #labels all external lines exiting the interaction
        for i in range(len(self.H.q_c)):
            
            if self.H.q_c[i] == 1:
                self.H_labels_c[i] = plabels_t.pop()
            if self.H.q_c[i] == -1:
                self.H_labels_c[i] = hlabels_t.pop()
                
        
        #label all internal lines
        for i in range(len(self.schema)):
            n = len(self.schema[i])
            for e in range(n):
                if self.schema[i][e] != None:
                    #Current index is connected
                    if self.T[i].q_c[e]==1:
                        self.T_labels[i][e] = plabels.pop()
                    if self.T[i].q_c[e]==-1:
                        self.T_labels[i][e] = hlabels.pop()
                    self.H_labels_a[self.schema[i][e]] = self.T_labels[i][e]
        
        #Identify the cluster-operators
        self.t_operators()
        
        #Identify the interaction operator
        self.v_operator()
        
        
        #Identify summation indices
        self.summation = ""
        for i in range(len(self.schema)):
            for e in range(len(self.schema[i])):
                if self.schema[i][e] != None:
                    self.summation += self.T_labels[i][e]


        #Connect the vertices in H_labels
        #Identify vertices in all operators

        
        
        
        self.count_loops()
        self.count_holes()
        self.count_equivalent_l()
        self.count_equivalent_t()
        self.calc_prefactor()
        self.find_permutations()


    def cancel_permutations(self):
        return 0
        
        
    def find_permutations(self):
        #Identify name, origin and orientation of outgoing lines
        orig = []
        for i in range(len(self.schema)):
            for e in range(len(self.schema[i])):
                if self.schema[i][e] == None:
                    orig.append([self.T_labels[i][e], self.T[i].q_c[e], "t"+str(i)])
        for i in range(len(self.H.q_c)):
            orig.append([self.H_labels_c[i], self.H.q_c[i], "v"])

        permutes = []
        #compare lines
        for i in range(len(orig)):
            for e in range(i,len(orig)):
                if i!=e and self.ineq_external(orig[i],orig[e]):
                    permutes.append(str(orig[i][0])+str(orig[e][0]))
        self.permutations = permutes
                    
                    
        return 0
    def ineq_external(self, line1, line2):
        #compare two external lines, decide wether they are inequivalent (permutations)
        ret = False
        if line1[1] == line2[1]: #equal orientation
            if line1[2] != line2[2]: #inequal origin #AND origins differ 
                ret = True
        return ret
        
    def calc_prefactor(self):
        self.prefactor = (0.5**self.n_equivalent_Ts)*(0.5**self.equiv_lines)*((-1)**(self.holes-self.loops))
    def count_equivalent_t(self):
        #count equivalent T-vertices
        t_config = []
        for i in range(len(self.T)):
            eq_h = 0
            eq_p = 0            
            for e in range(len(self.T[i].q_c)):
                if self.schema[i][e] != None:
                    if self.T[i].q_c[e] == -1:
                        eq_h += 1
                    if self.T[i].q_c[e] == 1:
                        eq_p += 1

            t_config.append([eq_p, eq_h])
        ignored = []
        t_eqs = 0
        for i in range(len(t_config)):
            for e in range(len(t_config)):
                if i!= e and i not in ignored:

                    if t_config[i] == t_config[e] and len(self.T[i].q_c) == len(self.T[e].q_c):
                        ignored.append(e)
                        t_eqs += 1
        self.n_equivalent_Ts = t_eqs 
        
    def count_equivalent_l(self):
        #Count equivalent lines
        eqs = 0
        for i in range(len(self.T)):
            eq_h = 0
            eq_p = 0            
            for e in range(len(self.T[i].q_c)):
                if self.schema[i][e] != None:
                    if self.T[i].q_c[e] == -1:
                        eq_h += 1
                    if self.T[i].q_c[e] == 1:
                        eq_p += 1
            eqs += eq_h/2 + eq_p/2
        self.equiv_lines = eqs
                
                
        
    def count_loops(self):
        lp = 0
        for i in range(len(self.T)):
            lp += (len(self.T[i].q_c) - sum(self.T[i].q_c))/2
        self.loops = lp
    def count_holes(self):
        nh = 0
        for i in range(len(self.T)):
            for e in range(len(self.T[i].q_c)):
                if self.T[i].q_c[e] == -1:
                    nh += 1
        for i in range(len(self.H.q_c)):
            if self.H.q_c[i] == -1:
                nh += 1
        self.holes = nh
                
    def t_operators(self):
        tt = []
        #print self.T_labels
        for i in range(len(self.T_labels)):
            #tt.append([])
            n = len(self.T_labels[i])/2
            a = ""
            j = ""
            for e in range(n):
                a += self.T_labels[i][e] 
                #a += self.T_labels[i][n-(e+1)] 
                #j += self.T_labels[i][n+e]
                j += self.T_labels[i][n+e]
            tt.append([n,a,j])
        self.t_ops = tt
    def v_operator(self):
        hh = []
        out = ""
        ins = ""
        for i in range(len(self.H.m_qc)):
            l = self.H.m_qc[i]
            if  l>0:
                l-=1

                #found a line passing through the interaction
                if self.H.q_a[l] == 1:
                    #Found a particle passing through
                    ins += self.H_labels_c[i]
                    out += self.H_labels_a[l]
                if self.H.q_a[l] == -1:
                    #Found a hole passing through
                    ins += self.H_labels_c[i]
                    out += self.H_labels_a[l]
                    
                    
        ######
        ######
        # This needs to happen in the correct ordeR!
        ######
        ######
                    
                    
                    
                    
        a_ignores = []            
        for i in range(len(self.H.m_qa)):
            l = self.H.m_qa[i]
            if  l<0 and l not in a_ignores:
                #found a line pair annihilating at the interaction
                l = abs(l) - 1
                out += self.H_labels_a[i]
                ins += self.H_labels_a[l]
                a_ignores.append(-(i+1))

        c_ignores = []
        for i in range(len(self.H.m_qc)):
            l = self.H.m_qc[i]
            if  l<0 and l not in c_ignores:
                #found a line pair created at the interaction
                l = abs(l) - 1
                out += self.H_labels_c[i]
                ins += self.H_labels_c[l]  
                c_ignores.append(-(i+1))                  
                                    
  
        self.h_op = [ins, out]
    def report(self):
        #print " Tlabels :", self.T_labels
        #print " Hlabelsa:", self.H_labels_a
        #print " Hlabelsc:", self.H_labels_c 
        #print "Equivavlent lines:", self.equiv_lines
        #print "Loops:", self.loops
        #print "Holes:", self.holes
        #print "Equivavlent Ts:",self.n_equivalent_Ts 
        #print "Prefactor   :", self.prefactor        
        #print "Permutations:", self.permutations
        #print "Summations  :", self.summation
        #print "Interaction :", self.h_op
        #print "Cluster     :", self.t_ops
        s = "prefactor: %f" % self.prefactor + "\n"
        s += "permutations: %s" % self.permutations+ "\n"
        s += "sum: %s" % self.summation+ "\n"
        s += "h_ops: %s" % self.h_op+ "\n"
        s += "t_ops: %s" % self.t_ops+ "\n"
        return s

                    

class cpp_func():
    def __init__(self, schema, target):
        self.schema = schema
        self.build(target)
    def permute(self, c, s1, s2):     
        c2 = []
        for i in c:
            c2.append(i)
        #print c2
        s1_i = None
        s2_i = None

        for i in range(len(c)):
            if c[i] == s1:
                s1_i = i
                c2[s1_i] = s2
            if c[i] == s2:
                s2_i = i
                c2[s2_i] = s1
            
        c3 = ""
        for i in c2:
            c3 += i

        return c + "-(" + c3 +")"
        
    def build(self, target):
        prefactor = self.schema.prefactor
        permutations = self.schema.permutations
        sums = self.schema.summation
        interaction = self.schema.h_op
        cluster = self.schema.t_ops
        
        #Begin building the function
        s = "double %s = 0.0\n" %target
        c = ""
        #Core function
        c1 = interaction[0][0]
        c2 = interaction[1][0]
        for i in range(len(interaction[0])-1):
            c1 += ","+interaction[0][i+1]
            c2 += ","+interaction[1][i+1]     
        c += "%s(%s)(%s)" % ("v"+str(len(interaction[0])),c1, c2) #Two particle interaction
        for i in range(len(cluster)):
            c1 = cluster[i][1][0]
            c2 = cluster[i][2][0]
            for e in range(len(cluster[i][1])-1):
                c1 += ","+cluster[i][1][e+1]
                c2 += ","+cluster[i][2][e+1]
            c += "*%s(%s)(%s)" % ("t"+str(cluster[i][0]),c1, c2) #Two particle interaction
        
        #Permute all possible
        for i in range(len(permutations)):
                c = self.permute(c, permutations[i][0], permutations[i][1])
                
            
        
        c = "%s += " %target + c
        c+=";\n"
        sms = []
        for i in range(len(sums)):
            indent = 4 * (len(sums)-i-1) * " "
            s,e = self.for_loop(sums[i])
            sms.insert(0,indent + s)
            sms.append(indent + e)
        indent = 4*(len(sums)) * " " 
        sms.insert(len(sms)/2, indent+c)
        
        #merge function
        C=""
        for i in sms:
            C+= i
            
        C = "double %s = 0.0\n"%target+ C 

        C +="double %s *= %f\n" % (target, prefactor)
        #Set up summations
        #print "\n", C
        self.cppfunction = "\n" + C
        self.target = target
        #print sms
            
    
    def for_loop(self, p):
        if p in "abcdefgh":
            #return 
            lps = "for(int %s = nElectrons; %s < nStates; %s ++){\n" % (p,p,p)
            lpe= "}\n"
        if p in "ijklmnop":
            #return 
            lps = "for(int %s = 0; %s < nElectrons; %s ++){\n" % (p,p,p)
            lpe = "}\n"
        return lps, lpe  

class symbolic_diag():
    def __init__(self, schema, target):
        self.schema = schema
        self.target = target
        self.build(target)
    def build(self, target):
        prefactor = self.schema.prefactor
        permutations = self.schema.permutations
        sums = self.schema.summation
        interaction = self.schema.h_op
        cluster = self.schema.t_ops
        
        pre_fact, pre_denom = self.factorize(prefactor)
        s = ""
        for i in permutations:
            s += "P(%s)" % i
        s += "\\frac{%i}{%i}" % (pre_fact, pre_denom)
        s += " \\sum_{%s}" % sums
        s += " \\langle %s || %s \\rangle" % (interaction[0],interaction[1])
        for i in cluster:
            s += " t_{%s}^{%s}" % (i[2], i[1])
        self.stringform = s
    
    def factorize(self, prefactor):
        ret = prefactor
        fract = 1
        while not float(ret) == int(ret):
            ret*=2
            fract*=2
            
        return ret, fract
     
class feynman_diag():
    #a diagrammatic representation of the operators
    def __init__(self, schema, target, pos, vis = True):
        self.schema = schema
        self.target = target
        self.pos = pos
        self.vis = vis
        
        
        #Parameters for plotting
        self.shift = .2 #horizontal shift of interaction vertices
        self.scale = 1.0 #size of diagrams
        self.split = .25 #angle of paired leaving vertices
        
        #building diagrams
        self.Build(target)
    def Build(self,target):
        scale = self.scale
        prefactor = self.schema.prefactor
        permutations = self.schema.permutations
        sums = self.schema.summation
        interaction = self.schema.h_op
        cluster = self.schema.t_ops
        schema = self.schema.schema
        px = self.pos[0]
        py = self.pos[1]
        
        #Create all vertices in the cluster operators
        #Find total number of vertices:
        nt = 0
        for i in range(len(cluster)):
            nt += cluster[i][0]
        #print "NT:", nt
        self.t_nodes = []
        self.exit_nodes = []
        t0 = -scale*nt/2 
        self.NT_tot = nt
       #print cluster
        for i in range(len(cluster)):
            self.t_nodes.append([])
            self.exit_nodes.append([])
            for e in range(cluster[i][0]):
                self.t_nodes[i].append(hnode([px+t0, py],  [1,-1],[]))
                #print px+t0, py
                self.exit_nodes[i].append([])
                self.exit_nodes[i][-1].append(hnode([px+t0-self.split*scale, py+2*scale],[],[]))
                self.exit_nodes[i][-1].append(hnode([px+t0+self.split*scale, py+2*scale],[],[]))
                t0 += scale

        
        #Set up H-vertices
        #counting hole lines passing through the interaction
        n_h_h = []
        for i in range(10):
            if i in self.schema.H.m_qc:
                n_h_h.append(self.schema.H.q_a[i-1])
        n_passing_holes = n_h_h.count(-1)
        n_passing_parts = n_h_h.count(1)
        #print n_h_h

        #count particle lines passing 
        n_h_h = 0
        for i in range(10):
            n_h_h += self.schema.H.m_qc.count(i)
        
                        
        #counting the line pairs annihilating at the interaction
        n_h_a = 0
        for i in range(10):
            n_h_a += self.schema.H.m_qa.count(-i)

        #counting the line pairs created at the interaction
        n_h_c = 0
        for i in range(10):
            n_h_c += self.schema.H.m_qc.count(-i)
        
        n_h_a /= 2
        n_h_c /= 2
        
        #print "Passing holes:", n_passing_holes
        #print "Passing particles:" , n_passing_parts
        #print "Annihil:", n_h_a
        #print "Created:", n_h_c
        
        self.h_nodes = []
        h0 = - .5*(n_h_a + n_passing_holes +n_passing_parts + n_h_c)*scale +  self.shift
        for i in range(n_h_a):
            #Add vertices that annihilate a pair of lines
            self.h_nodes.append(hnode([px+h0,py+1*scale], [],[1,-1]))
            h0 += scale 
        for i in range(n_h_c):
            #Add vertices that create a pair of lines
            self.h_nodes.append(hnode([px+h0,py+1*scale], [1,-1],[]))
            h0 += scale
        for i in range(n_passing_parts):
            self.h_nodes.append(hnode([px+h0,py+1*scale], [1],[1]))
            h0 += scale            
        for i in range(n_passing_holes):
            self.h_nodes.append(hnode([px+h0,py+1*scale], [-1],[-1]))
            h0 += scale
        #Done setting up the H-vertices
        
        
        #Grouping up lines
        TL = self.schema.T_labels
        HLa = self.schema.H_labels_a
        HLc = self.schema.H_labels_c   

        lines_t = [] #list of all connections in the diagram, both external and internal
        lines_t_internal = [] #list of all connections in the diagram, both external and internal
        #for i in range(len(self.schema)):
        lines_h_out = []
        
        
        HLc = self.schema.H_labels_c 
        #import mapping 
        mqc = self.schema.H.m_qc
        mqa = self.schema.H.m_qa
        
        for i in range(len(schema)):            
            lines_t.append([])
            lines_t_internal.append([])            
            for e in range(len(schema[i])):
                if schema[i][e] == None:
                    lines_t[i].append([TL[i][e], [i,e], None])
                if schema[i][e] != None:
                    lines_t[i].append([TL[i][e], [i,e], i])
        for i in range(len(HLc)):
            lines_h_out.append([HLc[i], i, None])

        self.lines_h_out = lines_h_out
        self.lines_t = lines_t
        self.lines_t_int = lines_t_internal
        
        
        
        
        
        self.Plot(target)
    def Plot(self, target):
        if self.vis:
            figure(figsize = (2, self.NT_tot/1.5), dpi = 80, edgecolor = 'k',facecolor='white')
            axis('off')
            axes().set_aspect('equal', 'datalim')
        hold('on')
        t_nodes = self.t_nodes
        h_nodes = self.h_nodes
        exit_nodes = self.exit_nodes

        
        #draw all nodes and horizontal lines
        pos = self.pos
        scale = self.scale
        cluster = self.schema.t_ops
        lines_t = self.lines_t
        msize =5
        c = 'black'
        for i in range(len(self.t_nodes)):
            if len(self.t_nodes[i])==1:
                nc = self.t_nodes[i][0] #current node
                #Draw 
                plot( [nc.x - 0.25*scale,nc.x + 0.25*scale],[nc.y, nc.y], color = c)
                plot(nc.x, nc.y,"o", color = c, markersize = msize)
            if len(self.t_nodes[i])>1:
                nc_next = self.t_nodes[i][0]                
                plot(nc_next.x, nc_next.y,"o", color = c, markersize = msize)
                for e in range(len(t_nodes[i])-1):
                    nc = self.t_nodes[i][e]
                    nc_next = self.t_nodes[i][e+1]
                    
                    plot( [nc.x,nc_next.x],[nc.y, nc_next.y], color = c)
                    plot(nc_next.x, nc_next.y,"o", color = c, markersize = msize)
                    
        
        #Draw interaction vertices
        if len(self.h_nodes)==1:
            nc = self.h_nodes[0]
            plot(nc.x,nc.y,"o", color = c, markersize = msize)
            plot( [nc.x,nc.x + 0.5*scale],[nc.y,nc.y], color = c, ls = "dotted")
        if len(self.h_nodes)>1:
            #connect horizontal vertices
            nc = self.h_nodes[0]
            plot(nc.x,nc.y,"o", color = c, markersize = msize)
            for e in range(len(self.h_nodes)-1):
                nc = self.h_nodes[e]
                nc_next = self.h_nodes[e+1]
                plot( [nc.x,nc_next.x],[nc.y, nc_next.y], color = c,ls = "dotted")
                plot(nc_next.x, nc_next.y,"o", color = c, markersize = msize)
                
        #Draw all connected and unconnected lines
        for i in range(len(cluster)): 
            for u in range(cluster[i][0]):
                p1_ind = self.get_l_index(lines_t[i],cluster[i][1][u])
                h1_ind = self.get_l_index(lines_t[i],cluster[i][2][u])   
                
                #print lines_t[i][p1_ind] , lines_t[i][h1_ind]                 
                if lines_t[i][p1_ind][2] == lines_t[i][h1_ind][2] == None:
                    #print "Two external unconnected lines!"
                    ncon(t_nodes[i][u],exit_nodes[i][u][0],order = 0, p_h = -1)
                    ncon(t_nodes[i][u],exit_nodes[i][u][1],order = 0, p_h = 1)
                    #print lines_t[i][p1_ind], lines_t[i][h1_ind]           
                if lines_t[i][p1_ind][2] == lines_t[i][h1_ind][2] != None:
                #if lines_t[i][p1_ind][2] == lines_t[i][h1_ind][2] == 0:
                    #print "Internal loop!"
                    #Scan for available loop
                    n = None
                    for t in range(len(h_nodes)):
                        #print "t:",t
                        if h_nodes[t].probe_loop_a():
                            #print "Found loop!"
                            #Found loop
                            h_nodes[t].connect_loop_a()
                            n = t
                            break
                    try:
                        hn = h_nodes[n]
                        
                            
                        ncon(t_nodes[i][u],hn,order = -1, p_h = 1)
                        ncon(hn,t_nodes[i][u],order = -1, p_h = 1)
                        
                    #Remove if causing trouble
                    except:
                        #connect to separate vertices
                        
                        #connect particle line
                        for t in range(len(h_nodes)):

                            if h_nodes[t].probe_line_a(-1):
                                #print "Found particle!"

                                h_nodes[t].connect_loop_a()
                                n = t
                                break
                        hn = h_nodes[n]
                        ncon(t_nodes[i][u],hn,order = 0, p_h = -1)
                        
                        #connect particle line
                        for t in range(len(h_nodes)):

                            if h_nodes[t].probe_line_a(1):
                                #print "Found particle!"

                                h_nodes[t].connect_loop_a()
                                n = t
                                break
                        hn = h_nodes[n]
                        ncon(t_nodes[i][u],hn,order = 0, p_h = 1)
                        
                            
                        #print "Could not connect loop."
        for i in range(len(cluster)): 
            for u in range(cluster[i][0]):
                p1_ind = self.get_l_index(lines_t[i],cluster[i][1][u])
                h1_ind = self.get_l_index(lines_t[i],cluster[i][2][u])   
                
                #print lines_t[i][p1_ind] , lines_t[i][h1_ind]          
                if lines_t[i][p1_ind][2] != lines_t[i][h1_ind][2]:
                    #print "One line in , one line out"
                    #print lines_t[i][p1_ind], lines_t[i][h1_ind]
                    
                    #draw line lines_t[i][p1_ind]
                    #begin = [self.t_pos[i][u],self.h_pos[hn/2]]
                    #print "HPOS:          ", self.h_pos, hn/2, hn
                    #print "TPOS:          ", self.t_pos,i,u
                    #begin = [self.t_pos[i][u],self.h_pos[hn/2]] #X-vector
                    #plot(begin, "o")
                    
                    if lines_t[i][p1_ind][2] == None:
                        #print "    Draw external particle"
                        ncon(t_nodes[i][u],exit_nodes[i][u][0],order = 0, p_h = 1)
                        #end = [self.pos[1],self.pos[1]+2]
                        #begin = [self.t_pos[i][u],self.t_pos[i][u]] #X-vector
                        #self.draw_lines(1, begin,end, 1)
                        #hn += 1                        
                    if lines_t[i][p1_ind][2] != None:
                        #print "    Draw internal particle"
                        #Scan for available particle annihilator
                        n = None
                        for t in range(len(h_nodes)):
                            if h_nodes[t].probe_line_a(1):
                                #Found loop
                                h_nodes[t].connect_line_a(1)
                                n =t
                                break
                        hn = h_nodes[n]
                        
                        ncon(t_nodes[i][u],hn,order = 0, p_h = 1)
                        
                        
                        #begin = [self.t_pos[i][u],self.h_pos[hn/2]] #X-vector
                        #end = [self.pos[1],self.pos[1]+1]
                        #self.draw_lines(1, begin,end, 1)
                        #hn += 1
                    
                    if lines_t[i][h1_ind][2] == None:
                        #print "    Draw external hole"
                        ncon(t_nodes[i][u],exit_nodes[i][u][1],order = 0, p_h =- 1)
                        #begin = [self.t_pos[i][u],self.t_pos[i][u]] #X-vector
                        #end = [self.pos[1],self.pos[1]+2]
                        #self.draw_lines(1, begin,end, -1)
                        #hn += 1
                    if lines_t[i][h1_ind][2] != None:
                        #print "    Draw internal hole"
                        
                        n = None
                        for t in range(len(h_nodes)):
                            if h_nodes[t].probe_line_a(-1):
                                #Found loop
                                h_nodes[t].connect_line_a(-1)
                                n = t
                                break
                        hn = h_nodes[n]
                        ncon(t_nodes[i][u],hn,order = 0, p_h = -1)
                        
                        
                        #begin = [self.t_pos[i][u],self.h_pos[hn/2]] #X-vector
                        #end = [self.pos[1],self.pos[1]+1]
                        #self.draw_lines(1, begin,end, -1)
                        #hn += 1
                    #hn += 1
        
        #Finally, we need to draw all lines exiting from the interaction
        for i in range(len(h_nodes)):
            if h_nodes[i].probe_line_c(1):
                #create particle leaving the interaction
                hn = h_nodes[i]
                ncon(hn, hnode([hn.x - scale*0.2, nc.y+1*scale],[],[]), order = 0, p_h = 1)
            if h_nodes[i].probe_line_c(-1):
                #create particle leaving the interaction
                hn = h_nodes[i]
                ncon(hn, hnode([hn.x + scale*0.2, nc.y+1*scale],[],[]), order = 0, p_h = -1)
        if self.vis:
            
            hold('off')
            show()


    def get_centered_pair(self, L, excluded):
        #pops centered line pairs from list (particle/hole)
        i1,i2 = None, None
        for i in range(len(L)-1):
            if L[i][0] in "abcdefgh" and L[i+1][0] in "ijklmnop":
                i1 = L[i]
                i2 = L[i+1]
                del(L[i+1])
                del(L[i])
                break              
        return i1, i2
        

    def draw_lines(self, order, begin, end, ph = 0):
        n1 = position(begin[0],begin[1])
        n2 = position(end[0], end[1])
        if order == 1:
            plot(begin, end, color = "black")
            half_x = begin[0] + (begin[1]-begin[0])/1.5
            half_y = end[0] + (end[1]-end[0])/1.5
            
            orient_x =  (begin[1]-begin[0])/2.0
            orient_y =  (end[1]-end[0])/2.0
            if ph == 1:
                draw_arrow([half_x, half_y],[orient_x, orient_y])
            if ph == -1:
                draw_arrow([half_x, half_y],[-orient_x, -orient_y])

        if order %2 == 0:
            #draw pairs of lines
            for i in range(order/2):
                ncon(n1,n2,-1,1)
                ncon(n2,n1,-1,-1)
            #print 0
        if order %2 != 0:
            #print "draw straight line and pairs"
            print ""
            
    def get_l_index(self, lines, target):
        ret = None
        for i in range(len(lines)):
            if lines[i][0] == target:
                ret = i
        return ret
    def get_t_index(self, To, target):
        ret = None
        for i in range(len(To)):
            for e in range(len(To[i])):
                if To[i][e] == target:
                    ret = [i,e]
        return ret[0], ret[1]
    
    def get_h_index(sefl,Ho,target):
        ret = None
        for i in range(len(Ho)):
            if Ho[i] == target:
                ret = i
        return ret

class position():
    def __init__(self, x, y):
        self.x = x
        self.y = y

class hnode():
    def __init__(self, pos, q_c, q_a):
        self.pos = pos
        self.x = pos[0]
        self.y = pos[1]
        self.q_c = q_c
        self.q_a = q_a
        self.setup()
    def setup(self):
        self.c_connections = [0,0]
        self.a_connections = [0,0]
        #particles, holes
        for i in self.q_c:
            if i == 1:
                self.c_connections[0] = None
            if i == -1:
                self.c_connections[1] = None
        for i in self.q_a:
            if i == 1:
                self.a_connections[0] = None
            if i == -1:
                self.a_connections[1] = None
    def probe_loop_a(self):
        #Probe annihilation for available loop connection
        ret = False
        if self.a_connections[0] == self.a_connections[1] == None:
            ret = True #A connection is possible
        return ret
    def probe_line_a(self, linetype):
        #Probe annihilation for 
        ret = False
        if linetype == -1:
            #scan for holw
            if self.a_connections[1] == None:
                ret = True
        if linetype == 1:
            #scan for particle
            if self.a_connections[0] == None:
                ret = True
        return ret
    def probe_line_c(self, linetype):
        #Probe annihilation for 
        ret = False
        if linetype == -1:
            #scan for holw
            if self.c_connections[1] == None:
                ret = True
        if linetype == 1:
            #scan for particle
            if self.c_connections[0] == None:
                ret = True
        return ret               

    def connect_loop_a(self):
        if self.probe_loop_a():
            self.a_connections[0] = 1
            self.a_connections[1] = 1
    def connect_line_a(self, linetype):
        if self.probe_line_a(linetype):
            if linetype == 1:
                self.a_connections[0] = 1
            if linetype == -1:
                self.a_connections[1] = 1
                
        
            
            
        
        
                
        

def nconnect(n1,n2,S,order="l0", p_h = None):
    N = 60

    Phx = (n1.x+n2.x)/2.0
    Phy = (n1.y+n2.y)/2.0
    
    lP = sqrt((n2.x-n1.x)**2 + (n2.y-n1.y)**2)
    dPx = (n2.x-n1.x)/lP    
    dPy = (n2.y-n1.y)/lP
    
    Cx = Phx - S*dPy
    Cy = Phy + S*dPx  
    
    
    lC = sqrt((S*dPy)**2 + (S*dPx)**2)
    #node(Phx,Phy, c="blue")
    #node(Cx,Cy, c="red")
    R = sqrt((Cx-n1.x)**2 + (Cy-n1.y)**2)

    lPC0 = sqrt(((Cx+R)-n1.x)**2 + (Cy - n1.y)**2)
    lPC1 = sqrt(((Cx+R)-n2.x)**2 + (Cy - n2.y)**2)

    dalpha = arccos((2*R**2 - lP**2)/(2.0*R**2))

    
    CPx = n1.x - Cx
    CPy = n1.y - Cy
    X,Y = 0,0
    if order == "0":
        #X = [n1.x, n2.x]
        #Y = [n1.y, n2.y]
        X = linspace(n1.x, n2.x, 6) #need to do this to get arrows working
        Y = linspace(n1.y, n2.y, 6)
    if order == "l0":
        if S<0:
            dalpha = 2*pi - dalpha
        A = linspace(0,dalpha, N)
        X,Y = rotate_v(CPx,CPy,A)
        X+=Cx
        Y+=Cy
    if order == "r0":
        if S>0:
            dalpha = 2*pi - dalpha      
        A = linspace(0,-dalpha, N)
        X,Y = rotate_v(CPx,CPy,A)
        X+=Cx
        Y+=Cy
        
    msize = 10
    if p_h == 1:
        draw_arrow([X[len(X)/2],Y[len(X)/2]], [-dPx,-dPy])
        #X[len(X)/2],Y[len(X)/2]
        #plot(X[len(X)/2],Y[len(X)/2], "^", color = "black", markersize = msize)
    if p_h == -1:
        draw_arrow([X[len(X)/2],Y[len(X)/2]], [dPx,dPy])
        #plot(X[len(X)/2],Y[len(X)/2], "v", color = "black", markersize = msize)
    plot(X,Y, color = "black")
    
def ncon(n1,n2,order = 0, p_h = None):
    if order == 0:
        nconnect(n1,n2,1,"0", p_h)
    if order > 0:
        nconnect(n1,n2,(-2+order),"l0", p_h)
    if order < 0:
        nconnect(n1,n2,(-2-order),"r0", p_h)

def draw_arrow(pos, point, s = .2, h = .1):
    #Draw an arrow at pos, pointing in the direction of point.
    #normalize direction
    p2 = sqrt(point[0]**2 + point[1]**2)
    point[0] /= p2
    point[1] /= p2
    
    #pi/2 degree rotation
    p_rotx, p_roty = point[1], -point[0]

    x0, y0 = pos[0], pos[1]
    x1, y1 = pos[0] - s*point[0], pos[1] - s*point[1]
    
    #plot the arrow
    plot([x0, x1+h*p_rotx],[y0, y1+h*p_roty], color = "black")
    plot([x0, x1-h*p_rotx],[y0, y1-h*p_roty], color = "black")
                
def rotate_v(x,y,alpha):
    ca = cos(alpha)
    sa = sin(alpha)
    return ca*x - sa*y, sa*x + ca*y            
    
        
def normal_ordered_hamiltonian():
    #Two particle hamiltonian
    F1 = Operator([1],[1])
    F2 = Operator([-1],[-1])
    F3 = Operator([],[1,-1])
    F4 = Operator([1,-1],[])
    
    V1 = Operator([1,1],[1,1])
    V2 = Operator([-1,-1],[-1,-1])
    V3 = Operator([1,-1],[1,-1])
    
    V4 = Operator([1],[1,1,-1])
    V5 = Operator([1,1,-1],[1])
    
    V7 = Operator([1,-1,-1],[-1])
    V6 = Operator([-1],[1,-1,-1])
    V9 = Operator([1,1,-1,-1],[])
    V8 = Operator([],[1,1,-1,-1])
    return [F1,F2,F3,F4,V1,V2,V3,V4,V5,V6,V7,V8,V9]                     
        
def cluster_operator(configuration):
    #configuration =[]
    T1= Operator([],[1,-1])
    T2= Operator([],[1,1,-1,-1])
    T3= Operator([],[1,1,1,-1,-1,-1])
    T4= Operator([],[1,1,1,1,-1,-1,-1,-1])
    T5= Operator([],[1,1,1,1,1,-1,-1,-1,-1,-1])
    T_list = [T1, T2, T3, T4, T5]
    ret = []
    for i in configuration:
        ret.append(T_list[i-1])
    return ret
                                    

class O():
    def __init__(self, H, T, name = "CC"):
        self.name = name
        self.H = H
        self.T = T
        self.diagrams = contract(H, T)
        self.Ndiag = len(self.diagrams)
        self.schematics = []
        self.stringforms = []
        self.cppcodes = []
        self.reports = []
        subindex = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","x","y","z"]
        for i in self.diagrams:
            self.schematics.append(schematic(i))
            self.schematics[-1].setup()
        self.N = len(self.schematics)
        for e in range(len(self.schematics)):
            
            c = self.schematics[e]
            self.reports.append(c.report())
            #name = self.name+"%i" %e
            sym = symbolic_diag(c,"%s" % (self.name+subindex[e]))
            self.stringforms.append(sym.stringform)
            CPPc = cpp_func(c,"%s" % (self.name+subindex[e]))
            self.cppcodes.append(CPPc.cppfunction)
    def E(self, i = None):
        E = self.H.E
        for i in range(len(self.T)):
            E += self.T[i].E
        return E
            
    def latex(self, i):
        ret = None
        try:
            ret =  self.stringforms[i]
        except:
            print "Invalid index given:", i
            print "There are only %i diagrams in this contraction." % self.N
        return ret
    def code(self, i):
        try:
            ret = self.cppcodes[i]
        except:
            print "Invalid index given:", i
            print "There are only %i diagrams in this contraction." % self.N
        return ret
    def diagram(self,i, pos = [0,0], vis = True):
        try:
            dia = feynman_diag(self.schematics[i], "D9", pos, vis)
        except:
            print "Invalid index given:", i
            print "There are only %i diagrams in this contraction." % self.N
    def report(self, i):
        ret = None
        try:
            ret = self.reports[i]
        except:
            print "Invalid index given:", i
            print "There are only %i reports in this contraction." % self.N
        return ret            

def list_mult(A, B):
    #multiply two lists
    ret = []
    for i in A:
        for e in B:
            ret.append(i+e)
    return ret


def list_pow(T, n):
    #Return list to the power of n
    th = T
    g = []
    for i in range(n-1):
        th = list_mult(th,T)
        g=th
    if n==1:
        g = T
    return g
            
def factorial(n):
    N = 1
    for i in range(1,n+1):
        N*=i
    return N

def expand_ansatz(T, n):
    #Taylor expand the exponential ansazt (exp(list(t-operatrs))) to a given truncation (n+1)
    expanded_eT = []
    prefactors  = []
    for i in range(1,n+1): #ignore i=0 -> uncorrelated energy
        elem = []
        prefactors.append(factorial(n)) # (**-1) inverted
        expT = list_pow(T,i)
        expanded_eT+=expT
    return expanded_eT

def combine_to_excitation(H,T,E, config = [1,0,0,0]):
    #Find combinations of operator elements in H, T that produce excitation level i
    skew = 0
    tex = []
    for h in range(len(H)):
        shift = 0
        for i in range(len(expT)):
            #print h, i
            stringname = "CC%i_%i" % (h,i)
            contr = O(H[h],expT[i], name = stringname)
            if contr.E() == E:
                for e in range(contr.Ndiag):
                    if config[0] == 1:
                        tx = contr.latex(e)
                        #Math(tx)
                        #print tx
                        tex.append(tx)
                    if config[1] == 1:
                        contr.diagram(e, [4*h,shift], True)
                    if config[1] == 2:
                        contr.diagram(e, [4*h,shift], False)
                    if config[3] == 1:
                        print contr.code(e)
                    if config[2] == 1:
                        print contr.reports[i]
                        
                    shift += 4
    if config[1] == 2:
        axes().set_aspect('equal', 'datalim')
        show()
        
    return tex
       
def combine_all(H,T, config = [1,0,0,0]):
    #Find combinations of operator elements in H, T that produce excitation level i
    
    for h in range(len(H)):
        shift = 0
        tex = []
        for i in range(len(expT)):
            #print h, i
            contr = O(H[h],expT[i])
            for e in range(contr.Ndiag):
                if config[0] == 1:
                    tx = contr.latex(e)
                    print tx
                    tex.append(tx)
                if config[0] == 2:
                    tx = contr.latex(e)
                    tex.append(tx)                    
                if config[1] == 1:
                    contr.diagram(e, [4*h,shift], True)
                if config[1] == 2:
                    contr.diagram(e, [4*h,shift], False)
                if config[3] == 1:
                    contr.code(e)
                    
                shift+= 4
    if config[1] ==2:
        axes().set_aspect('equal', 'datalim')
        show()
    return tex
    
T_2 = Operator([],[1,1,-1,-1])
H = normal_ordered_hamiltonian() #Including one- and two-particle interactions
expT = expand_ansatz([[T_2]],3)  #Taylor expand a list of lists to the 3rd order
contr = O(H[0],[T_2], name = "L")
print contr.report(0)
print contr.latex(0)
print "\n"

contr = O(H[4],[T_2], name = "L")
print contr.report(0)
print contr.latex(0)
print "\n"

contr = O(H[6],[T_2], name = "L")
print contr.report(0)
print contr.latex(0)
#print tx[3]