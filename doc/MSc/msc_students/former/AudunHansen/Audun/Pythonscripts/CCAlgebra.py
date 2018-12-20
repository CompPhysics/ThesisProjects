# <!-- collapse=True -->
from IPython.display import display, Math, Latex 

#from sympy.interactive import printing
#printing.init_printing()

from numpy import *
from itertools import *
from matplotlib.pyplot import *

class Operator():
    #Normal ordered operator for cluster algebra (diagrammatic)
    def __init__(self, q_create, q_annihilate):
        self.q_c = q_create
        self.q_a = q_annihilate
        self.diagrams = []
        self.contracted = []
        self.I = []
        self.T_operator = []
        self.T_vertices = []
        self.vertexform = []
        self.vertexform_locked = []
        self.labels_p = ['h','g','f','e','d','c','b','a']
        self.labels_h = ['p','o','n','m','l','k','j','i']
        self.enable_printing = False
        
        self.assess_excitation()
        
    def combine(self, ops, excitation = None):
        #Assosciate a list of T-operators (ops) to the current operator instance    
        #Find all possible ways to combine the operators using self.distinct_combinations() 
        T = []
        for i in ops:
            T.append(i.q_c)
            self.T_vertices.append(i.vertexform)
        self.T_operator = T
        
        #Finding acceptable combinations of internal contractions between the operators
        self.I = self.distinct_combinations(self.q_a, T)
        
        #Find excitation level of combination
        self.assess_excitation()
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
    
    def assess_contributions(self, excitation_level = None):
        #self.label_vertices()
        #returns the contribution to the CC-energy or amplitude eq. in symbolic form
        
        enable_printing = self.enable_printing #True #self.enable_printing

        
        self.loops = []
        self.holes = []
        self.equivalents = []
        self.stringforms = []
        for i in self.I:
            #for each distinct diagram, generate the expression for the contribution
            """
            Following S-B, ch. 10, this section evaluates the contribution to the CC-eqs. using the diagrammatic rules:
            (1)  Label internal and external lines with particle and hole labels
            (2)  Assosciate f(i,a) with every 1-p operator vertex
            (3)  Assosciate <lout rout!!lin rin> with every 2-p operator vertex
            (4)  Assosciate t(ab,ij) for every amplitude
            (5)  Sum over all internal lines
            (6)  Multiply by 1/2 for each pair of equivalent internal lines
            (7)  Multiply by 1/2 for each pair of equivalent T-vertices
            (8)  Include a phase factor (-)**(holes-loops)
            (9)  Not yet implemented
            (10) Not yet implemented
            """
            
            H_labels = []
            H_ins = []
            H_outs = []
            T_labels = []
            prefactor = 1.0 #this factor is multiplied to the sum and adjusted according to the rules laid out above
            
            #Comparing distribution of lines
            i_budget = self.itemcount(i) #internal distribution
            t_budget = self.itemcount(self.T_operator) #all lines from T
            t_budget_external = self.itemcount(self.T_operator) #distribution of external lines in T
            for e in range(len(t_budget_external)):
                #Subtracting internal lines 
                t_budget_external[e][0] -= i_budget[e][0]
                t_budget_external[e][1] -= i_budget[e][1]
            
            t_external_equiv = self.find_identical(t_budget_external) #len of this list yields number of identical external distribution in the Ts
            t_equiv = self.find_identical(i_budget) #len of this line yields the number of identical internal distributions in the Ts
            n_equivalent_t = 0
            for e in range(len(t_external_equiv)):
                if t_external_equiv[e] == t_equiv[e]:
                    n_equivalent_t += 1

            
            #Counting number of equivalent lines, holes and loops
            n_equi_lines = 0
            n_loops = 0
            n_loops_external = 0            
            n_holes = 0
            
            #Counting equivalent lines
            for e in i_budget:
                n_equi_lines += e[0]//2
                n_equi_lines += e[1]//2

            
            #Counting loops
            for e in range(len(i_budget)):
                n_loops += self.nloops(i_budget[e][0], i_budget[e][1])
                n_loops_external += self.nloops(t_budget_external[e][0], t_budget_external[e][1]) #So-called quasi-loops (S-B, ch. 10)

            
            #Counting holes
            for e in range(len(i_budget)):
                n_holes += i_budget[e][1] + t_budget_external[e][1]
            n_holes += self.q_c.count(-1)

            
            
            #labelling all lines
            
            #Setting up a mapping to keep track of connected lines in T
            t_mapping = []
            for e in range(len(self.T_operator)):
                t_mapping.append([])
                for u in range(len(self.T_operator[e])):
                    t_mapping[e].append(0) #a 0 implies an external line
            
            for e in range(len(i)):
                for u in range(len(i[e])):
                    linetype = i[e][u] #+1 = particle / -1 = hole
                    
                    for l in range(len(self.T_operator[e])):
                        if self.T_operator[e][l] == linetype:
                            if t_mapping[e][l] == 0:
                                #Connect and label line
                                t_mapping[e][l] = 1 #1 implies a connected line
                                break
                                

            #Labeling all looped lines:
            plabels = ['h','g','f','e','d','c','b','a']
            hlabels = ['p','o','n','m','l','k','j','i']
            
            i_connections = []
            i_labels = []
            for e in range(len(i)):
                i_connections.append([])
                i_labels.append([])                
                for u in range(len(i[e])):
                    i_connections[e].append(0)
                    i_labels[e].append(0)                    

            h_vertices = [] #list containing pairs of operators
            t_vertices = []
            for e in range(len(self.T_operator)):
                t_vertices.append([])

                
            #Connect all loops connecting at the hamiltonian
            for e in range(len(i)):
                for u in range(len(i[e])):
                    if i_connections[e][u] == 0:
                        linetype = i[e][u]
                        found = False
                        for ee in range(len(i)):
                            for uu in range(len(i[ee])):
                                if i_connections[e][u] == 0 and i_connections[ee][uu] == 0 and i[ee][uu] == -linetype:
                                    #Connect
                                    i_connections[e][u] = 2
                                    i_connections[ee][uu] = 2 #2 denotes a looped connection vertex
                                    if linetype == -1:
                                        i_labels[e][u] = hlabels.pop()
                                        i_labels[ee][uu] = plabels.pop()
                                        if e == ee:
                                            t_vertices[e].append([i_labels[e][u], i_labels[ee][uu]])
                                        if e != ee:
                                            t_vertices[e].append([i_labels[e][u], -1])
                                            t_vertices[ee].append([1,i_labels[ee][uu]])
                                    if linetype == 1:
                                        i_labels[e][u] = plabels.pop()
                                        i_labels[ee][uu] = hlabels.pop()
                                        if e == ee:
                                            t_vertices[e].append([i_labels[e][u], i_labels[ee][uu]])
                                        if e != ee:
                                            t_vertices[e].append([i_labels[e][u], -1])
                                            t_vertices[ee].append([1,i_labels[ee][uu]])
                                    #set in/out in vertex n
                                    h_vertices.append([i_labels[e][u],i_labels[ee][uu]])

                                        
                                    found = True
                                    break
                            if found:
                                break

            #Connect all lines passing through the hamiltonian
            for e in range(len(i)):
                for u in range(len(i[e])):
                    if i_connections[e][u] == 0:
                        linetype = i[e][u]
                        found = False
                        for ee in range(len(self.q_c)):
                            if self.q_c[ee] == linetype:
                                #add vertex, connect
                                i_connections[e][u] = 1
                                if linetype == 1:
                                    i_labels[e][u] = plabels.pop()
                                    t_vertices[e].append([i_labels[e][u], -1])
                                if linetype == -1:
                                    i_labels[e][u] = hlabels.pop()
                                    t_vertices[e].append([1, i_labels[e][u]])
                                h_vertices.append([i_labels[e][u], i_labels[e][u]])
                                
                                found = True
                            if found:
                                break
                                    
            #Label remaining lines in the T_operator
            for e in range(len(t_vertices)):
                for u in range(len(t_vertices[e])):
                    if t_vertices[e][u][0] == 1:
                       t_vertices[e][u][0] = plabels.pop()
                    if t_vertices[e][u][1] == -1:
                       t_vertices[e][u][1] = hlabels.pop()
            
            #Add unconnected vertices to T_operator
            for e in range(len(self.T_operator)):
                while len(self.T_operator[e])/2 > len(t_vertices[e]):
                    t_vertices[e].append([plabels.pop(), hlabels.pop()])
                                
                    
            stringform = '' #The stringform will contain the CC contribution in latex format         
            if enable_printing:                    
                print "=== Distinct contribution with excitation %i ===" % self.E
                print "Connection pattern:", i
            prefactor  = ((-1)**(n_holes - (n_loops+n_loops_external)))
            predivisor = (2**n_equi_lines)*(2**n_equivalent_t)
            if enable_printing:
                print "Multiplier:", prefactor, "/", predivisor

            stringform+=' \\frac{%i}{%i} ' % (prefactor, predivisor)
            summingover = ''
            
            summation_indices = []
            for e in range(len(h_vertices)):
                summation_indices.append(h_vertices[e][0])
                summation_indices.append(h_vertices[e][1])
            summation_indices2 = set(summation_indices) #only unique values
            #print summation_indices

            for e in summation_indices2:
                summingover+= e
            stringform+= '\sum_{%s} ' % summingover
            if enable_printing:
                print "Sum over   :", summation_indices
            
            H_tensor = ['','']
            
            for e in range(len(h_vertices)):
                H_tensor[1]+=h_vertices[e][0]
                H_tensor[0]+=h_vertices[e][1]
            if enable_printing:
                print "H-tensor   : <%s||%s> " % (H_tensor[0], H_tensor[1])
            
            stringform += '[%s|H|%s]' % (H_tensor[0], H_tensor[1])

            
            T_tensor = []
            for e in range(len(t_vertices)):
                T_tensor.append(['', ''])
                for u in range(len(t_vertices[e])):
                    T_tensor[e][1]+=t_vertices[e][u][0]
                    T_tensor[e][0]+=t_vertices[e][u][1]
            if enable_printing:
                print "T-tensor(s):"
            for e in T_tensor:
                if enable_printing:
                    print "             T(%s,%s)" % (e[0], e[1])
                stringform += 't_{%s}^{%s}  ' % (e[0], e[1])
            stringform += "      (excitation:%i) " % self.E
            if excitation_level == None:
                self.stringforms.append(stringform)
            else:
                if excitation_level == self.E:
                    self.stringforms.append(stringform)
                
                                
        if enable_printing:    
            print " "
        if len(self.I) == 0:
            #self.stringforms.append('No contributions found.')
            if enable_printing:
                print "=== No distinct contribution found ==="
                print " "
        #print self.I
        for e in range(len(self.I)):
            hh,tt = self.visualize(self.q_a, self.q_c, [0,e*2.6], self.I[e], T = self.T_operator, t = 0)
            self.diagrams.append([hh,tt])
        #self.hhvis, self.tvis = self.visualize(self.q_a, self.q_c, [0,0], self.I[0], T = None, t = 0)
        return 0

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
            
            
    
    def distinct_combinations(self,H,T):
        #Returns all possible combinations of H and T
        #I - all q-particle annihilation operators
        #T list of list with T operators. ex. [[-1,1],[-1,1]]  = T_1 T_1
       lenH  = len(H)
       lenT  = len(T)
       lenTi = []
       for i in range(lenT):
           lenTi.append(len(T[i]))
       H+=[0 for i in range(lenT-1)] #adding zeros to denote separations in the cluster-operators
       
       #Creating countlist for T-operator to keep track of q-operators in each clusteroperator
       T_budget = self.itemcount(T)
       
       
       #Create all permutations
       H_permuted = permutations(H)
       
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
    def printout(self):
        for i in self.stringforms:
            display(Math(r'%s' % i))
    def plot_diagrams(self):
        figure(figsize = (2, len(self.diagrams)+1), dpi = 80, edgecolor = 'k',facecolor='white')
        p2 =0.0
        for i in self.diagrams:
            pos2 = [0,p2]
            for e in i[0]:
                e.draw()
            for e in i[1]:
                for u in e:
                    u.draw()
            p2 += 1.8
        #axisbg='red'
        #set_cmap('hot')
        axis('off')
        axes().set_aspect('equal', 'datalim')
        show()


    def visualize(self,h_below, h_above, pos, I, T = None, t = 0):
        #create operator vxnode objects from Ob-object
        #NV = len(O.L)/2
        NV = (len(h_below) + len(h_above))/2
        Nbelow = len(h_below)
        Nabove = len(h_above)
        
        c_below = []
        for i in range(Nbelow):
            c_below.append(0) #A zero implies a "free" line below the interaction line
            
        c_above = []
        for i in range(Nabove):
            c_above.append(0) #A zero implies a "free" line above the interaction line
        ncount = NV
        #(1) Identify lines passing through the interaction line
        vnodes = []
        for i in range(Nabove):
            for e in range(Nbelow):
                if c_below[e] == 0 and c_above[i] == 0:
                    if h_above[i]== h_below[e]:
                        #Append the operator to vnodes
                        #print "Found a line passing through the interaction."
                        c_above[i] = 1
                        c_below[e] = 1
                        ncount -= 1
                        if h_above[i] == 1:
                            c1,c2,c3,c4 = [0,None],[1,None],[0,None],[1,None]
                            nd = vxnode([pos[0] + i, pos[1]+1.5], [c1,c2,c3,c4])
                            vnodes.append(nd)
                                            
                        if h_above[i] == -1:
                            c1,c2,c3,c4 = [1,None],[0,None],[1,None],[0,None]
                            nd = vxnode([pos[0] + i, pos[1]+1.5], [c1,c2,c3,c4])
                            vnodes.append(nd)
    
                    
        
        #(2) Identify lines annihilating at the interaction 
        for i in range(Nbelow):
            for e in range(Nbelow):
                if c_below[e] == 0 and c_below[i] == 0 and i!= e:
                    if h_below[i]== -1*h_below[e]:
                        #Append the operator to vnodes
                        #print "Found a line annihilating at the interaction."
                        c_below[i] = 1
                        c_below[e] = 1
                        ncount -= 1
                        c1,c2,c3,c4 = [0,None],[0,None],[1,None],[1,None]
                        nd = vxnode([pos[0] + i, pos[1]+1.5], [c1,c2,c3,c4])
                        vnodes.append(nd)
                
        #(3) Identify lines created at the interaction 
        for i in range(Nabove):
            for e in range(Nabove):
                if c_above[e] == 0 and c_above[i] == 0 and i!= e:
                    if h_above[i]== -1*h_above[e]:
                        #Append the operator to vnodes
                        #print "Found a line creating at the interaction"
                        c_above[i] = 1
                        c_above[e] = 1
                        ncount -= 1
                        c1,c2,c3,c4 = [1,None],[1,None],[0,None],[0,None]
                        nd = vxnode([pos[0] + i, pos[1]+1.5], [c1,c2,c3,c4])
                        vnodes.append(nd)
        #print "C_above:", c_above
        #print "C_below:", c_below
    
        
        if NV%2 != 0:
            print "Warning: non-binary operator."
    
    
        for i in range(len(vnodes)-1):
            if t == 0:
                vnodes[i].opconnect(vnodes[i+1].pos)
            if t == 1:
                vnodes[i].tconnect(vnodes[i+1].pos)
            
        tnodes = []
        #print "lenT:", len(T), T
        
        if T != None:
            
            for t in range(len(T)):
                tnodes.append([])
                
                for i in range(len(T[t])/2):
                    c1,c2,c3,c4 = [1,None],[1,None],[0,None],[0,None]
                    #nd = vxnode([pos[0] + i, pos[1]], [c1,c2,c3,c4])
                    tnodes[t].append(vxnode([pos[0] + i, pos[1]], [c1,c2,c3,c4]))
                        
                for i in range(len(tnodes[t])-1):             
                    tnodes[t][i].tconnect(tnodes[t][i+1].pos)
                
            p = 0
            for i in range(len(tnodes)):
                #print "Tlen:", len(tnodes[i])
                for e in range(len(tnodes[i])):
                    tnodes[i][e].pos[0] = pos[0] + p
                    p += 1
            for i in range(len((vnodes))):
                vnodes[i].pos[0] += .5


        
        #Contract  T
        #I contains a recipe for the contractions in the diagram.
        #Iterate over each element in I and match up the contractions

        #For every element in I, match up corresponding elements in H and T

        for i in range(len(tnodes)):

            for e in range(len(I[i])):

                cn = I[i][e]

                self.m = 0
                nn = 0
                cond = True
                while cond:
                    cond = self.probe(tnodes[i],vnodes,cn) #probe and perform a possible connection

                    nn += 1
                    if self.m != 0:

                        break
                    if nn>10:
                        break
        for i in range(len(tnodes)):
            pass              
        return [vnodes, tnodes]
            
    
    def probe(self,T,H,cn):
        cond = True
        if cn == -1:
            #hole line
            for i in range(len(T)):
                for e in range(len(H)):
                    if T[i].config[0][0] == 1 and T[i].config[0][1] == None:
                        if H[e].config[2][0] == 1 and H[e].config[2][1] == None:
                            #perform connection
                            H[e].config[2][1] = T[i].pos
                            T[i].config[0][0] = 0
                            self.m = 1
                            cn = 0
                            break
                if self.m == 1:
                    break
                            
        if cn == 1:
            #particle line
            for i in range(len(T)):
                for e in range(len(H)):

                    if T[i].config[1][0] == 1 and T[i].config[1][1] == None:
                        if H[e].config[3][0] == 1 and H[e].config[3][1] == None:
                            #perform connection
                            H[e].config[3][1] = T[i].pos
                            T[i].config[1][0] = 0
                        
                            self.m = 1

                            cn = 0
                            
                            break
                if self.m == 1:
                    break
        if self.m == 1:
            cond = False
        return cond   

    def nconnect(self,n1,n2,S,order="l0", p_h = None):
        N = 60
    
        if n1.x==n2.x and n1.y == n2.y:
            Cx = n1.x + S
            Cy = n1.y
            
            X = Cx + S*cos(linspace(0,2*pi,N))
            Y = Cy + S*sin(linspace(0,2*pi,N))
        #S = -1
        else:
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
                X = [n1.x, n2.x]
                Y = [n1.y, n2.y]
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
            draw_arrow([X[len(X)/2],Y[len(X)/2]], [-dPx,-dPy])
            #plot(X[len(X)/2],Y[len(X)/2], "v", color = "black", markersize = msize)
        plot(X,Y, color = "black")
    
    def ncon(self,n1,n2,order = 0, p_h = None):
        if order == 0:
            nconnect(n1,n2,1,"0", p_h)
        if order > 0:
            nconnect(n1,n2,(-2+order),"l0", p_h)
        if order < 0:
            nconnect(n1,n2,(-2-order),"r0", p_h)
    
    
def rotate_v(x,y,alpha):
    ca = cos(alpha)
    sa = sin(alpha)
    return ca*x - sa*y, sa*x + ca*y
    

class vxnode():
    def __init__(self, pos, config):
        #config = [[1,None],[1,None],[1,None],[1,None]]
        self.pos = pos
        self.hole_up = 0  #config[0] #Outgoing Q-particle creation operators
        self.part_up = 0  #config[1] #Outgoing Q-particle creation operators
        self.hole_down = 0#config[2] #Outgoing Q-particle creation operators
        self.part_down = 0#config[3] #Outgoing Q-particle creation operators
        self.config = config
        self.subline = False
        
        self.Opconnect = []
        self.Tconnect = []
        self.vconnect_h = []
        self.c = "black"
    def opconnect(self, pos):
        #connect horizontal to another vxnode
        self.Opconnect.append(pos)
    def tconnect(self, pos):
        #connect horizontal to another vxnode
        self.Tconnect.append(pos)        
    def ttype(self):
        #Draw a solid, horizontal line through the operator
        self.subline = True
    def draw(self, pos2 = None):
        if pos2 != None:
            self.pos[0] += pos2[0]
            self.pos[1] += pos2[1]            
            
        msize = 10
        c= self.c
        hold('on')
        sx = .4
        sy = 2.0
        #if len(self.Tconnect) != 0:
            #sy *= 1.5
            #sx *= 1.3
        config = self.config
        plot(self.pos[0],self.pos[1], ".", color = c,markersize = 10)
        if config[0][0] == 1:
            if config[0][1] == None:
                #Draw straight line hole up
                plot([self.pos[0], self.pos[0]-sx],[self.pos[1], self.pos[1]+sy],color = c)

                #print "DRAWING ARROW"
                draw_arrow([(self.pos[0] + self.pos[0]-sx)/2.0,(self.pos[1] + self.pos[1]+sy)/2.0], [sx,-sy])
                
                #plot((self.pos[0] + self.pos[0]-sx)/2.0,(self.pos[1] + self.pos[1]+sy)/2.0,"v", color = c, markersize = msize)                
            if config[0][1] != None:
                #connect to node config[0][1] in hole up manner
                order = -1
                #if config[1][1] != None:
                #    order = 0
                self.ncon(node(self.pos),node(config[0][1]), order, -1)
                
        if config[1][0] == 1:
            if config[1][1] == None:
                #Draw straight line particle up
                plot([self.pos[0], self.pos[0]+sx],[self.pos[1], self.pos[1]+sy],color = c)
                
                #print "DRAWING ARROW"
                draw_arrow([(self.pos[0] + self.pos[0]+sx)/2.0,(self.pos[1] + self.pos[1]+sy)/2.0], [sx,sy])
                #plot((self.pos[0] + self.pos[0]+sx)/2.0,(self.pos[1] + self.pos[1]+sy)/2.0,"^", color = c, markersize = msize)                
            if config[1][1] != None:
                #connect to node config[0][1] in particle up manner
                order = -1
                #if config[0][1] != None:
                #    order = 0
                self.ncon(node(config[1][1]),node(self.pos),order,1)
                
        if config[2][0] == 1:
            if config[2][1] == None:
                #Draw straight line hole down
                plot([self.pos[0], self.pos[0]-sx],[self.pos[1], self.pos[1]-sy],color = c)
                
                #print "DRAWING ARROW"
                draw_arrow([(self.pos[0] + self.pos[0]-sx)/2.0,(self.pos[1] + self.pos[1]-sy)/2.0], [-sx,-sy])
                #plot((self.pos[0] + self.pos[0]-sx)/2.0,(self.pos[1] + self.pos[1]-sy)/2.0,"v", color = c, markersize = msize)                
            if config[2][1] != None:
                #connect to node config[0][1] in hole down manner
                #print "Active"
                order = -1
                #if config[3][1] != None:
                #    order = 0
                self.ncon(node(config[2][1]), node(self.pos),order, -1)
                
        if config[3][0] == 1:
            if config[3][1] == None:
                plot([self.pos[0], self.pos[0]+sx],[self.pos[1], self.pos[1]-sy],color = c)
                
                #print "DRAWING ARROW"
                draw_arrow([(self.pos[0] + self.pos[0]+sx)/2.0,(self.pos[1] + self.pos[1]-sy)/2.0], [-sx,sy])
                #plot((self.pos[0] + self.pos[0]+sx)/2.0,(self.pos[1] + self.pos[1]-sy)/2.0,"^", color = c, markersize = msize)
                #Draw straight line particle down

            if config[3][1] != None:
                #connect to node config[0][1] in particle down manner
                order = -1
                #if config[2][1] != None:
                #    order = 0
                self.ncon(node(self.pos), node(config[3][1]), order, 1)
        for i in range(len(self.Opconnect)):
            plot([self.pos[0],self.Opconnect[i][0]],[self.pos[1],self.Opconnect[i][1]], ls = "dotted", color = c)
        
        for i in range(len(self.Tconnect)):
            plot([self.pos[0],self.Tconnect[i][0]],[self.pos[1],self.Tconnect[i][1]], color = c)
            
        if self.subline:
            plot([self.pos[0]-sx, self.pos[0]+sx], [self.pos[1], self.pos[1]], color = c)    
    def ncon(self,n1,n2,order = 0, p_h = None):
        if order == 0:
            self.nconnect(n1,n2,1,"0", p_h)
        if order > 0:
            self.nconnect(n1,n2,(-2+order),"l0", p_h)
        if order < 0:
            self.nconnect(n1,n2,(-2-order),"r0", p_h)
    def nconnect(self,n1,n2,S,order="l0", p_h = None):
        N = 60
    
        if n1.x==n2.x and n1.y == n2.y:
            Cx = n1.x + S
            Cy = n1.y
            
            X = Cx + S*cos(linspace(0,2*pi,N))
            Y = Cy + S*sin(linspace(0,2*pi,N))
        #S = -1
        else:
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
                X = [n1.x, n2.x]
                Y = [n1.y, n2.y]
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
            draw_arrow([X[len(X)/2],Y[len(X)/2]], [-dPx,-dPy])
            #plot(X[len(X)/2],Y[len(X)/2], "v", color = "black", markersize = msize)
        plot(X,Y, color = "black")

def draw_arrow(pos, point, s = .2, h = .1):
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

class node():
    def __init__(self, V, c= "black"):
        self.x = V[0]
        self.y = V[1]
        #plot(x,y,".", color = c,markersize = 15)
                
def normal_ordered_hamiltonian():
    F1 = Operator([1],[1])
    F2 = Operator([-1],[-1])
    F3 = Operator([1,-1],[])
    F4 = Operator([],[1,-1])
    
    V1 = Operator([1,1],[1,1])
    V2 = Operator([-1,-1],[-1,-1])
    V3 = Operator([1,-1],[1,-1])
    
    V4 = Operator([1,1,-1],[1])
    V5 = Operator([1],[1,1,-1])
    
    V6 = Operator([1,-1,-1],[-1])
    V7 = Operator([-1],[1,-1,-1])
    V8 = Operator([1,1,-1,-1],[])
    V9 = Operator([],[1,1,-1,-1])
    return [F1,F2,F3,F4,V1,V2,V3,V4,V5,V6,V7,V8,V9]

def cluster_operator(configuration):
    #configuration =[]
    T1= Operator([1,-1],[])
    T2= Operator([1,1,-1,-1],[])
    T3= Operator([1,1,1,-1,-1,-1],[])
    T4= Operator([1,1,1,1,-1,-1,-1,-1],[])
    T5= Operator([1,1,1,1,1,-1,-1,-1,-1,-1],[])
    T_list = [T1, T2, T3, T4, T5]
    ret = []
    for i in configuration:
        ret.append(T_list[i-1])
    return ret

def generate_all_combinations(H,T, excitation_level = None, printing = 0):
    for i in H:
        i.combine(T)
        i.assess_contributions(excitation_level)
        i.printout()


T = cluster_operator([2,1])
H = normal_ordered_hamiltonian()

generate_all_combinations(H,T,0)