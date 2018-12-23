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
        self.Emax = array([0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 29, 30, 32, 33, 34, 35, 36, 37, 38, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 51, 52, 53, 54, 56, 57, 58, 59, 61, 62, 64, 65, 66, 67, 68, 69, 70, 72, 73, 74, 75, 76, 77, 78, 80, 81, 82, 83, 84, 85, 86, 88, 89, 90, 91, 93, 94, 96, 97, 98, 99, 100, 101, 102, 104, 105, 106, 107, 108, 109, 110, 113, 114, 115, 116, 117, 118, 120, 121, 122, 123, 125, 126, 128, 129, 130, 131, 132, 133, 134, 136, 137, 138, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149, 150, 152, 153, 154, 155, 157, 158, 160, 161, 162, 163, 164, 165, 166, 168, 169, 170, 171, 172, 173, 174, 176, 177, 178, 179, 180, 181, 182, 184, 185, 186, 187, 189, 190, 192, 193, 194, 195, 196, 197, 198, 200, 201, 202, 203, 204, 205, 206, 208, 209, 210, 211, 212, 213, 214, 216, 217, 218, 219, 221, 222, 224, 225, 226, 227, 228, 229, 230, 232, 233, 234, 235, 236, 237, 238, 241, 242, 243, 244, 245, 246, 248, 249, 250, 251, 253, 254, 256, 257, 258, 259, 260, 261, 262, 264, 265, 266, 267, 268, 269, 270, 272, 273, 274, 275, 276, 277, 278, 280, 281, 282, 283, 285, 286, 288, 289, 290, 291, 292, 293, 294, 296, 297, 298, 299, 300, 301, 302, 304, 305, 306, 307, 308, 309, 310, 312, 313, 314, 315, 317, 318, 320, 321, 322, 323, 324, 325, 326, 328, 329, 330, 331, 332, 333, 334, 336, 337, 338, 339, 340, 341, 342, 344, 345, 346, 347, 349, 350, 352, 353, 354, 355, 356, 357, 358, 360, 361, 362, 363, 364, 365, 366, 369, 370, 371, 372, 373, 374, 376, 377, 378, 379, 381, 382, 384, 385, 386, 387, 388, 389, 390, 392, 393, 394, 395, 396, 397, 398, 400, 401, 402, 403, 404, 405, 406, 408, 409, 410, 411, 413, 414, 416, 417, 418, 419, 420, 421, 422, 424, 425, 426, 427, 428, 429, 430, 432, 433, 434, 435, 436, 437, 438, 440, 441, 442, 443, 445, 446, 449, 450, 451, 452, 453, 454, 456, 457, 458, 459, 460, 461, 462, 464, 465, 466, 467, 468, 469, 470, 472, 473, 474, 475, 477, 478, 481, 482, 483, 484, 485, 486, 488, 489, 490, 491, 492, 493, 494, 497, 498, 499, 500, 501, 502, 504, 505, 506, 507, 509, 510, 512, 513, 514, 515, 516, 517, 518, 520, 521, 522, 523, 524, 525, 526, 528, 529, 530, 531, 532, 533, 534, 536, 537, 538, 539, 541, 542, 544, 545, 546, 548, 549, 550, 552, 553, 554, 555, 556, 557, 558, 561, 562, 563, 565, 566, 568, 569, 570, 571, 573, 574, 576, 577, 578, 579, 580, 581, 582, 584, 585, 586, 587, 588, 589, 590, 593, 594, 595, 596, 598, 601, 602, 603, 605, 606, 609, 611, 612, 613, 614, 616, 617, 618, 619, 620, 621, 622, 625, 626, 627, 629, 630, 633, 635, 637, 638, 641, 642, 644, 645, 646, 648, 649, 650, 651, 652, 653, 654, 656, 657, 658, 659, 661, 662, 664, 666, 667, 670, 673, 674, 675, 677, 678, 680, 681, 683, 684, 685, 686, 689, 693, 694, 697, 698, 699, 701, 706, 707, 708, 710, 712, 713, 714, 716, 717, 718, 721, 722, 723, 724, 726, 729, 730, 731, 734, 737, 738, 739, 741, 745, 747, 748, 749, 750, 753, 755, 757, 758, 761, 766, 768, 769, 770, 771, 774, 776, 782, 785, 786, 792, 794, 801, 803, 805, 806, 809, 811, 813, 817, 819, 822, 829, 834, 836, 838, 842, 843, 844, 846, 854, 866, 867, 869, 873, 875, 881, 891, 902, 904, 906, 910, 918, 937, 939, 941, 947, 972, 974, 978])[N]
        #Nm = N #int(sqrt(N) + 1)
        self.Nm = N
        Nm = N
        #Creating the basis

        for x in range(-Nm, Nm):
            for y in range(-Nm, Nm):
                for z in range(-Nm,Nm):
                    e = x*x + y*y + z*z
                    if e <=self.Emax:
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
        #self.Epphh = zeros((self.Np**2, self.Nh**2))
        
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
        #self.states[p,0]*(p==q)
        E = self.bs.states[:,0]
        pp = arange(Np**2)
        hh = arange(Nh**2)
        a = pp%Np + Nh
        b = pp//Np+ Nh
        i = hh%Nh
        j = hh//Nh
        #A = kron(ones((Nh**2,1)), E[a] + E[b]).T
        #I = kron(ones((Np**2,1)), E[i] + E[j])
        #self.Epppp = I-A
        self.IJ = E[i] + E[j]
        self.AB = E[a] + E[b]

        t0 = clock()
        self.B2 = optiV(self.bs)
        t1 = clock()
        
        
        
        print "Time spent setting up interactions:", t1-t0
        self.Vhhhh = self.B2.Vhhhh.tocsr()
        self.Vpppp = self.B2.Vpppp.tocsr()
        print "density:", len(self.Vpppp.data)/float(Np*Np*Np*Np)
        
        #self.lVpppp = self.B2.lVpppp #a block compressed representation of Vpppp
        
        self.Vhhpp = self.B2.Vhhpp.tocsr()
        self.Vpphh = self.B2.Vpphh.tocsr()
        self.Vhpph = self.B2.Vhpph.tocsr()
        #self.Vphhp = self.B2.Vhpph.T.tocsr()
        self.Tpphh = self.B2.Vpphh.tocsr()
        #self.Tpphh = csr_matrix(self.Vpphh.toarray()/self.Epppp)
        self.divide_eps(self.Tpphh)
        




        #Aligned matrices for L3, Q2, Q3 and Q4 multiplications
        #self.VL3 = self.perm_ind_ib_aj2ai_bj(self.Vhpph)   #old, working ?  
        
        self.VL3 = self.perm_ind_ia_bj2bi_aj(self.Vhpph)
            
        self.VQ2 = self.perm_ind_ij_ab2ai_bj(self.Vhhpp)
        self.VQ3 = self.perm_ind_ij_ba2iab_j(self.Vhhpp)
        self.VQ4 = self.perm_ind_ij_ba2bji_a(self.Vhhpp)
        #self.test_L3()
        #self.compare_L3VT()
    def sparsity(self):
        return len(self.Vpppp.data)/float(self.Np**4)
    def test_L3(self):
        Np = self.Np
        Nh = self.Nh
        V = self.Vhpph.toarray()
        T = self.Vpphh.toarray()
        VT = zeros((Np**2, Nh**2))
        self.sL3()
        self.PL3 = self.L3 - self.perm_ind_ba_ij(self.L3) - self.perm_ind_ab_ji(self.L3) + self.perm_ind_ba_ji(self.L3)
        #self.PQ2 = self.Q2 - self.perm_ind_ab_ji(self.Q2) 
        #self.PQ3 = self.Q3 - self.perm_ind_ab_ji(self.Q3)
        #self.PQ4 = self.Q4 - self.perm_ind_ba_ij(self.Q4)
        SQ2 = self.PL3.toarray()
        ret = True
        for a in range(Np):
            for b in range(Np):
                for i in range(Nh):
                    for j in range(Nh):
                        val = SQ2[a + b*Np, i+j*Nh]
                        s = 0
                        for k in range(Nh):
                            for c in range(Np):
                                #for c in range(Np):
                                #    for d in range(Np):
                                s += V[k + b*Nh, c + j*Np]*T[a + c*Np,i + k*Nh]
                                s -= V[k + a*Nh, c + j*Np]*T[b + c*Np,i + k*Nh]
                                s -= V[k + b*Nh, c + i*Np]*T[a + c*Np,j + k*Nh]
                                s += V[k + a*Nh, c + i*Np]*T[b + c*Np,j + k*Nh]

                        VT[a + b*Np, i+j*Nh] = s
                        #if VT[a + b*Np, i+j*Nh] != SQ2[a + b*Np, i+j*Nh]:
                        if abs(VT[a + b*Np, i+j*Nh] - SQ2[a + b*Np, i+j*Nh])>10e-25:
                            print VT[a + b*Np, i+j*Nh]
                            print SQ2[a + b*Np, i+j*Nh]
                            print VT[a + b*Np, i+j*Nh] - SQ2[a + b*Np, i+j*Nh]
                            print "=========================="
                            print "= Inconsistencies found in L3"
                            print "=========================="
                            print " "
                            #ret = False
                            break
                    if ret == False:
                        break
                        
                if ret == False:
                    break
                    
            if ret == False:
                break
               
        print "Done!" 


    def comparison(self, Z1,Z2):
        Nx = len(Z1)
        Ny = len(Z1[0])
        ret = True
        for i in range(Nx):
            for e in range(Ny):
                if Z1[i,e] != Z2[i,e]:
                    ret = False
                    break
            if ret == False:
                break
        return ret
    def compare_L3VT(self):
        self.sL3()
        self.old_L3()
        figure(1)
        imshow(self.L3_VT)
        show()
        figure(2)
        L3 = self.L3.toarray()
        imshow(L3)
        show()                
      
        Nh = self.Nh
        Np = self.Np
        ret = True
        for i in range(Nh):
            for j in range(Nh):
                for a in range(Np):
                    for b in range(Np):
                        if L3[a + b*Np, i +j*Nh] != self.L3_VT[a + b*Np, i + j*Nh]:
                            ret = False
                            break
                    if ret == False: break
                if ret == False: break
            if ret == False: break
        if ret == False: 
            print " "
            print "******************************"
            print "* Product failed consistency test."
            print "******************************"
            print " "
        else:
            print " "
            print "******************************"
            print "* Success: Product is consistent."
            print "******************************"
            print " "                                         
    
    def test_VL3(self):
        Nh = self.Nh
        Np = self.Np
        ret = True
        for i in range(Nh):
            for j in range(Nh):
                for a in range(Np):
                    for b in range(Np):
                        if self.VL3[b + i*Np, a + j*Np] != self.bs.v(b +Nh, i, a+Nh , j):
                            print self.VL3[b + i*Np, a + j*Np], self.bs.v(b +Nh, i, a+Nh , j)
                            ret = False
                            #break
                    #if ret == False: break
                #if ret == False: break
            #if ret == False: break
        if ret == False: 
            print " "
            print "******************************"
            print "* VL3 failed consistency test."
            print "******************************"
            print " "
        else:
            print " "
            print "******************************"
            print "* Success: VL3 is consistent."
            print "******************************"
            print " "
        

                            
        


        
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
        self.PQ2 = self.Q2 - self.perm_ind_ab_ji(self.Q2) #- self.perm_ind_ba_ij(self.Q2) + self.perm_ind_ba_ji(self.Q2)
        self.PQ3 = self.Q3 - self.perm_ind_ab_ji(self.Q3)
        self.PQ4 = self.Q4 - self.perm_ind_ba_ij(self.Q4)
        
        #Sum all contributions
        
        self.Tpphh = self.Vpphh + 0*.5*(self.L1 + self.L2) + self.PL3 + .25*self.Q1 + self.PQ2 - 0*.5*self.PQ3 - .5* self.PQ4 #/self.Epphh 
        #self.Tpphh = csr_matrix((self.Vpphh + .5*(self.L1 + self.L2) + self.PL3 + .25*self.Q1 + self.PQ2 - .5*self.PQ3 - .5* self.PQ4)/self.Epppp)
        self.divide_eps(self.Tpphh)
        #self.Tpphh = self.Tpphh #.toarray()
        #self.sp_epsdiv(self.Tpphh)

        #calculate energy
        self.energy_s()
        
        #Update UI
        print "        Correlation energy:", self.C_energy

        #Update amplitudes (have been temporarily dense due to division above)
        self.Tpphh = csr_matrix(self.Tpphh)
    def divide_eps(self, M):
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        
        E = self.bs.states[:,0]
        #M.data/=(self.IJ[cols] - self.AB[rows])
        M.data/=(E[i]+E[j] - E[a+self.Nh] - E[b+self.Nh])

        
    def e0_(self):
        Np = self.Np
        Nh = self.Nh
        Vhhpp = self.Vhhpp.toarray()
        Tpphh = self.Tpphh.toarray()
        e0 = 0.0
        for i in range(Nh):
            for j in range(Nh):
                for a in range(Np):
                    for b in range(Np):    
                        e0 += Vhhpp[i+j*Nh, a+b*Np]*Tpphh[a + b*Np, i+j*Nh]   
        return .25*e0     
    def energy(self):
        Np = self.Np
        Nh = self.Nh                    
        C = self.Vhhpp.dot(self.Tpphh)
        N = len(C)
        self.C_energy = .25*sum(C.diagonal())
        #self.C_energy = .25*sum(C[range(0,N), range(0,N)])
        
    def energy_s(self):
        Np = self.Np
        Nh = self.Nh                    
        C = self.Vhhpp.dot(self.Tpphh)
        cols, rows = C.indices, self.unpack_indptr(C.indptr)
        #N = len(C)
        self.C_energy = .25*sum(C.data[cols==rows])
        #self.C_energy = .25*sum(C[range(0,N), range(0,N)])    
        
    def sp_epsdiv(self, M):
        #sparse matrix energy division
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        a,b,i,j = rows%self.Np, rows//self.Np,cols%self.Nh, cols//self.Nh
        #print self.bs.states[:,0][i] + self.bs.states[:,0][j] - self.bs.states[:,0][a] - self.bs.states[:,0][b]
        M.data/=(self.bs.states[:,0][i] + self.bs.states[:,0][j] - self.bs.states[:,0][a+Nh] - self.bs.states[:,0][b+Nh])
        
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
        
    def perm_ind_ia_bj2bi_aj(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        i,a,b,j = rows%self.Nh, rows//self.Nh,cols%self.Np, cols//self.Np
        return coo_matrix((M.data, (b + i*self.Np, a + j*self.Np)), shape=(self.Np*self.Nh, self.Np*self.Nh)).tocsr()

    def perm_ind_ib_aj2ai_bj(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        i,b,a,j = rows%self.Nh, rows//self.Nh,cols%self.Np, cols//self.Np
        return coo_matrix((M.data, (a + i*self.Np, b + j*self.Np)), shape=(self.Np*self.Nh, self.Np*self.Nh)).tocsr()
    
    def perm_ind_ib_aj2bi_aj(self,M):
        #Sparse permutations
        cols, rows = M.indices, self.unpack_indptr(M.indptr)
        i,b,a,j = rows%self.Nh, rows//self.Nh,cols%self.Np, cols//self.Np
        return coo_matrix((M.data, (b + i*self.Np, a + j*self.Np)), shape=(self.Np*self.Nh, self.Np*self.Nh)).tocsr()
        
        
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
        self.L2 = self.Tpphh.dot(self.Vhhhh) 
        #self.L2 = (self.Vhhhh.T.dot(self.Tpphh.T)).T

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

class optiV():
    #quick initialization of interactions and amplitudes
    def __init__(self, bs):
        self.bs = bs
        self.Nh = bs.nparticles
        self.Np = bs.nstates-self.Nh
        self.Nm = 2*bs.Nm + 2
        self.Nm2 = self.Nm**2

        self.sVpppp2()
        self.sVhhhh()
        self.sVhhpp()
        self.sVhpph()
        self.sVpphh()

    def kdplt(self,oneN, states, P, Q, N):
        #kdplus, elementwise
        #oneNp = ones((Np**2, 1), dtype = int)
        return kron(oneN, states[:,N][P] + states[:,N][Q])
    
    def kdplq(self,states, P, Q, N):
        #kdplus, elementwise
        return states[:,N][P] + states[:,N][Q]

    def sVpppp2(self):
        #Setting up V^ab_cd, include spin considerations
        Np = self.Np
        Nh = self.Nh
        Nm = self.Nm+self.Nm + 2
        print Nm
        Nm2= Nm**2
        bs = self.bs
        
        #Setting up p,q indices
        AB = arange(Np**2)
        P = AB%Np + Nh
        Q = AB//Np+ Nh
        
        # kdplus = 4*pi/self.L3 #(kp[:,1]+kq[:,1]==kr[:,1]+ks[:,1])*(kp[:,2]+kq[:,2]==kr[:,2]+ks[:,2])*(kp[:,3]+kq[:,3]==kr[:,3]+ks[:,3])*4*pi/self.L3#d_k+k k+k

        idents_pp = bs.states[P,1:4]+bs.states[Q,1:4]
        #print idents_pp
        idents_pp = sum(array([1,Nm, Nm2]) * idents_pp, 1)
        
        self.idents_pp = idents_pp
        #print idents_pp
        #map unique outcomes
        uniques = unique(idents_pp)
        self.lVpppp = []
        T0 = array([], dtype = int)
        T1 = array([], dtype = int)
        data = array([])
        for e in uniques:
            T = where(idents_pp==e)[0]
            O = ones(len(T), dtype = int)
            
            t0 = kron(O, T)
            t1 = kron(T, O)
            

            a = t0%Np  + Nh
            b = t0//Np + Nh
            c = t1%Np  + Nh
            d = t1//Np + Nh
            
            a_ = bs.states[a,4] #spin
            b_ = bs.states[b,4]
            c_ = bs.states[c,4]
            d_ = bs.states[d,4]
            
            dpr_dqs = (a_==c_)*(b_==d_)
            dps_pqr = (a_==d_)*(b_==c_)

            dpr_dqs[dps_pqr] = True
            
            #dpr_dqs = (a_+b_) == (c_ + d_)
            
            #dpr_dqs[a==c] = False
            #dpr_dqs[a==d] = False
            

 
            t0 = t0[dpr_dqs]
            t1 = t1[dpr_dqs]
            """
            a = t0%Np  + Nh
            b = t0//Np + Nh
            c = t1%Np  + Nh
            d = t1//Np + Nh
            
            a_ = bs.states[a,:] #full
            b_ = bs.states[b,:]
            c_ = bs.states[c,:]
            d_ = bs.states[d,:]
            
            data_ = self.bs.V2(a_,b_,c_,d_)
            data = append(data, data_[data_!=0])
            """
            T0 = append(T0, t0) #[data_!=0])
            T1 = append(T1, t1) #[data_!=0])
            
        
        a = T0%Np  + Nh
        b = T0//Np + Nh
        c = T1%Np  + Nh
        d = T1//Np + Nh
            
        a_ = bs.states[a,:] #full
        b_ = bs.states[b,:]
        c_ = bs.states[c,:]
        d_ = bs.states[d,:]
            
        data = self.bs.V2(a_,b_,c_,d_)

        #print "Occupancy (nonzeros in pppp):", count_nonzero(data)/float(len(data))
        self.Vpppp = coo_matrix((data[data!=0], (T0[data!=0], T1[data!=0])), shape=(Np**2, Np**2))


    def sVpppp(self):
        #Setting up V^ab_cd, include spin considerations
        Np = self.Np
        Nh = self.Nh
        Nm = self.Nm
        Nm2= self.Nm2
        bs = self.bs
        
        #Setting up p,q indices
        AB = arange(Np**2)
        P = AB%Np + Nh
        Q = AB//Np+ Nh

        idents_pp = bs.states[P,1:4]+bs.states[Q,1:4]
        idents_pp = sum(array([1,Nm, Nm2]) * idents_pp, 1)
        
        self.idents_pp = idents_pp
        
        #map unique outcomes
        uniques = unique(idents_pp)
        self.lVpppp = []
        T0 = array([], dtype = int)
        T1 = array([], dtype = int)
        data = array([])
        for e in uniques:
            T = where(idents_pp==e)[0]
            O = ones(len(T), dtype = int)
            
            t0 = kron(O, T)
            t1 = kron(T, O)

            a = t0%Np  + Nh
            b = t0//Np + Nh
            c = t1%Np  + Nh
            d = t1//Np + Nh
            
            a_ = bs.states[a,4] #spin
            b_ = bs.states[b,4]
            c_ = bs.states[c,4]
            d_ = bs.states[d,4]
            
            #dpr_dqs = (a_==c_)*(b_==d_)
            #dps_pqr = (a_==d_)*(b_==c_)

            #dpr_dqs[dps_pqr] = True
            
            dpr_dqs = (a_+b_) == (c_ + d_)
            
            #dpr_dqs[a==c] = False
            #dpr_dqs[a==d] = False

 
            t0 = t0[dpr_dqs]
            t1 = t1[dpr_dqs]
            
            a = t0%Np  + Nh
            b = t0//Np + Nh
            c = t1%Np  + Nh
            d = t1//Np + Nh
            
            a_ = bs.states[a,:] #full
            b_ = bs.states[b,:]
            c_ = bs.states[c,:]
            d_ = bs.states[d,:]
            
            data_ = self.bs.V2(a_,b_,c_,d_)
            data = append(data, data_[data_!=0])
            T0 = append(T0, t0[data_!=0])
            T1 = append(T1, t1[data_!=0])

        #print "Occupancy (nonzeros in pppp):", count_nonzero(data)/float(len(data))
        self.Vpppp = coo_matrix((data, (T0, T1)), shape=(Np**2, Np**2))
        
    def sVhhhh(self):
        #Setting up V^ab_cd, include spin considerations
        Np = self.Np
        Nh = self.Nh
        Nm = self.Nm
        Nm2= self.Nm2
        bs = self.bs
        
        #Setting up p,q indices
        IJ = arange(Nh**2)
        P = IJ%Nh 
        Q = IJ//Nh

        idents_hh = bs.states[P,1:4]+bs.states[Q,1:4]
        idents_hh = sum(array([1,Nm, Nm2]) * idents_hh, 1)
        
        self.idents_hh = idents_hh
        
        #map unique outcomes
        uniques = unique(idents_hh)
        
        T0 = array([], dtype = int)
        T1 = array([], dtype = int)
        data = array([])
        for e in uniques:
            T = where(idents_hh==e)[0]
            O = ones(len(T), dtype = int)
            
            t0 = kron(O, T)
            t1 = kron(T, O)

            k = t0%Nh
            l = t0//Nh
            i = t1%Nh
            j = t1//Nh
            
            k_ = bs.states[k,4] #spin
            l_ = bs.states[l,4]
            i_ = bs.states[i,4]
            j_ = bs.states[j,4]
            
            dpr_dqs = (k_==i_)*(l_==j_)
            dps_pqr = (k_==j_)*(l_==i_)

            dpr_dqs[dps_pqr] = True

            t0 = t0[dpr_dqs]
            t1 = t1[dpr_dqs]
            
            k = t0%Nh
            l = t0//Nh
            i = t1%Nh
            j = t1//Nh
            
            k_ = bs.states[k,:] #full
            l_ = bs.states[l,:]
            i_ = bs.states[i,:]
            j_ = bs.states[j,:]
            
            data_ = self.bs.V2(k_,l_,i_,j_)
            data = append(data, data_[data_!=0])
            T0 = append(T0, t0[data_!=0])
            T1 = append(T1, t1[data_!=0])

        #print "Occupancy (nonzeros in hhhh):", count_nonzero(data)/float(len(data))
        self.Vhhhh = coo_matrix((data, (T0, T1)), shape=(Nh**2, Nh**2))                        


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
                """
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
                """

        i = T0%Nh
        j = T0//Nh 
        a = T1%Np + Nh
        b = T1//Np + Nh
        
        a_ = bs.states[a,:]
        b_ = bs.states[b,:]
        i_ = bs.states[i,:]
        j_ = bs.states[j,:]
        data = bs.V(i_.T,j_.T,a_.T,b_.T)
        #print "Occupancy (nonzeros in hhpp):", count_nonzero(data)/float(len(data))

        self.Vhhpp = coo_matrix((data[data!=0], (T0[data!=0], T1[data!=0])), shape=(Nh**2, Np**2))
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
                """   
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
                """

        a = T0%Np + Nh
        b = T0//Np + Nh
        i = T1%Nh 
        j = T1//Nh
        
        a_ = bs.states[a,:]
        b_ = bs.states[b,:]
        i_ = bs.states[i,:]
        j_ = bs.states[j,:]
        data = bs.V(a_.T,b_.T,i_.T,j_.T)
        #print "Occupancy (nonzeros in pphh):", count_nonzero(data)/float(len(data))

        self.Vpphh = coo_matrix((data[data!=0], (T0[data!=0], T1[data!=0])), shape=(Np**2, Nh**2))
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
                T0 = append(T0, t0)
                T1 = append(T1, t1)  
                """
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
        #print "Occupancy (nonzeros in hpph):", count_nonzero(data)/float(len(data))
        self.Vhpph = coo_matrix((data[data!=0], (T0[data!=0], T1[data!=0])), shape=(Np*Nh, Nh*Np))
        
        #return self.lVhpph

    def ident(self,v):
        #A unique identifying integer for the momentum combinations
        return v[0] + v[1]*self.bs.Nm + v[2]*self.bs.Nm**2 + v[3]*self.bs.Nm**3    
    def composite_dot(self, M1,M2):
        #A function for dotting a composite (blockwise) M1 with ordinary array M2
        Z = zeros((len(M2), len(M2[0])))
        for i in range(len(M1)):
            Z[M1[i][1]] = dot(M1[i][0], M2[M1[i][1]])
        return Z

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
                            #print M[i + j*Nh, k+l*Nh], bs.v(i , j ,k, l)
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
                            break
                    if ret == False:
                        break
                if ret == False:
                        break
            if ret == False:
                break
                        
        if ret == False:
            print "Found %i errors in hhpp" % n
        else:
            print "*************************"
            print "* Vhhpp is consistent!  *"
            print "*************************"

bs = electronbasis(3,1.0,2)
print bs.nstates
print "%.4f" % (.5*bs.hfenergy(14)/14.0)
print "1.9434"

"""

Q = CCD(bs)
N = 6
dens = zeros(N)
nstates = zeros(N)
for i in range(N):
    bs = electronbasis(i+2,1.0,2)
    print "%.4f" % (2*bs.hfenergy(14)/14)
    print "1.9434"
    
    Q = CCD(bs)
    dens[i] = Q.sparsity()
    nstates[i] = bs.nstates

from matplotlib.pyplot import *
print dens, nstates
plot(nstates, dens)
xlabel("Number of states [$N_s$]")
ylabel("Density [fraction of nonzero elements]")
title("Density of pp-pp interaction matrix")
show()

"""


"""
Ns = bs.nstates
Nh = bs.nparticles
print Ns, Nh
val = 0
for a in range(Nh,Ns):
    for b in range(Nh, Ns):
        for i in range(Nh):
            for j in range(Nh):
                val += bs.v(i,j,a,b)*bs.v(a,b,i,j)*.25/(bs.states[i][0] + bs.states[j][0] -bs.states[a][0] - bs.states[b][0])


print val/Nh
"""                


"""                                                                                                   
t0 = clock()                        
tb = electronbasis(4,1.0,14)
t1 = clock()

B2 = optiV(tb)
t2 = clock()
print "C:", t2-t1

B3 = consistency_test()
#B3.test_Vpppp(B2.Vpppp.toarray(), tb)
#B3.test_Vpphh(B2.Vpphh.toarray(), tb)
B3.test_Vhhpp(B2.Vhhpp.toarray(), tb)
#B3.test_Vhhhh(B2.Vhhhh.toarray(), tb)
#B3.test_Vhpph(B2.Vhpph.toarray(), tb)



print "Time spent on initializing basis:", t1-t0



print "====="
print "Number of states   :", tb.nstates
print "Number of particles:", tb.nparticles
print "====="
t0 = clock()
Q = CCD(tb)
t1 = clock()
print "Time spent on initializing solver:", t1-t0



#print "Old energy:", Q.e0_()
#print self.Tpphh.toarray().max()
Q.energy_s()
print "Energy    :", Q.C_energy

t0 = clock()
#print B2.Vhhpp

a,b,c,d = 33,34,28,5
print "---", Q.Vpppp[a + b*Q.Np, c + d*Q.Np], "---", tb.v(a+Q.Nh, b+Q.Nh, c+Q.Nh, d+Q.Nh)


for i in range(20):
    Q.advance()

#print abs(Q.C_energy+0.392696592145699 )
"""