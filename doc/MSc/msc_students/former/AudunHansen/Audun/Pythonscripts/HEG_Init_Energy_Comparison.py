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
            self.states[i][0] *= 2*(pi**2)/self.L2 #Multiplying in the missing factors in the single particle energy
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



bs = electronbasis(2,1.0,2)

Ns = bs.nstates
Nh = bs.nparticles
n = 0
val = 0

def condition(a,b,i,j):
    ret = False
    kdx = (bs.states[a][1]+bs.states[b][1])**2
    kdy = (bs.states[a][2]+bs.states[b][2])**2
    kdz = (bs.states[a][3]+bs.states[b][3])**2
    kms = (bs.states[a][4]+bs.states[b][4])**2
    
    jdx = (bs.states[i][1]+bs.states[j][1])**2
    jdy = (bs.states[i][2]+bs.states[j][2])**2
    jdz = (bs.states[i][3]+bs.states[j][3])**2
    jms = (bs.states[i][4]+bs.states[j][4])**2
    
    if kdx+kdy+kdz+kms == jdx+jdy+jdz+jms:
        ret = True
    return ret
    
    
    
    
vals = zeros(((Ns-Nh)**2*Nh**2), dtype = float)
n = 0
for a in range(Nh,Ns):
    for b in range(Nh, Ns):
        for i in range(Nh):
            for j in range(Nh):
                vals[n] = bs.v(i,j,a,b)*bs.v(a,b,i,j)*.25/(bs.states[i][0] + bs.states[j][0] -bs.states[a][0] - bs.states[b][0])
                n += 1


print "Initialization energy   :", sum(vals)/Nh
print "Gustavs results         :", -0.0090095438138865
print "Absolute error          ;", (-0.00900954381388655 - sum(vals)/Nh)
print "Relative error to Gustav:", (-0.00900954381388655 - sum(vals)/Nh)/-0.00900954381388655
#bs.liststates()
      

