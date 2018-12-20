from numpy import *

from CCD_HEG_Sparse_Implementation_2 import *

class optiV2():
    #quick initialization of interactions and amplitudes
    def __init__(self, bs):
        self.bs = bs
        self.Nh = bs.nparticles
        self.Np = bs.nstates-self.Nh
        self.Nm = 2*bs.Nm + 2
        self.Nm2 = self.Nm**2

        #self.sVpppp3()
        #self.sVhhhh()
        #self.sVhhpp()
        #self.sVhpph()
        #self.sVpphh()

    def kdplt(self,oneN, states, P, Q, N):
        #kdplus, elementwise
        #oneNp = ones((Np**2, 1), dtype = int)
        return kron(oneN, states[:,N][P] + states[:,N][Q])
    
    def kdplq(self,states, P, Q, N):
        #kdplus, elementwise
        return states[:,N][P] + states[:,N][Q]
    def sVpppp3(self):
        #Setting up V^ab_cd, include spin considerations
        Np = self.Np
        Nh = self.Nh
        Nm = self.Nm+self.Nm + 2
        #print Nm
        Nm2= Nm**2
        bs = self.bs
        
        #Setting up p,q indices
        #AB = arange(Np**2)
        #P = AB%Np + Nh
        #Q = AB//Np+ Nh
        
        aa = zeros(int(Np*((Np+1)/2.0)), dtype = int)
        bb = zeros(int(Np*((Np+1)/2.0)), dtype = int)
        n = 0
        for a in range(Np):
            for b in range(a,Np):
                aa[n] = a
                bb[n] = b
                n += 1
                
        
        
        # kdplus = 4*pi/self.L3 #(kp[:,1]+kq[:,1]==kr[:,1]+ks[:,1])*(kp[:,2]+kq[:,2]==kr[:,2]+ks[:,2])*(kp[:,3]+kq[:,3]==kr[:,3]+ks[:,3])*4*pi/self.L3#d_k+k k+k

        idents_pp = bs.states[aa+Nh,1:4]+bs.states[bb+Nh,1:4]
        #print idents_pp
        idents_pp = sum(array([1,Nm, Nm2]) * idents_pp, 1)
        
        self.idents_pp = idents_pp
        #print idents_pp
        #map unique outcomes
        uniques = unique(idents_pp)
        self.lVpppp = []
        T0 = array([], dtype = int)
        T1 = array([], dtype = int)
        va = array([], dtype = int)
        vb = array([], dtype = int)
        vc = array([], dtype = int)
        vd = array([], dtype = int)
        data = array([])
        for e in uniques:
            T = where(idents_pp==e)[0]
            O = ones(len(T), dtype = int)
            
            t0 = kron(O, T)
            t1 = kron(T, O)
            

            a = aa[t0] #+Nh
            b = bb[t0] #+Nh
            c = aa[t1] #+Nh
            d = bb[t1] #+Nh
            
            va = append(va,a)
            vb = append(vb,b)
            vc = append(vc,c)
            vd = append(vd,d)
            #a = t0%Np  + Nh
            #b = t0//Np + Nh
            #c = t1%Np  + Nh
            #d = t1//Np + Nh
            
            """
            a_ = bs.states[a+Nh,4] #spin
            b_ = bs.states[b+Nh,4]
            c_ = bs.states[c+Nh,4]
            d_ = bs.states[d+Nh,4]
            
            dpr_dqs = (a_==c_)*(b_==d_)
            dps_pqr = (a_==d_)*(b_==c_)

            dpr_dqs[dps_pqr] = True
            
            #dpr_dqs = (a_+b_) == (c_ + d_)
            
            #dpr_dqs[a==c] = False
            #dpr_dqs[a==d] = False
            """

 
            #t0 = t0[dpr_dqs]
            #t1 = t1[dpr_dqs]
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
            #T0 = append(T0, t0) #[data_!=0])
            #T1 = append(T1, t1) #[data_!=0])
            
            
        a_ = bs.states[va+Nh,:] #full
        b_ = bs.states[vb+Nh,:]
        c_ = bs.states[vc+Nh,:]
        d_ = bs.states[vd+Nh,:]
        
        #Use symmetries
        data = self.bs.V2(a_,b_,c_,d_)
        
        data = append(data, data)
        
        na = append(va,vb)
        nb = append(vb,va)
        nc = append(vc,vd)
        nd = append(vd,vc)
        
        #Once more use symmetries
        data = append(data,-data)
        
        a = append(na,nc)
        b = append(nb,nd)
        c = append(nc,nb)
        d = append(nd,na)
        
        
        T0 = a + b*Np
        T1 = c + d*Np

        #print "Occupancy (nonzeros in pppp):", count_nonzero(data)/float(len(data))
        #print len(data[data!=0])
        self.VppppQuick = coo_matrix((data[data!=0], (T0[data!=0], T1[data!=0])), shape=(Np**2, Np**2))

    def sVpppp2(self):
        #Setting up V^ab_cd, include spin considerations
        Np = self.Np
        Nh = self.Nh
        Nm = self.Nm+self.Nm + 2
        #print Nm
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
        print len(data[data!=0])
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


bs = electronbasis(2, 1.0, 14)
V = optiV2(bs)
t0 = clock()
V.sVpppp3()
#t1 = clock()
#V.sVpppp2()
#t2 = clock()
#print "With symmetries:", t1-t0
#print "No symmetries  :", t2-t1


tst = consistency_test()

tst.test_Vpppp(V.VppppQuick.toarray(), bs)

#V.sVpppp2()
#imshow(V.Vpppp.toarray()-V.VppppQuick.toarray())
#imshow(V.VppppQuick.toarray())
#show()

"""
Np = 38

N = int(Np*((Np+1)/2.0))

#N = Np*Np/2 + Np/2

print N, Np**2

va = zeros(N)
vb = zeros(N)
n = 0
for a in range(Np):
    #va[:-a] = 
    for b in range(a,Np):
        va[n] = a
        vb[n] = b
        n+=1

print va
print vb
"""