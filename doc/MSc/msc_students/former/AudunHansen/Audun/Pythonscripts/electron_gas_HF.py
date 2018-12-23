# -*- coding: utf-8 -*-
"""
Electron system interaction integrals
"""

from numpy import *

class electron_system():
    def __init__(self, L, N):
        self.L = L #size of volume is L^3
        self.N = int(sqrt(N) + 1)
        
        #working in atomic units
        self.e_charge = 1 
        self.hbar = 1
        self.m = 1
        self.generate_state_list() #generate a list of state in the order of increasing energy
        self.precompute_factors() #precompute common factors
    def generate_state_list(self):
        self.K = []
        self.S = []
        for x in range(self.N):
            for y in range(self.N):
                for z in range(self.N):
                    self.K.append([x,y,z])
                    self.K.append([x,y,z])
                    self.S.append(-1)
                    self.S.append(1)                    
        #print self.K
        
    def precompute_factors(self):
        self.p2factor = (4*pi*self.e_charge**2)/float(self.L)**3
        self.p1factor = (self.hbar**2)/float(2*self.m)
        self.efactor =  self.p1factor*(2*pi/float(self.L))**2
    def v(self,k_p, m_sp, k_q, m_sq, k_r, m_sr, k_s, m_ss):
        
        outer = self.p2factor*self.kd(self.sm(k_p, k_q),self.sm(k_r, k_s))
        inner1 = self.kd(m_sp, m_sr)*self.kd(m_sq, m_ss)*(1- self.kd(k_p, k_r))/self.absdiff2(k_r, k_p)
        inner2 = self.kd(m_sp, m_ss)*self.kd(m_sq, m_sr)*(1- self.kd(k_p, k_s))/self.absdiff2(k_s, k_p)
        return outer*inner1*inner2
        
    def f(self, k_p, m_sp, k_q, m_sq):
        #return fock matrix element
        return self.p1factor*self.kd(k_p, k_q)*self.kd(m_sp, m_sq)+self.p1AS(k_p, m_sp, k_q, m_sq)
    def p1AS(self, k_p, m_sp, k_q, m_sq):
        ret = 0
        for i in range(self.N):
            for s in [-1,1]:
               ret += self.v(k_p, m_sp, self.K[i], self.S[i], k_q, m_sq, self.K[i], self.S[i]) 
        return ret
    def sp_energy(self, nx,ny,nz):
        return self.efactor*(nx**2 + ny**2 + nz**2)
        
    def kd(self,p,q):
        #kroenecker delta function
        return (p==q)*1
    def absdiff2(self,k_p, k_q):
        k_ret = 0
        for i in range(len(k_p)):
            k_ret+= (k_p[i]-k_q[i])**2
        if k_ret == 0:
            k_ret = 1
        return k_ret
    def sm(self, k_p, k_q):
        sm = []
        for i in range(len(k_p)):
            sm.append(k_p[i] + k_q[i])
        return sm
        
el = electron_system(1.0, 100)
#print el.f([0,1,2],-1,[1,2,0],1)


print  el.v([0,1,2],-1,[1,2,0],1,[0,1,2],-1,[1,2,0],1)
#print "Test"
    