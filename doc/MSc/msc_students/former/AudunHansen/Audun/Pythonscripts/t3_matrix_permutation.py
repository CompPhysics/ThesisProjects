from numpy import *
from matplotlib.pyplot import *
#shuffling matrices

Na = 3
Nb = 2
Nc = 4

Ni = 5
Nj = 2
Nk = 1

Nrows = Na*Nb*Nc
Ncols = Ni*Nj*Nk

Z = zeros((Nrows, Ncols), dtype = int)
for a in range(Na):
    for b in range(Nb):
        for c in range(Nc):
            for i in range(Ni):
                for j in range(Nj):
                    for k in range(Nk):
                        Z[a +b*Na + c*Na*Nb, i + j*Ni + k*Ni*Nk] = a +b*Na + c*Na*Nb+ Na*Nb*Nc*(i + j*Ni + k*Ni*Nk)

#print Z
"""
def permute(Ntrue, Nperm, Z):
    Z2 = Z*0
    Na, Nb, Nc, Ni, Nj, Nk = Nperm
    Np, Nq, Nr, Ns, Nt, Nu = Ntrue
    
    for a in range(Np):
        for b in range(Nq):
            for c in range(Nr):
                for i in range(Ns):
                    for j in range(Nt):
                        for k in range(Nu):
                            Z2[c +b*Np + a*Np*Nq, i + j*Ns + k*Ns*Nt] = Z[a +b*Na + c*Na*Nb,i + j*Ni + k*Ni*Nk] 
    return Z2
"""
    
def permute(Ntrue, Nperm, Z):
    Z2 = Z*0
    #Nperm = [0,1,2,3,4,5]
    
    Na, Nb, Nc, Ni, Nj, Nk = Ntrue
    Np, Nq, Nr, Ns, Nt, Nu = Ntrue[Nperm]
    print Ntrue
    print Ntrue[Nperm]
    for a in range(Na):
        for b in range(Nb):
            for c in range(Nc):
                for i in range(Ni):
                    for j in range(Nj):
                        for k in range(Nk):
                            p,q,r,s,t,u = array([a,b,c,i,j,k])[Nperm]
                            Z2[p +q*Np + r*Np*Nq, s + t*Ns + u*Ns*Nt] = Z[a +b*Na + c*Na*Nb,i + j*Ni + k*Ni*Nk] 
    return Z2    
    

"""
Z2 = zeros((Nrows, Ncols), dtype = int)
for a in range(Na):
    for b in range(Nb):
        for c in range(Nc):
            for i in range(Ni):
                for j in range(Nj):
                    for k in range(Nk):
                        Z2[c +b*Nc + a*Nb*Nc, i + j*Ni + k*Ni*Nk] = Z[a +b*Na + c*Na*Nb,i + j*Ni + k*Ni*Nk]
"""
Ns = array([Na,Nb,Nc,Ni,Nj,Nk])
Zab = permute(Ns ,array([1,0,2,3,4,5], dtype = int), Z)     
Zac = permute(Ns ,array([2,1,0,3,4,5], dtype = int), Z)     
Zij = permute(Ns ,array([1,0,2,4,3,5], dtype = int), Z)     
Zik = permute(Ns ,array([1,0,2,5,4,3], dtype = int), Z)  
Zacik = permute(Ns ,array([2,1,0,5,4,3], dtype = int), Z)
lv = linspace(Z.min(), Z.max(), 300)
         
figure(1)
title("Permuting blocks in the T3-amplitude")
subplot(3,2,1)
title("No permutation")
axis("off")
contourf(Z, levels = lv, cmap = "RdGy")

subplot(3,2,2)
axis("off")
title("P(ab)")
contourf(Zab, levels = lv, cmap = "RdGy")

subplot(3,2,3)
axis("off")
title("P(ac)")
contourf(Zac, levels = lv, cmap = "RdGy")

subplot(3,2,4)
axis("off")
title("P(ij)")
contourf(Zij, levels = lv, cmap = "RdGy")

subplot(3,2,5)
axis("off")
title("P(ik)")
contourf(Zik, levels = lv, cmap = "RdGy")

subplot(3,2,6)
axis("off")
title("P(ac)P(ik)")
contourf(Zacik, levels = lv, cmap = "RdGy")

show()


