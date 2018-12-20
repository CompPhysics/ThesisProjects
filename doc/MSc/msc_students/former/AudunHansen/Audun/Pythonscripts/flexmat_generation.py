
from itertools import *

def constr22(L, R):
    l = L
    r = R
    string = ""
    string += "sp_mat flexmat::%s_%s(){\n" %(l,r)
    string += "    if(N%s_%s == 0){\n" %(l,r)
    string += "        locations.set_size(vp.size(), 2);\n"
    string += "        locations.col(0) = v%s + v%s*iN%s;\n" % (L[0], L[1], L[0])
    string += "        locations.col(1) = v%s + v%s*iN%s;\n" % (R[0], R[1], R[0])
    string += "        V%s_%s = sp_mat(locations.t(), vValues, iN%s*iN%s, iN%s*iN%s);\n" % (l,r,L[0], L[1], R[0], R[1])
    string += "        N%s_%s = 1;\n" % (l,r)
    string += "        return V%s_%s;\n" % (l,r)
    string += "    }\n"
    string += "    else{\n"
    string += "        return V%s_%s;\n" % (l,r)
    string += "    }\n"
    string += "}\n"
    return string

def header22(L,R):
    string = ""
    string += "    int N%s_%s;\n" % (L,R)
    string += "    sp_mat V%s_%s;\n" % (L,R)
    string += "    sp_mat %s_%s();\n" %(L,R)
    return string
    
def constr13(L, R):
    l = L
    r = R
    string = ""
    string += "sp_mat flexmat::%s_%s(){\n" %(l,r)
    string += "    if(N%s_%s == 0){\n" %(l,r)
    string += "        locations.set_size(vp.size(), 2);\n"
    string += "        locations.col(0) = v%s;\n" % (L[0])
    string += "        locations.col(1) = v%s + v%s*iN%s + v%s*iN%s*iN%s;\n" % (R[0], R[1], R[0], R[2], R[0], R[1])
    string += "        V%s_%s = sp_mat(locations.t(), vValues, iN%s, iN%s*iN%s*iN%s);\n" % (l,r,L[0], R[0], R[1], R[2])
    string += "        N%s_%s = 1;\n" % (l,r)
    string += "        return V%s_%s;\n" % (l,r)
    string += "    }\n"
    string += "    else{\n"
    string += "        return V%s_%s;\n" % (l,r)
    string += "    }\n"
    string += "}\n"
    return string

def header13(L,R):
    string = ""
    string += "    int N%s_%s;\n" % (L,R)
    string += "    sp_mat V%s_%s;\n" % (L,R)
    string += "    sp_mat %s_%s();\n" %(L,R)
    return string    
    
def constr31(L, R):
    l = L
    r = R
    string = "\n"
    string += "sp_mat flexmat::%s_%s(){\n" %(l,r)
    string += "    if(N%s_%s == 0){\n" %(l,r)
    string += "        locations.set_size(vp.size(), 2);\n"
    string += "        locations.col(0) = v%s + v%s*iN%s + v%s*iN%s*iN%s;\n" % (L[0], L[1], L[0], L[2], L[0], L[1])
    string += "        locations.col(1) = v%s;\n" % (R[0])
    string += "        V%s_%s = sp_mat(locations.t(), vValues, iN%s*iN%s*iN%s,iN%s);\n" % (l,r,L[0], L[1], L[2], R[0])
    string += "        N%s_%s = 1;\n" % (l,r)
    string += "        return V%s_%s;\n" % (l,r)
    string += "    }\n"
    string += "    else{\n"
    string += "        return V%s_%s;\n" % (l,r)
    string += "    }\n"
    string += "}\n"
    return string

def header31(L,R):
    string = ""
    string += "    int N%s_%s;\n" % (L,R)
    string += "    sp_mat V%s_%s;\n" % (L,R)
    string += "    sp_mat %s_%s();\n" %(L,R)
    return string    

def reset_matrix(L,R):
    string = ""
    string += "    N%s_%s = 0;" % (L,R)
    return string

def update_from31(L,R):
    string = ""
    string += "void flexmat::update_as_%s_%s(sp_mat spC, int Np, int Nq, int Nr, int Ns){\n" % (L,R)
    string += "    iNp = Np;\n"
    string += "    iNq = Nq;\n"
    string += "    iNr = Nr;\n"
    string += "    iNs = Ns;\n"
    string += "    unpack_sp_mat H(spC);\n"
    string += "    v%s = conv_to<uvec>::from(floor(H.vT0/(iN%s*iN%s)));\n" %(L[2], L[0], L[1])
    string += "    v%s = conv_to<uvec>::from(floor((H.vT0 - v%s*iN%s*iN%s)/iN%s));\n" %(L[1], L[2], L[0], L[1], L[0])
    string += "    v%s = conv_to<uvec>::from(H.vT0 - v%s*iN%s*iN%s - v%s*iN%s);\n" %(L[0], L[2], L[0], L[1], L[1], L[0])
    string += "    v%s = conv_to<uvec>::from(H.vT1);\n" % R
    string += "    vValues = H.vVals;\n"
    string += "    deinit();\n"
    return string

def update_from13(L,R):
    string = "\n"
    string += "void flexmat::update_as_%s_%s(sp_mat spC, int Np, int Nq, int Nr, int Ns){\n" % (L,R)
    string += "    iNp = Np;\n"
    string += "    iNq = Nq;\n"
    string += "    iNr = Nr;\n"
    string += "    iNs = Ns;\n"
    string += "    unpack_sp_mat H(spC);\n"
    string += "    v%s = conv_to<uvec>::from(floor(H.vT1/(iN%s*iN%s)));\n" %(R[2], R[0], R[1])
    string += "    v%s = conv_to<uvec>::from(floor((H.vT1 - v%s*iN%s*iN%s)/iN%s));\n" %(R[1], R[2], R[0], R[1], R[0])
    string += "    v%s = conv_to<uvec>::from(H.vT1 - v%s*iN%s*iN%s - v%s*iN%s);\n" %(R[0], R[2], R[0], R[1], R[1], R[0])
    string += "    v%s = conv_to<uvec>::from(H.vT0);\n" % L
    string += "    vValues = H.vVals;\n"
    string += "    deinit();\n"
    string += "}\n"
    return string
 
def header_from31(L,R):
    string = "\n"
    string += "    void update_as_%s_%s(sp_mat spC, int Np, int Nq, int Nr, int Ns);\n" % (L,R)
    string += "    void update_as_%s_%s(sp_mat spC, int Np, int Nq, int Nr, int Ns);\n" % (R,L)
    return string
     
 
string = "\n"
head = "\n"
d = permutations(["p","q","r", "s"])
for i in d:
    string += update_from13(i[0],i[1]+i[2]+i[3])
    string += update_from31(i[0]+i[1]+i[2],i[3])
    head += header_from31(i[0],i[1]+i[2]+i[3])
    

#d = permutations(["p","q","r", "s"])
#for i in d:
#    print reset_matrix(i[0],i[1]+i[2]+i[3]) 
string += head
f = open("flexmat_updates.txt", "w")
f.write(string)
f.close()

"""
string = "\n"
d = permutations(["p","q","r", "s"])
for i in d:
    print reset_matrix(i[0]+i[1], i[2]+i[3])
d = permutations(["p","q","r", "s"])
for i in d:
    print reset_matrix(i[0],i[1]+i[2]+i[3]) 
d = permutations(["p","q","r", "s"]) 
for i in d:
    print reset_matrix(i[0]+i[1]+i[2],i[3])    
"""    

"""
for i in d:
    string += header22(i[0]+i[1], i[2]+i[3])
d = permutations(["p","q","r", "s"])
for i in d:
    string += header13(i[0],i[1]+i[2]+i[3])
d = permutations(["p","q","r", "s"])
for i in d:
    string += header31(i[0]+i[1]+i[2],i[3])
"""    
#f = open("flexmat_perms.txt", "w")
#f.write(string)
#f.close()