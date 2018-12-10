#< ab  || cd  > = [[ a,b ] , [ c,d ]]
#A script to find the optimal alignment of diagrams used in the CCDT t3 amplitude equation



def perm(a, i,e):

    ai= a[1][e]
    ae = a[1][i]
    api = a[3][e]
    ape = a[3][i]
    a[1][i] = ai
    a[1][e] = ae
    a[3][i] = api
    a[3][e] = ape

def perm2(a, i,e):

    ai= a[0][e]
    ae = a[0][i]
    api = a[2][e]
    ape = a[2][i]
    a[0][i] = ai
    a[0][e] = ae
    a[2][i] = api
    a[2][e] = ape

def align(a, b, c, left_indices):
    #1. assign all left_indices in a[0], a[2]
    a_ = [[],[],[],[]]
    b_ = [[],[],[],[]]
    c_ = [[],[],[],[]]
    for i in range(len(a[0])):
        #bra
        if a[0][i] in ["d", "e", "f"]:
            #move to ket 
            a_[1].append(a[0][i])
            a_[3].append(a[2][i])
        if a[0][i] in ["a", "b", "c"]:
            #keep in bra
            a_[0].append(a[0][i])
            a_[2].append(a[2][i])
        
        #ket
        if a[1][i] in ["i", "j", "k"]:
            #move to bra
            a_[0].append(a[1][i])
            a_[2].append(a[3][i])
        if a[1][i] in ["l", "m", "n"]:
            #keep in ket
            a_[1].append(a[1][i])
            a_[3].append(a[3][i])    
        
    
    #2. assign all left indices in b to a_[1]
    for i in range(len(b[0])):
        if b[0][i] in a_[1]:
            b_[0].append(b[0][i])
            b_[2].append(b[2][i])
        if b[0][i] not in a_[1]:
            b_[1].append(b[0][i])
            b_[3].append(b[2][i])  
    for i in range(len(b[0])):       
        if b[1][i] in a_[1]:
            b_[0].append(b[1][i])
            b_[2].append(b[3][i])        
        if b[1][i] not in a_[1]:
            b_[1].append(b[1][i])
            b_[3].append(b[3][i])
    #ensure correct order in a[1]
    #a_temp = a_
    print b_
    print a_
    for i in range(len(a_[1])):
        if a_[1][i] != b_[0][i]:
            for e in range(len(a_[1])):
                if a_[1][e] == b_[0][i]:
                    perm(a_, e,i)

    
    #3. align c to b_[1]
    for i in range(len(c[0])):
        if c[0][i] in b_[1]:
            c_[0].append(c[0][i])
            c_[2].append(c[2][i])
        if c[0][i] not in b_[1]:
            c_[1].append(c[0][i])
            c_[3].append(c[2][i])  
    for i in range(len(c[0])):       
        if c[1][i] in b_[1]:
            c_[0].append(c[1][i])
            c_[2].append(c[3][i])        
        if c[1][i] not in b_[1]:
            c_[1].append(c[1][i])
            c_[3].append(c[3][i])

    for i in range(len(c_[0])):
        if b_[1][i] != c_[0][i]:

            for e in range(len(c_[0])):
                if c_[0][e] == b_[1][i]:

                    perm2(c_, i,e)   
    #print "A:", a_
    #print "B:", b_
    #print "C:", c_
    return a_,b_,c_

 
                

    
def diagsort(a,c):
    #align diagram to the T3 amplitude
    nr = {"a": "p", "b": "q","c": "r", "i": "s","j": "t", "k": "u" }
    
    retrs = "update_as_"
    for i in range(len(a[0])):
        retrs += nr[a[0][i]]
    retrs += "_"
    for i in range(len(c[1])):
        retrs += nr[c[1][i]]
    return retrs
    

            
    
    #align to t3 amp
        

def setup(a,b,c):
    #assign general indices pqrs
    a = [a[0], a[1], [],[]]
    b = [b[0], b[1], [],[]]
    c = [c[0], c[1], [],[]]
    indx = "pqrstu"
    n = 0
    for i in range(len(a[0])):
        a[2].append(indx[n])
        n+= 1
    for i in range(len(a[1])):
        a[3].append(indx[n])
        n+= 1
    
    n = 0
    for i in range(len(b[0])):
        b[2].append(indx[n])
        n+= 1
    for i in range(len(b[1])):
        b[3].append(indx[n])
        n+= 1
    
    n = 0
    for i in range(len(c[0])):
        c[2].append(indx[n])
        n+= 1
    for i in range(len(c[1])):
        c[3].append(indx[n])
        n+= 1
        
    #identify left indices
    left_indices = []
    for i in range(len(a[0])):
        if a[0][i] in ["a", "b", "c"]:
            left_indices.append(a[0][i])
        if a[1][i] in ["i", "j", "k"]:
            left_indices.append(a[1][i])           
    
    a,b,c = align(a,b,c, left_indices)
    
    """
      
    #align indices in a,b  
    diag = [[],[]]
    ap = [[],[]]
    bp = [[],[]]
    cp = [[],[]]
    #1. identify open lines in a

    for i in range(len(a)):
        if a[0][i] in ["d", "e", "f"]:
            diag[0].append(a[0][i])
            ap[0].append(a[0][i])
            #a_s.append(A[0][i])
        if a[1][i] in ["i", "j", "k"]:
            diag[0].append(a[1][i])
            ap[0].append(a[1][i])
            #a_s.append(A[1][i])
        if a[0][i] not in ["d", "e", "f"]:
            ap[1].append(a[0][i])
        if a[1][i] not in ["l", "m", "n"]:
            ap[1].append(a[1][i])
    
    #align closed lines in a-b
    
    for i in range(len(ap[1])):
        pass
        
    

    a_s = "."
    b_s = "."
    c_s = "."  
    """     
    #2. use internal lines from a to form first part of b
    return a,b,c
        

def generate_t2t2(v,t2,t3):
    #measure "level of alignment" of existing tensors
    #we ideally want it to begin with abc, and end with ijk
    #contractions occur over lmn and def
    t3ind = 0
    
    contractions = ["l","m","d","e"]
    
    
    #begin by evaluate where to place the t3 amplitudes
    for i in range(len(t3[0])):
        if t3[0][i] in ["a", "b"]:
            t3ind += 1
        if t3[1][i] in ["i", "j"]:
            t3ind -= 1
            
    #inspect if t2 has a preferred placement
    for i in range(len(t2[0])):
        if t2[0][i] in ["a", "b"]:
            t3ind += 1
        if t2[1][i] in ["i", "j"]:
            t3ind -= 1   
    #print t3ind
    
    if t3ind >= 0:
        #place t3 first
        a,b,c = setup(t3, v, t2)
        #a = t3
        t3str = "t3."
        for i in range(len(a[2])):
            t3str += a[2][i]
        t3str += "_"
        for i in range(len(a[3])):
            t3str += a[3][i]
        t3str += "()"
        
        t2str = "t2."
        for i in range(len(c[2])):
            t2str += c[2][i]
        t2str += "_"
        for i in range(len(c[3])):
            t2str += c[3][i]
        t2str += "()"
        
        vint = "vhhpp."
        for i in range(len(b[2])):
            vint += b[2][i]
        vint += "_"
        for i in range(len(b[3])):
            vint += b[3][i]
        vint += "()"
        
        matmult = t3str + "*" + vint + "*" + t2str
        
        
    else:
        #place t3 last  
        a,b,c = setup(t2, v, t3)
        
        t2str = "t3."
        for i in range(len(a[2])):
            t2str += a[2][i]
        t2str += "_"
        for i in range(len(a[3])):
            t2str += a[3][i]
        t2str += "()"
        
        t3str = "t2."
        for i in range(len(c[2])):
            t3str += c[2][i]
        t3str += "_"
        for i in range(len(c[3])):
            t3str += c[3][i]
        t3str += "()"
        
        vint = "vhhpp."
        for i in range(len(b[2])):
            vint += b[2][i]
        vint += "_"
        for i in range(len(b[3])):
            vint += b[3][i]
        vint += "()"
        
        matmult = t2str + "*" + vint + "*" + t3str
    #print matmult
    retstr = diagsort(a,c)
    strng = retstr + "("  + matmult + ")"
    #print a
    #print b        
    #print c      
    return a, b, c, strng        
    

def generate(v,t2,t3):
    #measure "level of alignment" of existing tensors
    #we ideally want it to begin with abc, and end with ijk
    #contractions occur over lmn and def
    t3ind = 0
    
    contractions = ["l","m","d","e"]
    
    
    #begin by evaluate where to place the t3 amplitudes
    for i in range(len(t3[0])):
        if t3[0][i] in ["a", "b", "c"]:
            t3ind += 1
        if t3[1][i] in ["i", "j", "k"]:
            t3ind -= 1
            
    #inspect if t2 has a preferred placement
    for i in range(len(t2[0])):
        if t2[0][i] in ["a", "b", "c"]:
            t3ind += 1
        if t2[1][i] in ["i", "j", "k"]:
            t3ind -= 1   
    #print t3ind
    
    if t3ind >= 0:
        #place t3 first
        a,b,c = setup(t3, v, t2)
        #a = t3
        t3str = "t3."
        for i in range(len(a[2])):
            t3str += a[2][i]
        t3str += "_"
        for i in range(len(a[3])):
            t3str += a[3][i]
        t3str += "()"
        
        t2str = "t2."
        for i in range(len(c[2])):
            t2str += c[2][i]
        t2str += "_"
        for i in range(len(c[3])):
            t2str += c[3][i]
        t2str += "()"
        
        vint = "vhhpp."
        for i in range(len(b[2])):
            vint += b[2][i]
        vint += "_"
        for i in range(len(b[3])):
            vint += b[3][i]
        vint += "()"
        
        matmult = t3str + "*" + vint + "*" + t2str
        
        
    else:
        #place t3 last  
        a,b,c = setup(t2, v, t3)
        
        t2str = "t3."
        for i in range(len(a[2])):
            t2str += a[2][i]
        t2str += "_"
        for i in range(len(a[3])):
            t2str += a[3][i]
        t2str += "()"
        
        t3str = "t2."
        for i in range(len(c[2])):
            t3str += c[2][i]
        t3str += "_"
        for i in range(len(c[3])):
            t23tr += c[3][i]
        t3str += "()"
        
        vint = "vhhpp."
        for i in range(len(b[2])):
            vint += b[2][i]
        vint += "_"
        for i in range(len(b[3])):
            vint += b[3][i]
        vint += "()"
        
        matmult = t2str + "*" + vint + "*" + t3str
    #print matmult
    retstr = diagsort(a,c)
    strng = retstr + "("  + matmult + ")"
    #print a
    #print b        
    #print c      
    return a, b, c, strng



def tex_pre(v,t2,t3):
    tx = " \\sum_{"
    for i in range(len(v[0])):
        tx += v[0][i] + v[1][i]
    tx += "} "
    tx += "\\langle %s %s \\vert \\vert %s %s \\rangle " % (v[0][0], v[0][1], v[1][0], v[1][1])
    tx += "t^{%s %s}_{%s %s}" % (t2[0][0], t2[0][1], t2[1][0], t2[1][1])
    #tx += "t^{%s %s %s}_{%s %s %s} " % (t3[0][0], t3[0][1], t3[0][2], t3[1][0], t3[1][1],t3[1][2])
    tx += "t^{%s %s}_{%s %s} " % (t3[0][0], t3[0][1], t3[1][0], t3[1][1])
    return tx

def tex_aligned(a,b,c):
    tx = " \\sum_{"
    for i in b[0]:
        tx+=i
    tx += "}"
    
    tx += " \\sum_{"
    for i in b[1]:
        tx+=i
    tx += "}"
    
    
    
    tx += " t^{"
    for i in a[0]:
        tx += i
    tx += "}_{"
    
    for i in a[1]:
        tx += i
    tx += "} \\langle "
    for i in b[0]:
        tx += i
    tx += "\\vert \\vert "
    for i in b[1]:
        tx += i
    tx +="\\rangle t^{"
    for i in c[0]:
        tx += i
    tx += "}_{"
    
    for i in c[1]:
        tx += i
    tx += "} "
    return tx

def gen_entry(v,t2,t3):
    #Generate table entry for diagram given by t2,t3,v
    tx1 =  tex_pre(v,t2,t3)
    a,b,c, strng = generate_t2t2(v,t2,t3)
    tx2 = tex_aligned(a,b,c)
    return "$$ " + tx1 + " \\rightarrow " + tx2 + "$$" , strng

v = [["l"],["d"]]
t2 = [["a","d"],["i","j"]]
t3 = [["b","c"],["l","k"]]

ltx, strng = gen_entry(v,t2,t3)

print ltx
print strng


def realign_diagram(d1,d2):
    n = {a:p, b:q, c:r, i:s, j:t, k:u}
    
