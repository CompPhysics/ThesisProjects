from sympy import *



def arrange(l):
    c = l.strip(" ")
    c = c.split("+")
    
    print c
    newlist = []
    for i in c:
        inew = i.strip(" ")
        newlist_c = []
        for e in range(len(inew)):
            
            en = inew[e]
            if en in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "*"]:
                pass
            else:
                if inew[e+1:e+3] == "**":
                    for i in range(int(inew[e+3])):
                        newlist_c.append(en)
                    #newlist.append(newlist_c)
                else:
                    newlist_c.append(en)
        if len(newlist_c) != 0:
            newlist.append(newlist_c)   
    #print newlist
    #print replace_list(["x","y","z"], [1,2,3], newlist)
    #print newlist

    return newlist

def replace_list(s, r, L):
    for i in range(len(s)):
        for e in range(len(L)):
            for u in range(len(L[e])):
                if L[e][u] == s[i]:
                    L[e][u] = r[i]
    return L

def listexpand(n, T):
    To  = [Symbol("x"), Symbol("y"), Symbol("z"), Symbol("r")]
    
    Tn = To[0]
    for i in range(len(T)-1):
        Tn += To[i+1]
    print Tn
    F = Symbol("1")
    for i in range(n+1):
        F += Tn**i
    
        
        
    expanded = str(expand(F))
    #expanded = F
    print expanded
    a = arrange(expanded)
    #print a
    Tns = str(Tn).split("+")
    for i in range(len(Tns)):
        Tns[i] = Tns[i].strip(" ")
    print "Tns:", Tns
    b = replace_list(Tns, T, a)
    #print b
    return b
    

print listexpand(3, ["1","2","3"])

#arrange(str(expand(f**10)))



