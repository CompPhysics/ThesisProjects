from numpy import *


def permute(order, o1, o2):
    #permutes two element
    p1 = order[o1]
    p2 = order[o2]
    order2  = order
    order2[o1] = p2
    order2[o2] = p1
    
    return order2

def seq_permute(order,o1,o2):
    #sequential permutations
    order2 = order
    for i in range(len(o1)):
        order2 = permute(order2, o1[i], o2[i])
    return order2

def perm_gen(label, order):
    tx = ""
    tx = label + "." + order[0]+order[1]+order[2]+"_"+order[3]+order[4]+order[5]+"()"
    return tx
    
def tri_perm(label, order, o1, o2, o3):
    tx = perm_gen(label, order)
    tx += "-" + perm_gen(label, permute(order, o1,o2)) + "-" + perm_gen(label, permute(order, o1,o3))
    return tx

def tri_split(label, order, o1, o2, o3, o4, o5, o6):
    tx = perm_gen(label, order)
    tx += "-" + perm_gen(label, permute(order, o4,o5))
    tx += "-" + perm_gen(label, permute(order, o4,o6))
    tx += "-" + perm_gen(label, permute(order, o1,o2))
    
    tx += "+" + perm_gen(label, seq_permute(order, [o1,o4],[o2,o5]))
    tx += "+" + perm_gen(label, seq_permute(order, [o1,o4],[o2,o6]))
    
    tx += "-" + perm_gen(label, permute(order, o1,o3))
    
    tx += "+" + perm_gen(label, seq_permute(order, [o1,o4],[o3,o5]))
    tx += "+" + perm_gen(label, seq_permute(order, [o1,o4],[o3,o6]))
    
    return tx

order = ["p", "q", "r", "s", "t", "u"]
label = "t2b"

#print seq_permute(order, [0,3], [1,4])

#print perm_gen(label, permute(order, 0,1))
#print tri_perm(label, order, 5,4,3)

print tri_split(label, order, 3,4,5,2,0,1)


