class qnode():
    #Object for keeping track of vertices
    def __init__(self, qa):
        #self.pos = pos
        self.connected = 0
        self.qa = qa #creation or annihilation +1/-1
        self.group = []
    def connect(i):
        self.connected = i

class qop():
    #a collection of vertices -- in principle an operator
    def __init__(self,qnodes):
        self.qnodes = qnodes
        self.pos = [0,0]
        self.setpos(self.pos)
    def setpos(self, pos):
        self.pos = pos
        for i in range(len(self.qnodes)):
            self.qnodes[i].qpos = self.pos[0] + i
            

class qvert():
    #A vertex in the operator with input, output
    def __init__(self, qin, qout):
        self.qin = qin   #qnode in
        self.qout = qout #qnode out
        self.qpos = 0.0


def CreateOperatorMap(o_cop):
    O = o_cop
    excluded = []
    olist = []
    for i in range(len(O)):
        if O[i] == +1:
            #find -1
            for e in range(len(O)):
                if e!=i:
                    if O[e] == -1:
                        if e not in excluded:
                            if i not in excluded:
                                excluded.append(i)
                                excluded.append(e)
                                olist.append(qvert(qnode(1), qnode(-1)))
    return qop(olist)
                              

print CreateOperatorMap([1,1,1,-1,-1,-1]).qnodes
