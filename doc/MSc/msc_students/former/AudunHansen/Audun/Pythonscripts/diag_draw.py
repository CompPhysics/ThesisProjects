from numpy import *
from matplotlib.pyplot import *


def nconnect(n1,n2,S,order="l0", p_h = None):
    N = 60

    Phx = (n1.x+n2.x)/2.0
    Phy = (n1.y+n2.y)/2.0
    
    lP = sqrt((n2.x-n1.x)**2 + (n2.y-n1.y)**2)
    dPx = (n2.x-n1.x)/lP    
    dPy = (n2.y-n1.y)/lP
    
    Cx = Phx - S*dPy
    Cy = Phy + S*dPx  
    
    
    lC = sqrt((S*dPy)**2 + (S*dPx)**2)
    #node(Phx,Phy, c="blue")
    #node(Cx,Cy, c="red")
    R = sqrt((Cx-n1.x)**2 + (Cy-n1.y)**2)

    lPC0 = sqrt(((Cx+R)-n1.x)**2 + (Cy - n1.y)**2)
    lPC1 = sqrt(((Cx+R)-n2.x)**2 + (Cy - n2.y)**2)

    dalpha = arccos((2*R**2 - lP**2)/(2.0*R**2))

    
    CPx = n1.x - Cx
    CPy = n1.y - Cy
    X,Y = 0,0
    if order == "0":
        #X = [n1.x, n2.x]
        #Y = [n1.y, n2.y]
        X = linspace(n1.x, n2.x, 6) #need to do this to get arrows working
        Y = linspace(n1.y, n2.y, 6)
    if order == "l0":
        if S<0:
            dalpha = 2*pi - dalpha
        A = linspace(0,dalpha, N)
        X,Y = rotate_v(CPx,CPy,A)
        X+=Cx
        Y+=Cy
    if order == "r0":
        if S>0:
            dalpha = 2*pi - dalpha      
        A = linspace(0,-dalpha, N)
        X,Y = rotate_v(CPx,CPy,A)
        X+=Cx
        Y+=Cy
        
    msize = 10
    if p_h == 1:
        draw_arrow([X[len(X)/2],Y[len(X)/2]], [-dPx,-dPy])
        #X[len(X)/2],Y[len(X)/2]
        #plot(X[len(X)/2],Y[len(X)/2], "^", color = "black", markersize = msize)
    if p_h == -1:
        draw_arrow([X[len(X)/2],Y[len(X)/2]], [dPx,dPy])
        #plot(X[len(X)/2],Y[len(X)/2], "v", color = "black", markersize = msize)
    plot(X,Y, color = "black")
    
def ncon(n1,n2,order = 0, p_h = None):
    if order == 0:
        nconnect(n1,n2,1,"0", p_h)
    if order > 0:
        nconnect(n1,n2,(-2+order),"l0", p_h)
    if order < 0:
        nconnect(n1,n2,(-2-order),"r0", p_h)

def draw_arrow(pos, point, s = .2, h = .1):
    #Draw an arrow at pos, pointing in the direction of point.
    #normalize direction
    p2 = sqrt(point[0]**2 + point[1]**2)
    point[0] /= p2
    point[1] /= p2
    
    #pi/2 degree rotation
    p_rotx, p_roty = point[1], -point[0]

    x0, y0 = pos[0], pos[1]
    x1, y1 = pos[0] - s*point[0], pos[1] - s*point[1]
    
    #plot the arrow
    plot([x0, x1+h*p_rotx],[y0, y1+h*p_roty], color = "black")
    plot([x0, x1-h*p_rotx],[y0, y1-h*p_roty], color = "black")
                
def rotate_v(x,y,alpha):
    ca = cos(alpha)
    sa = sin(alpha)
    return ca*x - sa*y, sa*x + ca*y            

def nodedraw(p):
    plot(p.x, p.y,"o", color = "black", markersize = 5)

class node():
    def __init__(self, x,y):
        self.x = x
        self.y = y

def node_connect(n1, n2, feature, label = None):
    if feature == 0:
        #amplitude/solid operator
        plot([n1.x, n2.x], [n1.y, n2.y], color = "Black")
        if label:
            text((n2.x-n1.x)/2.0 + label[0],(n2.y-n1.y)/2.0, label[1])
    if feature == 1:
        #interaction vertex
        plot([n1.x, n2.x], [n1.y, n2.y], "--",color = "Black")
    if feature == 2:
        #particle
        plot([n1.x, n2.x], [n1.y, n2.y], color = "Black")
        draw_arrow([n1.x + (n2.x-n1.x)/2.0,n1.y + (n2.y-n1.y)/2.0], [n2.x-n1.x,n2.y-n1.y])
        if label:
            text(n1.x + (n2.x-n1.x)/2.0 + label[0],n1.y + (n2.y-n1.y)/2.0, label[1])
    if feature == 3:
        #hole
        plot([n1.x, n2.x], [n1.y, n2.y], color = "Black")
        draw_arrow([n1.x + (n2.x-n1.x)/2.0,n1.y + (n2.y-n1.y)/2.0], [-(n2.x-n1.x),-(n2.y-n1.y)])
        if label:
            text(n1.x + (n2.x-n1.x)/2.0 + label[0],n1.y +(n2.y-n1.y)/2.0, label[1])
    if feature == 4:
        #particle+hole
        ncon(n1, n2, -1, 1)
        ncon(n2, n1, -1, 1)
        if label:
            text(n1.x + (n2.x-n1.x)/2.0 - label[0],n1.y + (n2.y-n1.y)/2.0, label[1])
            text(n1.x + (n2.x-n1.x)/2.0 + label[0],n1.y + (n2.y-n1.y)/2.0, label[2])
    if feature == 5:
        #connect all nodes in list n1 using solid vertex
        for i in range(len(n1)-1):
            node_connect(n1[i], n1[i+1], 0)
            nodedraw(n1[i])
        nodedraw(n1[-1])
    if feature == 6:
        #connect all nodes in list n1 using dotted vertex
        for i in range(len(n1)-1):
            node_connect(n1[i], n1[i+1], 1)
            nodedraw(n1[i])
        nodedraw(n1[-1])
        
def Cr(n):
    return "\\hat{a}_%s^\\dagger" %n
    
def An(n):
    return "\\hat{a}_%s" %n

def simple_lines():
    spread = .2
    diagname = "$\\vert \\Phi^{a}_{i} \\rangle$"
    figure(figsize = (4,3), dpi = 80, edgecolor = 'k',facecolor='white')
    plot([-.1,1.1], [-.1,1.1], color = "white")
    axis('off')
    title(diagname)
    axes().set_aspect('equal', 'datalim')
    hold("on")
    
    #interaction = [node(0,1), node(1,1)]
    above = [node(-spread, 2), node(spread, 2), node(1-spread, 2), node(1+spread, 2)]
    below = [node(-spread, 0), node(spread, 0), node(1-spread, 0), node(1+spread, 0)]
    node_connect(below[0], above[0], 2, [-.2, "a"])
    #node_connect(below[3], above[3], 3, [.2, "j"])
    node_connect(below[3], above[3], 3, [-.2, "i"])
    #node_connect(below[2], above[2], 2, [.2, "b"])
    show()

def contraction():
    spread = .4
    #diagname = "$\\sum_{abij} \\langle ai \\vert \\vert bj \\rangle \\hat{a}_i^\\dagger$"
    #"\hat{F}_N = "
    #diagname ="$\\sum_{ij} f_{ij}  \\hat{a}_j \\hat{a}_i^{\\dagger}$" 
    #diagname ="$\\sum_{ab} f_{ab}  \\hat{a}_a^{\\dagger} \\hat{a}_b$" 
    #diagname ="$\\sum_{ai} f_{ia}  \\hat{a}_i^{\\dagger} \\hat{a}_a$"
    #diagname ="$\\sum_{ai} f_{ai}  \\hat{a}_i \\hat{a}_a^{\\dagger}$"
    diagname ="$\\sum_{bc} f_{bc}  \\hat{a}_b^{\\dagger} \\hat{a}_c$" 
    
    figure(figsize = (3,4), dpi = 80, edgecolor = 'k',facecolor='white')
    plot([-.1,2.2], [-.1,2.2], color = "white")
    axis('off')
    title("Example of contraction")
    axes().set_aspect('equal', 'datalim')
    hold("on")
    
    interaction = [node(.5,1), node(1,1)]
    above = [node(-spread, 2), node(spread, 2), node(1-spread, 2), node(1+spread, 2)]
    below = [node(-spread, 0), node(spread, 0), node(1-spread, 0), node(1+spread, 0)]
    
    SDbl = [node(-spread, -2), node(spread, -2), node(1-spread, -2), node(1+spread, -2)]
    SDab = [node(-spread, -.5), node(spread, -.5), node(1-spread, -.5), node(1+spread, -.5)]
    
    #node_connect(interaction, above, 6)
    node_connect(interaction[0], interaction[1], 1)
    nodedraw(interaction[1])
    
    #node_connect(below[0], above[0], 3)
    node_connect(interaction[1], above[-1], 2, [-.2, "b"])
    #node_connect(interaction[1], above[-2], 3, [-.2, "i"])
    
    #node_connect(interaction[1], below[-1], 2, [-.2, "i"])
    node_connect(interaction[1], below[-1], 3, [-.2, "c"])
    #node_connect(interaction[1], below[-1], 3)
    
    node_connect(SDbl[1], SDab[1], 3, [-.2, "i"])
    node_connect(SDbl[3], SDab[3], 2, [-.2, "a"])
    #node_connect(interaction[1], below[-1], 3)
    
    text(2, 0, "=") 
    h = -.5
    v = .2
    interaction = [node(2.6 +v,.5+ h), node(3 + v,.5 +h)]
    above = [node(2 + v -spread, 1.5+h), node(2 +spread + v, 1.5+h), node(3-spread + v, 1.5+h), node(3+spread +v, 1.5+h)]
    below = [node(2 + v -spread, -.5+h), node(2 +spread + v, -.5+h), node(3-spread + v, -.5+h), node(3+spread +v, -.5+h)]
    
    node_connect(interaction[0], interaction[1], 1)
    nodedraw(interaction[1])
    
    #node_connect(below[0], above[0], 3)
    node_connect(interaction[1], above[-1], 2, [-.2, "b"])
    #node_connect(interaction[1], above[-2], 3, [-.2, "i"])
    
    #node_connect(interaction[1], below[-1], 2, [-.2, "i"])
    node_connect(interaction[1], below[-1], 3, [-.2, "a"])
    node_connect(below[1], above[1], 2, [-.2, "i"])

    #node_connect(n1,n2,2, [-.2, "a"])
    #ncon(n1, n2, -1, 1)
    #ncon(n2, n1, -1, 1)
    #nodedraw(n1)
    show()

def onebody_op():
    spread = .4
    #diagname = "$\\sum_{abij} \\langle ai \\vert \\vert bj \\rangle \\hat{a}_i^\\dagger$"
    #"\hat{F}_N = "
    #diagname ="$\\sum_{ij} f_{ij}  \\hat{a}_j \\hat{a}_i^{\\dagger}$" 
    diagname ="$\\sum_{bc} f_{bc}  \\hat{a}_b^{\\dagger} \\hat{a}_c$" 
    #diagname ="$\\sum_{ai} f_{ia}  \\hat{a}_i^{\\dagger} \\hat{a}_b$"
    #diagname ="$\\sum_{ai} f_{ai}  \\hat{a}_i \\hat{a}_a^{\\dagger}$"
    
    figure(figsize = (3,4), dpi = 80, edgecolor = 'k',facecolor='white')
    plot([-.1,2.2], [-.1,2.2], color = "white")
    axis('off')
    title(diagname)
    axes().set_aspect('equal', 'datalim')
    hold("on")
    
    interaction = [node(.5,1), node(1,1)]
    above = [node(-spread, 2), node(spread, 2), node(1-spread, 2), node(1+spread, 2)]
    below = [node(-spread, 0), node(spread, 0), node(1-spread, 0), node(1+spread, 0)]
    #node_connect(interaction, above, 6)
    node_connect(interaction[0], interaction[1], 1)
    nodedraw(interaction[1])
    
    #node_connect(below[0], above[0], 3)
    #node_connect(interaction[1], above[-1], 2, [-.2, "a"])
    #node_connect(interaction[1], above[-2], 3, [-.2, "i"])
    
    node_connect(interaction[1], below[-1], 2, [-.2, "i"])
    node_connect(interaction[1], below[-2], 3, [-.2, "a"])
    #node_connect(interaction[1], below[-1], 3)
    
    

    #node_connect(n1,n2,2, [-.2, "a"])
    #ncon(n1, n2, -1, 1)
    #ncon(n2, n1, -1, 1)
    #nodedraw(n1)
    show()

def normal_hamilt():
    #draw_arrow([0,0], [1,1])
    spread = .3
    #diagname = "$\\sum_{abij} \\langle ai \\vert \\vert bj \\rangle \\hat{a}_i^\\dagger$"
    
    a = "a"
    b = "b"
    c = "c"
    d = "d"
    i = "i"
    j = "j"
    k = "k"
    l = "l"

    #diagname = "$\\sum_{abcd} \\langle ab\\vert \\hat{g} \\vert cd \\rangle \\{ %s %s %s %s \\}$" % (Cr(a), Cr(b), An(d), An(c))
    #diagname ="$\\sum_{ijkl} \\langle ij\\vert \\hat{g} \\vert kl \\rangle \\{ %s %s %s %s \\}$" % (Cr(i), Cr(j), An(l), An(k))
    #diagname ="$\\sum_{aibj} \\langle ij\\vert \\hat{g} \\vert bj \\rangle \\{ %s %s %s %s \\}$" % (Cr(a), Cr(i), An(j), An(b))
    #diagname ="$\\sum_{abci} \\langle ab\\vert \\hat{g} \\vert ci \\rangle \\{ %s %s %s %s \\}$" % (Cr(a), Cr(b), An(i), An(c))
    #diagname ="$\\sum_{iajk} \\langle ia\\vert \\hat{g} \\vert jk \\rangle \\{ %s %s %s %s \\}$" % (Cr(i), Cr(a), An(k), An(j))
    #diagname ="$\\sum_{aibc} \\langle ai\\vert \\hat{g} \\vert bc \\rangle \\{ %s %s %s %s \\}$" % (Cr(a), Cr(i), An(b), An(c))
    #diagname ="$\\sum_{ijka} \\langle ij\\vert \\hat{g} \\vert ka \\rangle \\{ %s %s %s %s \\}$" % (Cr(i), Cr(j), An(a), An(k))
    #diagname ="$\\sum_{abij} \\langle ab\\vert \\hat{g} \\vert ij \\rangle \\{ %s %s %s %s \\}$" % (Cr(a), Cr(b), An(j), An(i)) 
    diagname ="$\\sum_{ijab} \\langle ij\\vert \\hat{g} \\vert ab \\rangle \\{ %s %s %s %s \\}$" % (Cr(i), Cr(j), An(b), An(a))
    
    
    figure(figsize = (1.0,3.1), dpi = 80, edgecolor = 'k',facecolor='white')
    plot([-.4,1.4], [-.3,1.4], color = "white")
    axis('off')
    #title(diagname)
    text(-.5,1.4,diagname, size = 13)
    text(-.3,-.3,"+", size = 20)
    text(.2,-.3,"-", size = 20)
    text(.7,-.3,"+", size = 20)
    text(1.2,-.3,"-", size = 20)
    axes().set_aspect('equal', 'datalim')
    hold("on")
    
    interaction = [node(0,1), node(1,1)]
    above = [node(-spread, 2), node(spread, 2), node(1-spread, 2), node(1+spread, 2)]
    below = [node(-spread, 0), node(spread, 0), node(1-spread, 0), node(1+spread, 0)]
    node_connect(interaction, above, 6)    
    node_connect(interaction[0], below[0], 3,[-.2, ""])  #a out
    node_connect(interaction[0], below[1], 2,[-.2, ""])  #i out
    node_connect(interaction[1], below[-2], 3,[-.2, ""]) #b out
    node_connect(interaction[1], below[-1], 2,[-.2, ""]) #j out
    
    #node_connect(interaction[0], below[0], 2,[-.2, "i"])  #a out
    #node_connect(interaction[1], below[-2], 2,[-.2, "j"])  #i out
    #node_connect(interaction[1], below[-1], 3,[-.2, "a"]) #b out
    #node_connect(interaction[1], below[-1], 2,[-.2, "j"]) #j out
    
    #node_connect(interaction[0], below[0], 2) 
    #node_connect(interaction[1], below[-1], 3)
    
    

    #node_connect(n1,n2,2, [-.2, "a"])
    #ncon(n1, n2, -1, 1)
    #ncon(n2, n1, -1, 1)
    #nodedraw(n1)
    show()
    
def normal_hamilt2():
    #draw_arrow([0,0], [1,1])
    spread = .3
    #diagname = "$\\sum_{abij} \\langle ai \\vert \\vert bj \\rangle \\hat{a}_i^\\dagger$"
    
    a = "a"
    b = "b"
    c = "c"
    d = "d"
    i = "i"
    j = "j"
    k = "k"
    l = "l"

    #diagname = "$\\sum_{abcd} \\langle ab\\vert \\hat{g} \\vert cd \\rangle \\{ %s %s %s %s \\}$" % (Cr(a), Cr(b), An(d), An(c))
    #diagname ="$\\sum_{ijkl} \\langle ij\\vert \\hat{g} \\vert kl \\rangle \\{ %s %s %s %s \\}$" % (Cr(i), Cr(j), An(l), An(k))
    #diagname ="$\\sum_{aibj} \\langle ij\\vert \\hat{g} \\vert bj \\rangle \\{ %s %s %s %s \\}$" % (Cr(a), Cr(i), An(j), An(b))
    #diagname ="$\\sum_{abci} \\langle ab\\vert \\hat{g} \\vert ci \\rangle \\{ %s %s %s %s \\}$" % (Cr(a), Cr(b), An(i), An(c))
    #diagname ="$\\sum_{iajk} \\langle ia\\vert \\hat{g} \\vert jk \\rangle \\{ %s %s %s %s \\}$" % (Cr(i), Cr(a), An(k), An(j))
    #diagname ="$\\sum_{aibc} \\langle ai\\vert \\hat{g} \\vert bc \\rangle \\{ %s %s %s %s \\}$" % (Cr(a), Cr(i), An(b), An(c))
    #diagname ="$\\sum_{ijka} \\langle ij\\vert \\hat{g} \\vert ka \\rangle \\{ %s %s %s %s \\}$" % (Cr(i), Cr(j), An(a), An(k))
    #diagname ="$\\sum_{abij} \\langle ab\\vert \\hat{g} \\vert ij \\rangle \\{ %s %s %s %s \\}$" % (Cr(a), Cr(b), An(j), An(i)) 
    diagname ="$\\sum_{ijab} \\langle ij\\vert \\hat{g} \\vert ab \\rangle \\{ %s %s %s %s \\}$" % (Cr(i), Cr(j), An(b), An(a))
    
    
    figure(figsize = (1.0,3.1), dpi = 80, edgecolor = 'k',facecolor='white')
    #plot([-.5,1.4], [-.0,1.4], color = "white")
    axis('off')
    #title(diagname)
    #text(-.5,1.4,diagname, size = 13)
    #text(-.3,-.3,"+", size = 20)
    #text(.2,-.3,"-", size = 20)
    #text(.7,-.3,"+", size = 20)
    #text(1.2,-.3,"-", size = 20)
    axes().set_aspect('equal', 'datalim')
    hold("on")
    
    interaction = [node(0,1), node(2,1)]
    above = [node(-spread, 2), node(spread, 2), node(1-spread, 2), node(1+spread, 2)]
    below = [node(-spread, -1), node(spread, -1), node(2-spread, -1), node(2+spread, -1)]
    node_connect(interaction, above, 6)    
    node_connect(interaction[0], below[0], 3,[-.2, ""])  #a out
    node_connect(interaction[0], below[1], 2,[-.2, ""])  #i out
    #node_connect(interaction[1], below[-2], 3,[-.2, ""]) #b out
    #node_connect(interaction[1], below[-1], 2,[-.2, ""]) #j out
    
    #node_connect(interaction[0], below[0], 2,[-.2, "i"])  #a out
    #node_connect(interaction[1], below[-2], 2,[-.2, "j"])  #i out
    #node_connect(interaction[1], below[-1], 3,[-.2, "a"]) #b out
    #node_connect(interaction[1], below[-1], 2,[-.2, "j"]) #j out
    
    #node_connect(interaction[0], below[0], 2) 
    #node_connect(interaction[1], below[-1], 3)
    
    

    #node_connect(n1,n2,2, [-.2, "a"])
    #ncon(n1, n2, -1, 1)
    #ncon(n2, n1, -1, 1)
    #nodedraw(n1)
    show()
 
def sign(node, s, sc):
    if s==0:
        text(node.x-.09-.2, node.y+sc, "-", size = 17)
    else: 
        text(node.x+.09-.2, node.y+sc, "+", size = 17)

def clusters():
    spread = .4
    diagname = "$\hat{T}_1$"
    figure(figsize = (4,1.5), dpi = 80, edgecolor = 'k',facecolor='white')
    plot([-1,1], [-1,2.4], color = "white")
    axis('off')
    #title(diagname)
    axes().set_aspect('equal', 'datalim')
    hold("on")
    
    #interaction = [node(0,1), node(1,1)]
    #above = [node(-spread, 2), node(spread, 2), node(1-spread, 2), node(1+spread, 2)]
    #below = [node(0, 0), node(spread, 0), node(1-spread, 0), node(1+spread, 0)]
    
    above = [node(-spread, 2), node(spread, 2), node(1-spread, 2), node(1+spread, 2), node(2-spread, 2), node(2+spread, 2), node(3-spread, 2), node(3+spread, 2)]
    below = [node(0, 0), node(1, 0), node(2, 0), node(3, 0)]
    
    text(1,-1, "$\hat{T}_4$")
    
    node_connect(below[0], above[0], 2, [-.4, ""])
    node_connect(below[0], above[1], 3, [.3, ""])
    
    sign(above[0], 1, .2)
    sign(above[1], 0, .2)
    
    node_connect(below[1], above[2], 2, [-.4, ""])
    node_connect(below[1], above[3], 3, [.3, ""])
    sign(above[2], 1, .2)
    sign(above[3], 0, .2)
    
    node_connect(below[2], above[4], 2, [-.4, ""])
    node_connect(below[2], above[5], 3, [.3, ""])
    sign(above[4], 1, .2)
    sign(above[5], 0, .2)
    
    node_connect(below[3], above[6], 2, [-.4, ""])
    node_connect(below[3], above[7], 3, [.3, ""])
    sign(above[6], 1, .2)
    sign(above[7], 0, .2)
    
    #node_connect(below[3], above[3], 3, [.2, "j"])
    

    plot([-.1,3.1], [0,0], color = (0,0,0), linewidth = 2)
    show()

normal_hamilt2()
#simple_lines()
#contraction()
#clusters()