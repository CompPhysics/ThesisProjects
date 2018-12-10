from numpy import *
from matplotlib.pyplot import *
#unique configs

Np = linspace(0,2000, 2001)
U = Np*(Np+1)*(Np+2)/6
Uf = Np**3
plot(Np, U)
hold("on")
plot(Np,Uf)
show()

#possibility: the number of configs corresponding to triple amplitude makes this easy to perform. ?