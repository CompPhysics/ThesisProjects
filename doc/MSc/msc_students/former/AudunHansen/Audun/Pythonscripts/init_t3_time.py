from numpy import *
from matplotlib.pyplot import *
from scipy.interpolate import *

#estimating time needed to initialize t3 amplitudes

Ns = [3, 4, 5, 6, 7,8,9,10,11]
Ns2 = array([54, 66, 114,162,186,246,294, 342,358])
t = array([0.25, 0.42,2.0, 5.2,8.0,15,24,41,51])

tt = polyfit(Ns2, log(t), 1)
print tt

figure(1)
plot(Ns2, t)
hold("on")
plot(Ns2, exp(tt[1])*exp(tt[0]*Ns2))
legend(["init time", "$%.2f exp(%.2f*N_s)$" % (exp(tt[1]), tt[0])]) #, )
title("Measured t3 amplitude initialization time")
xlabel("Number of states, Ns")
ylabel("Initialization time [s]")
show()

print exp(tt[1])*exp(tt[0]*1000)/3600.0
