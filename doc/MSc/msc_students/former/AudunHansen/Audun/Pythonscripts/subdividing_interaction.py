from numpy import *
from matplotlib.pyplot import *

figure(1)
hold("on")
c = (.3,.3,.5)
c2 = (1,1,1)
a = 2.0
x = array([0,16])
x1 = array([0,0])
x2 = array([16,16])

plot([-.2,16.2], [-.2,16.2], color = c2) 

plot(x1,x, "-", color = c, linewidth = a)
plot(x2,x, "-", color = c, linewidth = a)

plot(x,x1, "-", color = c, linewidth = a)
plot(x,x2, "-", color = c, linewidth = a)

plot([2,2], x, "--", color = c)
plot([6,6], x, "--", color = c)
plot([10,10], x, "--", color = c)

plot(x,[14,14], "--", color = c)
plot(x,[10,10], "--", color = c)
plot(x,[6,6], "--", color = c)

text(-2,3, "$\langle pp|$")
text(-2,8, "$\langle ph|$")
text(-2,12, "$\langle hp|$")
text(-2,15, "$\langle hh|$")

text(0,17, "$|hh \\rangle$")
text(3,17 ,"$|hp \\rangle$")
text(7,17 ,"$|ph \\rangle$")
text(12,17, "$|pp \\rangle$")

tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off
    
tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels along the bottom edge are off
title('The full interaction matrix with regions')
show()