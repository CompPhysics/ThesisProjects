import numpy as np
from  matplotlib import pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.units as units
import matplotlib.ticker as ticker
rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Double-well potential']})
font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
m = 0.067
homega = 3.0
Lx = 50.0
x = np.linspace(-120.0, 120.0)
v = 0.5*m*homega*(x*x-2*Lx*abs(x)+Lx*Lx)

plt.plot(x, v, 'b-')
plt.title(r'{\bf Double-well potential for $L_x=50$ [nm]}', fontsize=20)     
plt.text(-100, 350, r'Parameters: $m*=0.067m_e$, $\hbar\omega=3.0$ [meV]', fontdict=font)
plt.text(-100, 300, r'$L_x=50$ [nm]', fontdict=font)
plt.xlabel(r'$x$ [nm]',fontsize=20)
plt.ylabel(r'$V(x,0)$ [MeV]',fontsize=20)

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.savefig('double.pdf', format='pdf')
