#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')
from matplotlib.mlab import load
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)


f = "omg"
t = 491
n = 1024
columns = 23

data = load(f)
tt = data[:,0]

fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twiny()
new_ticks = np.array([i for i in xrange(0,1001,200)])
def tick_function(x):
    return x

#for i in range(1,20):
for i in range(1,14):
    lr = 0.05 * i
    ax.plot(tt, data[:,i+3], label=r'$'+str(lr)+'$')
    ax.grid(color='gray', linestyle='-', linewidth=0.5)
    ax.set_yscale('log')

ax.set_xlabel(r'$T/T_{\rm rlx}(T=0)$')
ax.set_ylabel(r'Lagrange radii')
ax.set_ylim(10**(-1.5), 10)
#ax.legend(loc='upper left',ncol=2)
#plt.show()
ax2.set_xticks(new_ticks)
ax2.set_xticklabels(new_ticks)
ax2.set_xlabel(r'$T$ ($N-$body units)')

#fig.savefig('lagrange.eps', format='eps')
fig.savefig('lagrange.pdf', format='pdf')
