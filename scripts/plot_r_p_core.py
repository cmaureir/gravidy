#!/usr/bin/env python

from matplotlib.mlab import load
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np

rc('text', usetex=True)

tntime = []
tttime = []
radius = []
density = []
f = "r_p_2.out"
f = open(f)

for line in f.readlines():
    if line == None: break
    line = line.split()
    line = filter(None,line)

    tn = float(line[0])
    tt = float(line[1])
    r  = float(line[2])
    p  = float(line[3])

    tntime.append(tn)
    tttime.append(tt)
    radius.append(r)
    density.append(p)

radius  = np.array(radius)/radius[0]
density = np.array(density)/density[0]

fig = plt.figure(1)

ax1 = fig.add_subplot(212)
ax2 = fig.add_subplot(211)
ax3 = ax2.twiny()


ax1.plot(tttime, radius,  label=r'$r_{c}$', color='red')
ax2.plot(tttime, density, label=r'$\rho_{c}$', color='red')

tnbu = [100,200,300,400]
ax3.set_xticks(tnbu)
ax3.set_xticklabels(tnbu)
ax2.xaxis.set_visible(False)

ax3.set_xlabel(r"$T$($N-$body units)", fontsize=20)
ax1.set_xlabel(r"$T/T_{\rm rlx}(T=0)$", fontsize=20)

ax1.set_ylabel(r'$r_{c}/r_{c}(T = 0)$',    fontsize=20)
ax2.set_ylabel(r'$\rho_{c}/\rho_{c}(T = 0)$', fontsize=20)

ax1.set_yscale('log')
ax2.set_yscale('log')

ax1.set_ylim(0.1,1)
ax2.set_ylim(1,1000)

#ax1.set_xlim(0,305)
#ax2.set_xlim(0,305)
ax1.grid(True, which='both')
ax2.grid(True, which='both')
ax3.grid(True, which='both')
plt.savefig("r_p_core.pdf", format='pdf')
#plt.show()
