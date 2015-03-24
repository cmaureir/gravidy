import sys
from matplotlib import use, rcParams, rc, ticker
use('Agg')
import matplotlib.pyplot as plt
import numpy as np
rcParams.update({'figure.autolayout': True})
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Palatino'], 'size': 8})
rc('figure', figsize=(8,9))


#fig, ax = plt.subplots()
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

files = []
files.append('n1_t5000_e0p99')
files.append('n3_t5000_e0p99')
info = [['-', r'$PEC$'],
        ['--', r'$P(EC)^3$']]
for i, f in enumerate(files):
    data = np.genfromtxt(f)

    time = data[:,0]/np.pi
    semimajor = data[:,1]
    ecc = data[:,2]
    distance = data[:,3]
    energy = data[:,4]

    mk, lb = info[i%len(files)]
    ax1.plot(time[::20], semimajor[::20], linestyle=mk, label=lb)
    ax2.plot(time[::40], ecc[::40],       linestyle=mk, label=lb)
    ax3.plot(time[::20], energy[::20],    linestyle=mk, label=lb)


ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')

ax1.legend(loc='lower right', fontsize=12)
ax2.legend(loc='lower right', fontsize=12)
ax3.legend(loc='lower right', fontsize=12)

#ax.legend(loc='upper right', fontsize=12)
#ax.set_xlim(0,60)
#ax.set_ylabel(r'$\Delta E$', fontsize=20)
#ax.set_ylabel(r'$\Delta a$', fontsize=20)
ax1.set_ylabel(r'$\Delta a$', fontsize=20)
ax2.set_ylabel(r'$\Delta e$', fontsize=20)
ax3.set_ylabel(r'$\Delta E$', fontsize=20)

ax1.set_title(r'Binary $a = 1,\ e = 0.1$', fontsize=20)
ax3.set_xlabel(r'$T(\Omega_{0})^{-1}$', fontsize=20)
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
fig.savefig('fig.pdf', format='pdf')
#fig.show()
