#!/usr/bin/env python

from matplotlib.mlab import load
import matplotlib.pyplot as plt

files = []
files.append('dt_1k_84')
#files.append('dt_1k_191')
files.append('dt_1k_214')
#files.append('dt_1k_268')
files.append('dt_1k_498')
#files.append('dt_1k_565')
files.append('dt_1k_581')
#files.append('dt_1k_651')
files.append('dt_1k_743')


#markers = ['--','_','-.']
#markers = ['o-','*-', 'x-', '<-', '>-', '^-', 'p-']

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_yscale('log', linewidth=2)
ax.grid(color='gray')
ax.set_ylabel(u'Timestep of particle $i$, $dt_{i}$ (N-body units)', fontsize=20)
ax.set_xlabel(u'Integration time (N-body units)', fontsize=20)
ax.set_title(u'Timesteps evolution $(1024\ particles)$', fontsize=20)


i = 0
for f in files:
    print('Reading file', f)
    data = load(f)
    t  = data[:,0] # Data from column 2
    dt = data[:,2] # Data from column 5
    ax.plot(t, dt, label=u'$i='+f.replace('dt_1k_','')+'$')
    print(len(t), len(dt))
    i += 1

ax.legend(loc='upper right')
#ax.set_ylim(1e-6, 1e-2)
plt.show()
