#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import matplotlib.mlab as mlab

files = []
files.append('1024_dt_0t.out')
files.append('1024_dt_1t.out')
files.append('1024_dt_2t.out')
files.append('1024_dt_3t.out')

fig = plt.figure()

plots = len(files)
i = 0
for name in files:
    print(name)
    dt = []
    f = open(name)
    for line in f.readlines():
        line = line.split()
        dt.append(int(line[3]))

    subplot = plots*100 + 10 + (i+1)
    ax = fig.add_subplot(subplot)

    n, bins, patches = ax.hist(dt, facecolor='green', alpha=0.75)
    bincenters = 0.5*(bins[1:]+bins[:-1])
    y = mlab.normpdf( bincenters, 100, 50)
    ax.plot(bincenters, y, 'r--', linewidth=1)

    ax.grid(color='gray', linestyle='-', linewidth=0.5)
    ax.set_ylabel('Amount of timesteps at $t='+str(i)+'$')
    #ax.set_title(u'Timesteps distribution (1024)\n($t = 0$, $t = 1$ and $t = 2$)')
    ax.set_xlabel('Exponent of the Timesteps value ($2^x$)')
    ax.set_xlim(-14, 0)
    ax.set_ylim(0,400)

    i += 1


plt.show()
#plt.savefig('1024_dt_1_2.png')
