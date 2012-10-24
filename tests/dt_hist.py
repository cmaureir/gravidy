#!/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import matplotlib.mlab as mlab
from matplotlib.mlab import load

f = open('dt_0.00125')
data = load(f)
times = data[:,1]
times0 = times[:128]
times1 = times[128:]

fig = plt.figure()

# t = 0, timesteps distribution
ax = fig.add_subplot(211)
n, bins, patches = ax.hist(times0, bins=len(times0)*10, facecolor='green', alpha=0.75)
bincenters = 0.5*(bins[1:]+bins[:-1])
y = mlab.normpdf( bincenters, 100, 50)
l = ax.plot(bincenters, y, 'r--', linewidth=1)

# t = 1, timesteps distribution
ax2 = fig.add_subplot(212)
n, bins, patches = ax2.hist(times1, bins=len(times1)*10, facecolor='green', alpha=0.75)
bincenters = 0.5*(bins[1:]+bins[:-1])
y = mlab.normpdf( bincenters, 100, 50)
l = ax2.plot(bincenters, y, 'r--', linewidth=1)

ax.grid(color='gray', linestyle='-', linewidth=0.5)
ax2.grid(color='gray', linestyle='-', linewidth=0.5)

ax.set_ylabel('Amount of timesteps at $t=0$')
ax2.set_ylabel('Amount of timesteps at $t=1$')

ax.set_title(u'Timesteps distribution\n$(t = 0$ and $t = 1)$')
ax2.set_xlabel('Timesteps values')

ax.set_xscale('log')
ax2.set_xscale('log')

ax.set_ylim(0,80)
ax2.set_ylim(0,80)
plt.show()

