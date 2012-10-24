#!/bin/env python

f = open('128.times')

timesteps = []
tmp = []

for line in f.readlines():
    line = line.split()
    part_id = int(line[0])
    part_dt = float(line[1])

    if part_id == 0:
        timesteps.append(tmp)
        tmp = []

    tmp.append(part_dt)

timesteps = timesteps[1:]



# Now the timesteps list will have
# all the dt for each iteration.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import matplotlib.mlab as mlab

fig = plt.figure()

# t = 0, timesteps distribution
ax = fig.add_subplot(211)
n, bins, patches = ax.hist(timesteps[0], bins=len(timesteps[0])*10, facecolor='green', alpha=0.75)
bincenters = 0.5*(bins[1:]+bins[:-1])
y = mlab.normpdf( bincenters, 100, 50)
l = ax.plot(bincenters, y, 'r--', linewidth=1)

# t = 1, timesteps distribution
ax2 = fig.add_subplot(212)
n, bins, patches = ax2.hist(timesteps[-1], bins=len(timesteps)*10, facecolor='green', alpha=0.75)
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

ax.set_ylim(0,40)
ax2.set_ylim(0,40)
plt.show()

