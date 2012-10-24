#!/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr


# Test information
#   Number of particles
#   Execution time
particles = [128,256,512,1024,2048, 4096, 8192, 16384]
cpu_times = [0.220, 1.103, 4.213, 19.802, 90.494, 429.872, 1867.562, 8803.11]

# Scales to use the BigO slopes
t  = np.array([i      for i in range(1,particles[-1])])
t2 = np.array([0.005*i for i in range(1,particles[-1])])
t3 = np.array([0.005*i for i in range(1,particles[-1])])

tt  = t[:4200]
tt3 = t3[:4200]

# Values of the nlogn slope
nlogn = [0.005*i*np.log(i) for i in t]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.grid(color='gray', linestyle='-', linewidth=0.5)
ax.set_ylabel('Execution time (s)', fontsize=20)
ax.set_xlabel('Particles', fontsize=20)

# Execution time
ax.plot(particles, cpu_times, 'o-', label='GraviDy')

# BigO notation plots
ax.plot(t, nlogn, '--')
ax.plot(t, t2*t2, '--')
ax.plot(tt, tt3**3, '--')

ax.text(1500, 4000, 'n^3',  fontsize=20)
ax.text(15000, 4000, 'n^2', fontsize=20)
ax.text(12000, 1000,  'n log n', fontsize=20)

# Limits according the test results
ax.set_ylim(0,9000)
ax.set_xlim(0,particles[-1])
ax.legend(loc='center')
ax.set_title('Execution time')
ax.grid(True)

from XKCDify import *
XKCDify(ax, xaxis_loc=0.0, yaxis_loc=1.0,
        xaxis_arrow='+-', yaxis_arrow='+-',
        expand_axes=True)
plt.show()

