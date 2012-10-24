#!/bin/env python

# Now the timesteps list will have
# all the dt for each iteration.

import numpy as np
from matplotlib.mlab import load
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import matplotlib.mlab as mlab
files = ['1024_all_191', '1024_all_56', '1024_all_84']
i = 1
for f in files:
    print('Reading file', f)
    data = load(f)
    t  = data[:,0]
    rx, ry, rz = data[:,2],  data[:,3],  data[:,4]
    vx, vy, vz = data[:,5],  data[:,6],  data[:,7]
    ax, ay, az = data[:,8],  data[:,9],  data[:,10]
    jx, jy, jz = data[:,11], data[:,12], data[:,13]

    particle = f.replace('1024_all_','')
    fig = plt.figure(i)
    fig.text(0.5,0.975, particle, horizontalalignment='center', verticalalignment='top')

    # Positions
    pr = fig.add_subplot(221)
    pr.set_title(u'Positions $r$')
    pr.plot(t, rx, label=u'$rx$')
    pr.plot(t, ry, label=u'$ry$')
    pr.plot(t, rz, label=u'$rz$')
    pr.legend(loc='upper right')
    pr.grid(True)

    # Velocities
    pv = fig.add_subplot(222)
    pv.set_title(u'Velocities $(v)$')
    pv.plot(t, vx, label=u'$vx$')
    pv.plot(t, vy, label=u'$vx$')
    pv.plot(t, vz, label=u'$vx$')
    pv.legend(loc='upper right')
    pv.grid(True)

    # Accelerations
    pa = fig.add_subplot(223)
    pa.set_title(u'Accelerations $(a)$')
    pa.plot(t, ax, label=u'$ax$')
    pa.plot(t, ay, label=u'$ay$')
    pa.plot(t, az, label=u'$az$')
    pa.legend(loc='upper right')
    pa.grid(True)

    # Jerks
    pj = fig.add_subplot(224)
    pj.set_title(u'Jerks $(j)$')
    pj.plot(t, jx, label=u'$jx$')
    pj.plot(t, jy, label=u'$jy$')
    pj.plot(t, jz, label=u'$jz$')
    pj.legend(loc='upper right')
    pj.grid(True)

    pa.set_xlabel(u'Time ($N$-body units)')
    pj.set_xlabel(u'Time ($N$-body units)')
    i += 1



plt.show()
