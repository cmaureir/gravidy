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
    rx = data[:,2]
    ry = data[:,3]
    rz = data[:,4]
    vx = data[:,5]
    vy = data[:,6]
    vz = data[:,7]
    ax = data[:,8]
    ay = data[:,9]
    az = data[:,10]
    jx = data[:,11]
    jy = data[:,12]
    jz = data[:,13]


    fig_rv = plt.figure(i)
    fig_rv.text(0.5,0.975,f.replace('1024_all_','Particle '), horizontalalignment='center', verticalalignment='top')
    # Positions
    prx = fig_rv.add_subplot(321)
    prx.plot(t, rx, label=u'$rx$')
    prx.set_title(u'Position (rx)')

    pry = fig_rv.add_subplot(323)
    pry.plot(t, rx, label=u'$rx$')
    pry.set_title(u'Position (ry)')

    prz = fig_rv.add_subplot(325)
    prz.plot(t, rx, label=u'$rx$')
    prz.set_title(u'Position (rz)')

    # Velocities
    pvx = fig_rv.add_subplot(322)
    pvx.plot(t, vx, label=u'$vx$')
    pvx.set_title(u'Velocity (vx)')

    pvy = fig_rv.add_subplot(324)
    pvy.plot(t, vy, label=u'$vx$')
    pvy.set_title(u'Velocity (vy)')

    pvz = fig_rv.add_subplot(326)
    pvz.plot(t, vz, label=u'$vx$')
    pvz.set_title(u'Velocity (vz)')

    i += 1

    fig_aj = plt.figure(i)
    fig_aj.text(0.5,0.975,f.replace('1024_all_','Particle '), horizontalalignment='center', verticalalignment='top')
    # Positions
    pax = fig_aj.add_subplot(321)
    pax.plot(t, ax, label=u'$ax$')
    pax.set_title(u'Acceleration (ax)')

    pay = fig_aj.add_subplot(323)
    pay.plot(t, ay, label=u'$ay$')
    pay.set_title(u'Acceleration (ay)')

    paz = fig_aj.add_subplot(325)
    paz.plot(t, az, label=u'$az$')
    paz.set_title(u'Acceleration (az)')

    # Velocities
    pjx = fig_aj.add_subplot(322)
    pjx.plot(t, jx, label=u'$jx$')
    pjx.set_title(u'Jerk (jx)')

    pjy = fig_aj.add_subplot(324)
    pjy.plot(t, jy, label=u'$jy$')
    pjy.set_title(u'Jerk (jy)')

    pjz = fig_aj.add_subplot(326)
    pjz.plot(t, jz, label=u'$jz$')
    pjz.set_title(u'Jerk (jz)')

    i += 1



plt.show()
