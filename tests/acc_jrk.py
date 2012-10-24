#!/usr/bin/env python

from matplotlib.mlab import load
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax.grid(color='gray')
ax2.grid(color='gray')
ax.set_ylabel(u'Acceleration, $a_{83}$ (N-body units)', fontsize=20)
ax2.set_ylabel(u'Jerk, $a^{(1)}_{83}$ (N-body units)', fontsize=20)
ax2.set_xlabel(u'Integration time (N-body units)', fontsize=20)
ax.set_title(u'Acceleration and Jerk of the particle 83\n$a_{83}$ and $a^{(1)}_{83}$', fontsize=20)


f = open('83-acc_jrk')
data = load(f)
t  = data[:,0] # Data from column 1
acc_x = data[:,1] # Data from column 2
acc_y = data[:,2] # Data from column 3
acc_z = data[:,3] # Data from column 4

jrk_x = data[:,4] # Data from column 2
jrk_y = data[:,5] # Data from column 3
jrk_z = data[:,6] # Data from column 4

ax.plot(t, acc_x, label=u'$a_{x}$')
ax.plot(t, acc_y, label=u'$a_{y}$')
ax.plot(t, acc_z, label=u'$a_{z}$')

ax2.plot(t, jrk_x, label=u'$a^{(1)}_{x}$')
ax2.plot(t, jrk_y, label=u'$a^{(1)}_{y}$')
ax2.plot(t, jrk_z, label=u'$a^{(1)}_{z}$')

ax.legend(loc='upper right')
ax2.legend(loc='upper right')
plt.show()
