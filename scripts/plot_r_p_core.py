#!/bin/env python

from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import numpy as np
import array
from lib.utils import *
from lib.file_utils import *

# Using epsilon = 0.0001
# Using eta = 0.01

def get_center_of_density(rx,ry,rz,m):
    nj = 6
    p = 0.0
    p_c = np.array([0.0,0.0,0.0])
    dsum = 0.0
    for i in range(0,n):
        d = array.array('d')
        for j in range(0,n):
            if i == j: continue
            drx = rx[j] - rx[i]
            dry = ry[j] - ry[i]
            drz = rz[j] - rz[i]
            r = np.sqrt(drx*drx + dry*dry + drz*drz)
            d.append(r)
        r_sort = np.argsort(d)
        radius = d[r_sort[nj-1]]
        aa = (nj-1) * m[i]
        bb = (4.0 * np.pi *  radius**3)/3.0
        p = aa/bb
        dsum += p
        p_c += np.array([rx[i], ry[i], rz[i]]) * p
    p_c /= dsum
    return p_c

def get_core_radius(rx, ry, rz, m, mp, c):
    # Empty distances array
    d = np.zeros(n)
    # Calculating the distances related to the center of density
    for i in range(n):
        drx = rx[i] - c[0]
        dry = ry[i] - c[1]
        drz = rz[i] - c[2]
        r   = np.sqrt(drx**2 + dry**2 + drz**2)
        d[i] = r
    core_mass = 0.0
    i = 0
    d_sort = np.argsort(d)
    for i in range(n):
        if core_mass  > mp:
            i -= 1
            break
        core_mass += m[d_sort[i]]
    return d[d_sort[i]], i

#############################################################################

t = 10
n = 1024
columns = 15
t_rlx = 20.247

f = get_filename()

summary, snapshot = read_gravidy_output(f, n, t)

# Reading all the positions and velocities
m  = snapshot[:,2]
rx, ry, rz = snapshot[:,3], snapshot[:,4], snapshot[:,5]

# Splitting the data by unit time
m  = np.split(m,  t)
rx = np.split(rx, t)
ry = np.split(ry, t)
rz = np.split(rz, t)

mc = 0.05
tntime = []
tttime = []
radius = []
density = []

print('Time  radius density')
for i in range(0,t):
    c = get_center_of_density(rx[i], ry[i], rz[i], m[i])
    r, nc = get_core_radius(rx[i], ry[i], rz[i], m[i], mc, c)
    p = (nc * m[i][0]) / (4 * np.pi * r**3)

    #print i,i/t_rlx, r, p
    tntime.append(i)
    tttime.append(float(i)/t)
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

#ax1.set_ylim(0.1,1)
#ax2.set_ylim(1,1000)

#ax1.set_xlim(0,305)
#ax2.set_xlim(0,305)
ax1.grid(True, which='both')
ax2.grid(True, which='both')
ax3.grid(True, which='both')
#plt.savefig("r_p_core.pdf", format='pdf')
plt.show()
