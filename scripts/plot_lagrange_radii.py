#!/bin/env python
from lib.utils import *
from lib.file_utils import *

################
f = get_filename()
t    = 10
n    = 1024
t_rh = 20.2
step = 0.1
################

summary, snapshot = read_gravidy_output(f, n, t)


times = []



# Splitting columns of the data
m  = snapshot[:,2]
rx = snapshot[:,3]
ry = snapshot[:,4]
rz = snapshot[:,5]

# Splitting the data by unit time
m  = np.split(m,  t)
rx = np.split(rx, t)
ry = np.split(ry, t)
rz = np.split(rz, t)

tmp_r = []
for i in range(t):
    c = get_center_of_density(n, rx[i], ry[i], rz[i],m[i])
    rc = []
    radius, radius_sort = get_radius(n, rx[i], ry[i], rz[i], c)
    rc = get_core_radius(n,radius, radius_sort, m[i], step)

    times.append(float(i)/t_rh)
    tmp_r.append(rc)

radii = []

for j in range(0,len(times)-1):
    radii.append([i[j] for i in tmp_r])

## Plot
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

fig = plt.figure(1)

ax1 = fig.add_subplot(111)
ax1.set_xlabel(r'$T / T_{\rm rlx} (T = 0)$')
ax1.set_ylabel(r'Lagrange radii')

for ii in range(0,len(times)-1):
    ax1.plot(times, radii[ii],  '-', color='red')

ax1.set_yscale('log')
plt.tight_layout()
plt.show()
