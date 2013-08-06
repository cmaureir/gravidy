#!/bin/env python
from lib.utils import *
from lib.file_utils import *

################
f    = '1k.in'
t    = 9
n    = 1024
t_rh = 20.2
################

summary, snapshot = read_gravidy_output(f, n, t)

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

for i in range(t):
    c = get_center_of_density(n, rx[i], ry[i], rz[i],m[i])
    rc = []
    radius, radius_sort = get_radius(n, rx[i], ry[i], rz[i], c)
    rc = get_core_radius(n,radius, radius_sort, m[i], 0.1)
    print(i/t_rh, c, rc)
