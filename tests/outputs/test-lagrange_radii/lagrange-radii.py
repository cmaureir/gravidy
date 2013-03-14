#!/bin/env python
import numpy as np
import array

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
    d_sort = np.argsort(d)
    core_mass = 0.0
    rr = []
    for i in range(n):
        if core_mass > mp:
            mp += 0.05
            rr.append(d[d_sort[i-1]])
        #    print("#", i-1, mp, core_mass,rr)
        core_mass += m[d_sort[i]]

    return rr

t = 1001
n = 1024
columns = 15
t_rh = 59.689350
files = []
#files.append('files/radii_4k.out')
files.append('1k_1000t.out')

m = np.zeros(t*n)
rx, ry, rz = np.zeros(n), np.zeros(n), np.zeros(n)
vx, vy, vz = np.zeros(n), np.zeros(n), np.zeros(n)

j = 0
# Reading files and getting the data
for f in files:
    data   = np.fromfile(file=open(f), dtype=float, sep=" ").reshape((t*n, columns))

    # Reading all the positions and velocities
    #m  = data[:,1]
    m  = [[0.0009765625 for i in range(n)] for j in range(t)]
    rx, ry, rz = data[:,2], data[:,3], data[:,4]

    # Splitting the data by unit time
    #m  = np.split(m,  t)
    rx = np.split(rx, t)
    ry = np.split(ry, t)
    rz = np.split(rz, t)

    for i in range(t):
        c = get_center_of_density(rx[i], ry[i], rz[i],m[i])
        rc = []
        mc = 0.05
        rc = get_core_radius(rx[i], ry[i], rz[i], m[i], mc, c)
        print(i/t_rh, c, rc)
    j += 1
