#!/usr/bin/env python
import csv
import operator
from math import sqrt, log
from numpy import array, sort

# Obtained from "rot"
dx, dy, dz = -2.030242657157e-01, -7.161864541739e-02, -1.099350671471e-01



m = []
rx = []
ry = []
rz = []
vx = []
vy = []
vz = []

# Reading files

def import_text(filename, separator):
    for line in csv.reader(open(filename), delimiter=separator,
                           skipinitialspace=True):
        if line:
            yield line

for data in import_text('fort.10.mod', '/'):
    t = [x for x in filter(None, data[0].split(' '))]
    m.append(float(t[0]))
    rx.append(float(t[1]))
    ry.append(float(t[2]))
    rz.append(float(t[3]))
    vx.append(float(t[4]))
    vy.append(float(t[5]))
    vz.append(float(t[6]))

# End

N = len(m)
M = sum(m)

distances = {}

for i in range(N):
    tmp = sqrt((dx - rx[i])**2 + (dy - ry[i])**2 + (dz - rz[i])**2)
    distances[i] = tmp

sorted_d = sorted(distances.items(), key=operator.itemgetter(1))

half_m = 0
j = 0
rhx = rhy = rhz = 0
for i in range(N):
    if half_m == M/2:
        j = i
        break
    half_m = half_m + m[sorted_d[i][0]]


## Starting relaxation time calculation
# half-mass radius
R_h = sorted_d[j-1][1]
print('R_h:',R_h)

# All the bodies with the same mass
G = 1
n = 128
m = 0.007812500000000
mean_m = m

R_h3 = R_h ** 3
B = 0.11;

# Constant
C = (n * R_h3)/mean_m

# B = 1/log_lambda ~ 0.11
t_relax = 0.138 * sqrt(C) * B
print('Relaxation time:',t_relax, '[nbody units]')

M = 1
energy_tot = -0.285491129674071

# nbody units
Rv = (-G * m * m) / (4 * energy_tot)

# 1 length unit     = 4.964492e-03 pc

# physics units
Rv = Rv * 4.964492e-03

# Transformation from nbody units to Myr
t_star = 14.94 * sqrt(Rv**3.0)/M

print('Star time:', t_star, '[Myr]')

