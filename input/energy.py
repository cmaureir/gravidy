#!/usr/bin/env python
import math
import sys


f = open(sys.argv[1])

m = []

rx = []
ry = []
rz = []

vx = []
vy = []
vz = []

m_total = 0.0

for line in f.readlines():
    line = line.replace('\t','')
    line = line.replace('\n','')
    values = line.split(" ")
    values = filter(None,values)
    m.append(float(values[0]))
    rx.append(float(values[1]))
    ry.append(float(values[2]))
    rz.append(float(values[3]))
    vx.append(float(values[4]))
    vy.append(float(values[5]))
    vz.append(float(values[6]))
    m_total += float(values[0])

epot_total = 0
ekin_total = 0

epot = 0
ekin = 0

j = 0
for i in range(len(m)):
    epot = 0
    for j in range(len(m)):
        if i != j:
            rxx = rx[j] - rx[i]
            ryy = ry[j] - ry[i]
            rzz = rz[j] - rz[i]
            epot -= (m[i] * m[j]) / math.sqrt(rxx*rxx + ryy*ryy + rzz*rzz + 1e-8)
    ekin = 0.5 * m[i] * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])

    ekin_total += ekin
    epot_total += 0.5* epot

print("Epot:", epot_total)
print("Ekin:", ekin_total)
print("Etotal:", ekin_total+epot_total)


G = 1
n = len(m)
M = m_total

Rv = (-G * M*M  ) / (2 * (epot_total))
U_t = math.sqrt( (Rv * Rv * Rv) / G * M)
t_cr = 2 * math.sqrt(2) * U_t
t_relax = (t_cr * n) / (22 * math.log(0.4 * n))
print("Rv: "  , Rv)
print("U_t: " , U_t)
print("t_cr: ",  t_cr)
print("t_relax: ", t_relax)
