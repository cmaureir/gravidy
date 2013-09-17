#!/usr/bin/env python
from math import sqrt, log, ceil

r = []
v = []
m = []


f = open('04-nbody-p1024_m1.in')

for line in f.readlines():
    line = line.split()
    rx = float(line[1])
    ry = float(line[2])
    rz = float(line[3])
    vx = float(line[4])
    vy = float(line[5])
    vz = float(line[6])

    m.append(float(line[0]))
    r.append([rx, ry, rz])
    v.append([vx, vy, vz])

n = len(r)
acc = [[0.0, 0.0, 0.0] for i in range(n)]
jrk = [[0.0, 0.0, 0.0] for i in range(n)]
e = 1e-4

for i in range(n):
    for j in range(n):
        if i == j:
            continue
        rx = r[j][0] - r[i][0]
        ry = r[j][1] - r[i][1]
        rz = r[j][2] - r[i][2]

        vx = v[j][0] - v[i][0]
        vy = v[j][1] - v[i][1]
        vz = v[j][2] - v[i][2]

        r2 = rx*rx + ry*ry + rz*rz + e*e
        rinv = 1/sqrt(r2)
        r2inv = rinv  * rinv
        r3inv = r2inv * rinv
        r5inv = r2inv * r3inv
        mr3inv = r3inv * m[j]
        mr5inv = r5inv * m[j]

        acc[i][0] += rx * mr3inv
        acc[i][1] += ry * mr3inv
        acc[i][2] += rz * mr3inv

        jrk[i][0] += vx * mr3inv + (3 * vx * rx * rx) * mr5inv
        jrk[i][1] += vy * mr3inv + (3 * vy * ry * ry) * mr5inv
        jrk[i][2] += vz * mr3inv + (3 * vz * rz * rz) * mr5inv

dt = []
eta = 0.01
for i in range(n):
    a2 = acc[i][0]**2 + acc[i][1]**2 + acc[i][2]**2
    j2 = jrk[i][0]**2 + jrk[i][1]**2 + jrk[i][2]**2

    tmp = eta * sqrt(a2/j2)
    tmp = 2**(ceil(log(tmp,2))-1)

    if tmp < 1e-7:
        tmp = 1e-7
    elif tmp > 0.125:
        tmp = 0.125

    dt.append(tmp)

how_many = set(dt)
print(len(how_many))

amount = [0 for i in range(len(how_many))]
k = 0
for i in how_many:
    for j in dt:
        if i == j:
            amount[k] += 1
    k += 1

for i,j in zip(how_many, amount):
    print(i, j)
