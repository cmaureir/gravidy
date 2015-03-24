import numpy as np
import sys

e = float(sys.argv[1])

r1 = np.zeros(3)
r2 = np.zeros(3)
v1 = np.zeros(3)
v2 = np.zeros(3)

m1 = 1.0                     # set masses
m2 = 1.0
mtot = m1+m2

m1 = m1/mtot                    # normalise to mtot = 1
m2 = m2/mtot

q = m1/m2                       # mass ratio

atot = 1.
a1 = atot / (1. + q)         # initial positions (a1 + a2 = a)
a2 = q * a1

a = np.array([a1,a2])
G = 1.


r = a * (1. - e)

#xbh = np.array([r[0], -r[1]])
r1[0] = -r[0]
r2[0] = r[1]

v = np.sqrt((1.+e)/(1.-e))

v1[1] = v * m2
v2[1] = v * -m1

print(m1, r1, v1)
print(m2, r2, v2)



G = 1
mu_std = G * (m1 + m2)
mu = (m1 * m2) / (m1 + m2)
com_v = (v1 * m1 + v2 * m2)/(m1+m2)

r = r2 - r1
v = v2 - v1


j = np.cross(r, v)
rr = np.linalg.norm(r)
vv = np.linalg.norm(v)
jj = np.linalg.norm(j)

rinv = 1.0/rr

j2 = jj*jj
v2 = vv*vv
m1m2 = m1 * m2

espec = v2 * 0.5 - mu_std/rr
semimajor = -mu_std / (2*espec)
ecc = np.sqrt(1.0+2.0*espec*j2/(mu_std*mu_std))

print("espec=", espec)
print("bind=", 0.5 * mu * v2 - (m1m2 * G * rinv))
print("a=", semimajor)
print("ecc=", ecc)
print("CoMv=", com_v)
