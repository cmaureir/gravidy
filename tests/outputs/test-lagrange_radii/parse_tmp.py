#!/usr/bin/env python

from matplotlib.mlab import load
from matplotlib import rc
import matplotlib.pyplot as plt

rc('text', usetex=True)
fig = plt.figure()
ax = fig.add_subplot(111)



time = []
radius = []
f = "tmp_radii"
f = open(f)
for line in f.readlines():
    if line == None: break
    line = line.strip()
    line = line.split(",")

    t = float(line[0].replace('(',''))
    r = float(line[4])
    time.append(t)
    radius.append(r)
#    print(t,r)

ax.plot(time, radius, label=u'r', color='red')

plt.show()
