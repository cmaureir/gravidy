#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)


# Fixed etas values
etas = [0.001, 0.003, 0.005, 0.007, 0.01, 0.03, 0.05, 0.07, 0.1]

ttime = dict()
rel_e = dict()

data_n = []
import subprocess
import re

p1 = subprocess.Popen(['ls', '-1','-I','*.py'], stdout=subprocess.PIPE)

for directory in p1.stdout:
    directory = directory.strip()
    if directory == None or directory == "": break
    d = int(directory)
    data_n.append(d)
    ttime[d], rel_e[d]  = np.zeros(len(etas)), np.zeros(len(etas))
    p2 = subprocess.Popen(['ls', '-1',directory], stdout=subprocess.PIPE)
    pos = 0
    for filename in p2.stdout:
        filename = directory + "/" + filename.strip()
        #s = re.search('e[0-9].[0-9]+',filename).group(0).replace('e','')
        f = open(filename)
        f.readline()
        f.readline()
        line = f.readline()
        line = line.strip().split()
        ttime[d][pos]  = float(line[7])
        rel_e[d][pos]  = float(line[9])
        pos += 1

f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
ax = [ax1, ax2, ax3, ax4, ax5, ax6]
for i,j in zip(ax,data_n):
    y = ttime[j]
    i.plot(etas, y, '*-',label=r'$n='+str(j)+'$')
    i.legend(loc='upper right')
    i.set_yscale('log')
    i.set_ylim(10**0,10**5)
    i.grid(True)

ax3.set_ylabel(r'Clock time ($t$)', fontsize=18)
ax5.set_xlabel(r'$\eta$', fontsize=15)
ax6.set_xlabel(r'$\eta$', fontsize=15)
#f1.set_ylabel(r'$\Delta E / E_{\rm t=0}', fontsize=15)
#f1.set_xlabel(r'eta Parameter $(\eta)$', fontsize=15)
##f1.set_xlim(0,0.055)
##f1.set_ylim(-(10**1),10**5)
plt.show()
