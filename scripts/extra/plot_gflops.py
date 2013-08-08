#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


N = [1024,2048,4096,8192,16384,32768]
gflop = [18.60, 33.65, 51.61, 87.84, 111.30, 143.20]

n = np.array(N)

fig = plt.figure()
f1 = fig.add_subplot(111)

f1.plot(n, gflop, 'o-',  color='red')

f1.set_ylabel(r'Performance [Gflops]', fontsize=15)
f1.set_xlabel(r'$N$', fontsize=15)
f1.set_xlim(10**2.5,4*10**4)
f1.set_ylim(10**1,10**3)
f1.legend(loc='upper left',ncol=1)
f1.axhline(linewidth=2,y=515.2, color='b')

f1.grid(True, which='both')
f1.set_yscale('log')
#f1.set_xscale('log')
plt.show()
