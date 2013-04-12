#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)


N = [1024,2048,4096,8192,16384,32768]
cpu_time = [
12.4800,
56.1800,
275.7700,
1159.9501,
5534.8301,
24973.0605]
]

gpu_time = [
5.66,   # peor 0.01
16.69,  # peor 0.01
65.37,  # peor 0.01
203.84, # peor 0.001
872.20, # mejor 0.001
3695.2800] # mejor 0.001

#N = [1024,2048,4096,8192,16384,32768,65536]

n = np.array(N)
acc = np.array(cpu_time) / np.array(gpu_time)

fig = plt.figure()
f1 = fig.add_subplot(111)
f1.plot(n, acc, 'o-')
f1.set_ylabel(r'Acceleration (CPU/GPU)', fontsize=15)
f1.set_xlabel(r'$N$', fontsize=15)
#f1.set_xlim(900,x_upper)
f1.set_ylim(0,10)
#f1.legend(loc='lower right',ncol=2)

f1.grid(True, which='both')
#f1.set_yscale('log')
#f1.set_xscale('log')
plt.show()
