#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

data_n = []
import subprocess

x_upper = 10**5
y_upper = 10**4
n = np.array([1024, 2048,4096,8192,16384,32768])
nn = np.array([i for i in range(1,x_upper)])
fig = plt.figure()
f1 = fig.add_subplot(111)
time = np.array([1.21, 3.22, 9.45, 23.31, 82.63, 275.53])
f1.plot(n , time                 , 'o-')
f1.plot(nn     , 3e-8*(nn**3)   , '--'   , color='green')
f1.plot(nn     , 9e-7*(nn**2)    , '--'   , color='cyan')
f1.plot(nn     , 1e-4*nn*np.log(nn) , '--'   , color='red')
f1.set_ylabel(r'Clock time $(s)$', fontsize=15)
f1.set_xlabel(r'$N$', fontsize=15)
f1.set_xlim(1e3,x_upper)
f1.set_ylim(1,y_upper)

f1.text(10**4.7, 10**(1.5),r'$N\log N$', style='italic',
        bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
f1.text(10**(3.5), 10**(3.5),r'$N^{3}$', style='italic',
        bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
f1.text(10**(4.4), 10**(3),r'$N^{2}$', style='italic',
        bbox={'facecolor':'cyan', 'alpha':0.5, 'pad':10})


f1.grid(True)
f1.set_yscale('log')
f1.set_xscale('log')


plt.savefig("test-time.pdf", format='pdf')
#plt.show()
