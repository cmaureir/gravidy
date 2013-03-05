#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)


# Fixed alpha values
alpha      = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

data_n = []
import subprocess

p1 = subprocess.Popen(['ls', '-1','-I','*.py'], stdout=subprocess.PIPE)

time = []
for directory in p1.stdout:
    directory = directory.strip()
    if directory == None or directory == "": break
    d = int(directory)
    data_n.append(d)
    p2 = subprocess.Popen(['ls', '-1',directory], stdout=subprocess.PIPE)
    pos = 0
    for filename in p2.stdout:
        filename = directory + "/" + filename.strip()
        f = open(filename)
        f.readline()
        f.readline()
        line = f.readline()
        line = line.strip().split()
        if pos == 4:
            time.append(float(line[7]))
        pos += 1

x_upper = 10**5
y_upper = 10**5
n = np.array(data_n)
nn = np.array([i for i in range(1,x_upper)])
fig = plt.figure()
f1 = fig.add_subplot(111)
f1.plot(data_n , time                 , 'o-')
f1.plot(nn     , 0.0000048*(nn**2)    , '--'   , color='cyan')
f1.plot(nn     , 0.00000003*(nn**3)   , '--'   , color='green')
f1.plot(nn     , 0.0007*nn*np.log(nn) , '--'   , color='red')
f1.set_ylabel(r'Clock time $[sec]$', fontsize=15)
f1.set_xlabel(r'$N$', fontsize=15)
f1.set_xlim(900,x_upper)
f1.set_ylim(5,y_upper)
#f1.legend(loc='lower right',ncol=2)

f1.text(10**4, 10**(1.5),r'$N\log N$', style='italic',
        bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
f1.text(10**(3.5), 10**(3.5),r'$N^{3}$', style='italic',
        bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})
f1.text(10**(4.3), 10**(3),r'$N^{2}$', style='italic',
        bbox={'facecolor':'cyan', 'alpha':0.5, 'pad':10})


f1.grid(True, which='both')
f1.set_yscale('log')
f1.set_xscale('log')
plt.show()
