#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)


def print_lists():
    global cpu_time
    global gpu_time
    global cpu_ite
    global gpu_ite
    print(cpu_time)
    print(gpu_time)
    print(cpu_ite)
    print(gpu_ite)

# Fixed alpha values
alpha      = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
cpu_ite  = dict()
gpu_ite  = dict()
gpu_time = dict()
cpu_time = dict()

data_n = []
import subprocess

p1 = subprocess.Popen(['ls', '-1','-I','*.py','-I','*.gpu'], stdout=subprocess.PIPE)

for directory in p1.stdout:
    directory = directory.strip()
    if directory == None or directory == "": break
    d = int(directory)
    data_n.append(d)
    cpu_ite[d], gpu_ite[d]  = np.zeros(len(alpha)), np.zeros(len(alpha))
    gpu_time[d], cpu_time[d] = np.zeros(len(alpha)), np.zeros(len(alpha))
    p2 = subprocess.Popen(['ls', '-1',directory], stdout=subprocess.PIPE)
    pos = 0
    for filename in p2.stdout:
        filename = directory + "/" + filename.strip()
        f = open(filename)
        f.readline()
        f.readline()
        line = f.readline()
        line = line.strip().split()
        cpu_ite[d][pos]  = float(line[2])
        gpu_ite[d][pos]  = float(line[3])
        gpu_time[d][pos] = float(line[6])
        cpu_time[d][pos] = float(line[7])-float(line[6])
        pos += 1

#cpu_ite    = np.array([30016, 29952, 29951, 29925, 27903, 27110, 24610, 0])
#gpu_ite    = np.array([1985, 2049, 2050, 2076, 4098, 4891, 7391, 32001])
#
#gpu_time   = np.array([685.07, 698.26, 698.35, 703.21, 1079.38, 1226.08, 1690.54, 5235.25])
#total_time = np.array([1579.38, 1581.23, 1570.49, 1550.01, 1459.49, 1546.96, 1860.31, 5244.68])

fig = plt.figure()
f1 = fig.add_subplot(111)
for i in data_n:
    y = cpu_time[i] + gpu_time[i]
    #y = cpu_time[i]
    #print(y)
    f1.plot(alpha, y, '*-',label=r'$N='+str(i)+'$')
f1.set_ylabel(r'Clock time $(s)$', fontsize=15)
f1.set_xlabel(r'Alpha Parameter $(\alpha)$', fontsize=15)
#f1.set_xlim(0,0.055)
#f1.set_ylim(-(10**1),10**5)
#f1.legend(loc='lower right',ncol=2)
f1.legend(loc='lower right',ncol=2)
f1.grid(True, which='both')
f1.set_yscale('log')
#f1.set_xscale('log')
#f1.annotate('Best', xy=(alpha[4], total_time[4]), xytext=(-20,20),
#            textcoords='offset points', ha='center', va='bottom',
#            bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),
#            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',
#                            color='red'))
plt.show()
