#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


N = [1024,2048,4096,8192,16384,32768]
cpu_time = [
   12.98,
   61.32,
  282.98,
 1227.40,
 5542.35,
26383.71]

omp_time = [
  8.19,
 34.94,
162.64,
682.56,
3227.91,
15076.40] # CAMBIAr

cpu_gpu_time = [
   3.57,
  13.42,
  54.28,
 208.91,
 904.82,
3722.92]

mpi1_time = [
   6,
  14,
  51,
 105,
 364,
1247]

mpi2_time = [
   2,
   7,
  27,
  64,
 317,
1145]


gpu_time = [
  1.21,
  3.22,
  9.45,
 23.31,
 82.63,
275.53]

n = np.array(N)
acc1 = np.array(cpu_time) / np.array(omp_time)
acc2 = np.array(cpu_time) / np.array(cpu_gpu_time)
acc3 = np.array(cpu_time) / np.array(mpi1_time)
acc4 = np.array(cpu_time) / np.array(mpi2_time)
acc5 = np.array(cpu_time) / np.array(gpu_time)

fig = plt.figure()
f1 = fig.add_subplot(111)

f1.plot(n, acc1, 'o-',  color='red'  , label=r'(OpenMP)/CPU')
f1.plot(n, acc2, 'bs-', color='blue' , label=r'(CPU + GPU)/CPU')
f1.plot(n, acc3, 'v-', color='yellow' , label=r'(MPI-1)/CPU')
f1.plot(n, acc4, 'd-', color='orange' , label=r'(MPI-2)/CPU')
f1.plot(n, acc5, 'g^-', color='green', label=r'GPU/CPU')

f1.set_ylabel(r'Acceleration', fontsize=15)
f1.set_xlabel(r'$N$', fontsize=15)
#f1.set_xlim(900,x_upper)
f1.set_ylim(0.5,10**(3))
f1.legend(loc='upper left',ncol=1)

f1.grid(True, which='both')
f1.set_yscale('log')
#f1.set_xscale('log')
plt.show()
