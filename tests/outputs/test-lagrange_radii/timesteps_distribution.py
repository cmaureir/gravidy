#!/bin/env python

from matplotlib.mlab import load
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import matplotlib.mlab as mlab
from matplotlib import rc
from math import log
import numpy as np
import array

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def get_exponent(data):
    return np.array([int(np.log(i)/np.log(2)) for i in data])
t = 701
n = 4096
columns = 3
t_rh = 59.689350
files = []
files.append('files/4k_dt.out')

m = np.zeros(t*n)
rx, ry, rz = np.zeros(n), np.zeros(n), np.zeros(n)
vx, vy, vz = np.zeros(n), np.zeros(n), np.zeros(n)

j = 0
# Reading files and getting the data
for f in files:
    data   = np.fromfile(file=open(f), dtype=float, sep=" ").reshape((t*n, columns))

    it = data[:,0]
    dt = data[:,2]

    it = np.split(it, t)
    dt = np.split(dt, t)

    dt[0] = get_exponent(dt[0])
    dt[1] = get_exponent(dt[1])
    dt[60] = get_exponent(dt[60])

    bins_amount = 14

    fig = plt.figure()
    color = '#62A9FF'

    # Simple
    ax1 = fig.add_subplot(211)
    ax1.hist(dt[0], bins=bins_amount, facecolor=color, alpha=0.9)

    ax3 = fig.add_subplot(212)
    ax3.hist(dt[60], bins=bins_amount, facecolor=color, alpha=0.9)

    ax1.grid(color='gray', linestyle='-', linewidth=0.5)
    ax3.grid(color='gray', linestyle='-', linewidth=0.5)

    ax3.set_xlabel(r'Power of $2$ Time-step', fontsize=14)

    ax1.set_ylim(0,900)
    ax3.set_ylim(0,900)
    ax1.set_xlim(-18,0)
    ax3.set_xlim(-18,0)

    fig.text(0.05 , 0.55 , r'Particles' , horizontalalignment='center' ,
            verticalalignment='top' , rotation=90 , fontsize=15)
    fig.text(0.93 , 0.85 , r'$T/[ T_{\rm rlx}(T=0) ] = 0$', horizontalalignment='center' ,
            verticalalignment='top' , rotation=-90, fontsize=15)
    fig.text(0.93 , 0.4 , r'$T/[ T_{\rm rlx}(T=0) ] = 1$', horizontalalignment='center' ,
            verticalalignment='top' , rotation=-90, fontsize=15)

    plt.show()
