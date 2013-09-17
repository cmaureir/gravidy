#!/bin/env python
from lib.utils import *
from lib.file_utils import *


f = get_filename()

################
t    = 10
n    = 1024
t_rh = 20.2
step = 0.1
################

summary, snapshot = read_gravidy_output(f, n, t)

time = summary[:,1]
cum_error =  summary[:,11]


## Plot
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

fig = plt.figure(1)

ax1 = fig.add_subplot(111)
ax1.set_ylabel(r'$|\frac{E_{\rm tot}(t) - E_{\rm tot}(0)}{E_{\rm tot}(0)}|$')
ax1.set_xlabel(r'$T [NBU]$')

ax1.plot(time, cum_error,  '-', color='red')

ax1.set_yscale('log')
plt.tight_layout()
plt.show()
