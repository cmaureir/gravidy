#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


etas = [0.001, 0.003, 0.005, 0.007, 0.01, 0.03, 0.05, 0.07, 0.1]

n = { 1024:[[
                3.5300,
                2.0500,
                1.7500,
                1.4100,
                1.1800,
                0.6800,
                0.5700,
                0.5100,
                0.4300
                ],[
                3.069933e-07,
                3.068267e-07,
                3.055612e-07,
                3.053740e-07,
                3.093654e-07,
                3.603642e-07,
                6.057363e-07,
                1.574985e-06,
                1.894220e-06
                ]],
      2048:[[
                9.5200,
                5.5600,
                4.3000,
                3.7400,
                3.0700,
                1.8800,
                1.4100,
                1.2400,
                1.1400
                ],[
                6.166337e-08,
                6.161530e-08,
                6.072894e-08,
                6.177926e-08,
                5.905466e-08,
                1.610511e-08,
                4.991711e-07,
                1.015560e-06,
                1.543131e-06
                ]],
      4096:[[
                29.6600,
                17.0800,
                13.3200,
                11.2600,
                9.3800 ,
                5.7300 ,
                4.3600 ,
                3.8400 ,
                3.2200
                ],[
                1.081358e-08,
                1.044539e-08,
                9.508120e-09,
                5.306622e-09,
                2.781821e-09,
                2.350845e-08,
                9.428625e-08,
                6.214629e-07,
                1.933086e-06
                ]],
      8192:[[
                71.4700,
                41.3400,
                32.2600,
                27.5000,
                22.9800,
                13.5100,
                10.4900,
                 8.8600,
                 7.6000
                ],[
                 6.244570e-08,
                 6.160458e-08,
                 6.058871e-08,
                 5.688093e-08,
                 5.149197e-08,
                 3.391105e-08,
                 3.215287e-07,
                 2.610568e-07,
                 3.699080e-06
                ]],
     16384:[[
                273.4300,
                149.5700,
                115.7900,
                 97.6200,
                 82.6300,
                 47.6700,
                 37.0400,
                 31.7200,
                 26.6700
                ],[
                2.784691e-07,
                2.770324e-07,
                2.740612e-07,
                2.681489e-07,
                2.568019e-07,
                2.754555e-07,
                1.129782e-06,
                6.297768e-07,
                8.294093e-05
                ]],
     32768:[[
                868.9500,
                500.5500,
                387.5900,
                330.9700,
                275.3600,
                160.7700,
                124.3600,
                105.8200,
                89.9400
                ],[
                8.210526e-07,
                8.231823e-07,
                8.276378e-07,
                8.338627e-07,
                8.460826e-07,
                8.148497e-07,
                2.810966e-07,
                7.650663e-05,
                2.910441e-04
                ]]
                }

f1, ((f1_ax1, f1_ax2), (f1_ax3, f1_ax4), (f1_ax5, f1_ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
f1_ax = [f1_ax1, f1_ax2, f1_ax3, f1_ax4, f1_ax5, f1_ax6]

i = 0
for j in sorted(n.keys()):
    ttime = n[j][0]
    f1_ax[i].plot(etas, ttime, '*-', label=r'$N='+str(j)+'$')
    f1_ax[i].legend(loc='upper right')
    f1_ax[i].set_yscale('log')
    #f1_ax[i].set_xscale('log')
    f1_ax[i].set_ylim(10**(-1),10**4)
    f1_ax[i].set_xlim(0,0.11)
    f1_ax[i].grid(True)
    i += 1

f1_ax3.set_ylabel(r'Clock time $[sec]$', fontsize=16)
f1_ax5.set_xlabel(r'$\eta$', fontsize=16)
f1_ax6.set_xlabel(r'$\eta$', fontsize=16)


f2, ((f2_ax1, f2_ax2), (f2_ax3, f2_ax4), (f2_ax5, f2_ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
f2_ax = [f2_ax1, f2_ax2, f2_ax3, f2_ax4, f2_ax5, f2_ax6]

i = 0
for j in sorted(n.keys()):
    energy = n[j][1]
    f2_ax[i].plot(etas, energy, '*-', label=r'$N='+str(j)+'$')
    f2_ax[i].legend(loc='lower right')
    f2_ax[i].set_yscale('log')
    #f2_ax[i].set_xscale('log')
    f2_ax[i].set_ylim(10**-9,10**-3)
    #f2_ax[i].set_xlim(0,0.11)
    f2_ax[i].grid(True)
    i += 1

f2_ax3.set_ylabel(r'$\Delta E / E_{\rm t=0}$', fontsize=16)
f2_ax5.set_xlabel(r'$\eta$', fontsize=16)
f2_ax6.set_xlabel(r'$\eta$', fontsize=16)




plt.show()
