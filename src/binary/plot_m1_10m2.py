import sys
from matplotlib import use, rcParams, rc, ticker
use('Agg')
import matplotlib.pyplot as plt
import numpy as np
rcParams.update({'figure.autolayout': True})
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Palatino'], 'size': 8})
rc('figure', figsize=(8,9))


#fig, ax = plt.subplots()
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

ax_1 = ax1.twinx()
ax_2 = ax2.twinx()
ax_3 = ax3.twinx()

files = []
files.append('m1_10m2_pec1_2000t_e0p2')
files.append('m1_10m2_pec2_2000t_e0p2')
files.append('m1_10m2_pec3_2000t_e0p2')
files.append('m1_10m2_pec4_2000t_e0p2')

files.append('m1_10m2_pec1_2000t_e0p9')
files.append('m1_10m2_pec2_2000t_e0p9')
files.append('m1_10m2_pec3_2000t_e0p9')
files.append('m1_10m2_pec4_2000t_e0p9')

files.append('m1_10m2_pec1_2000t_e0p99')
files.append('m1_10m2_pec2_2000t_e0p99')
files.append('m1_10m2_pec3_2000t_e0p99')
files.append('m1_10m2_pec4_2000t_e0p99')

info = [
        ['-' ,'#000000',  r'$PEC$'],
        ['-','#ff0000', r'$P(EC)^2$'],
        ['-','#00ff00', r'$P(EC)^3$'],
        ['-','#0000ff', r'$P(EC)^4$']
       ]
for i, f in enumerate(files):
    data = np.genfromtxt(f)

    time = data[:,0]/2.0
    #semimajor = data[:,1]
    #ecc = data[:,2]
    #distance = data[:,3]
    energy = data[:,4]
    cpu_time = data[:,5]

    mk, cl, lb = info[i%4]

    if "e0p2" in f:
        ax_1.plot(time[::1], cpu_time[::1],  color=cl, linestyle='--', label=lb)
        ax1.plot(time[::1],  energy[::1],    color=cl, linestyle=mk)
    elif "e0p99" in f:
        ax_3.plot(time[::1], cpu_time[::1],  color=cl, linestyle='--', label=lb)
        ax3.plot(time[::1],  energy[::1],    color=cl, linestyle=mk)
    elif "e0p9" in f:
        ax_2.plot(time[::1], cpu_time[::1],  color=cl, linestyle='--', label=lb)
        ax2.plot(time[::1],  energy[::1],    color=cl, linestyle=mk)

ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')

ax_1.legend(loc='lower right', fontsize=10)
ax_2.legend(loc='lower right', fontsize=10)
ax_3.legend(loc='lower right', fontsize=10)

ax1.text(200, 1e-8, r'$e_{0} = 0.2$', color='black', fontsize=12,
        bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=0.8'))
ax1.text(530, 1e-6, r'$n = 2,3,4$', fontsize=14)
ax1.text(530, 3e-5, r'$n = 1$',     fontsize=14)

ax2.text(200, 1e-9, r'$e_{0} = 0.9$', color='black', fontsize=12,
        bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=0.8'))
ax2.text(530, 1e-8, r'$n = 2,3,4$', fontsize=14)
ax2.text(530, 2e-5, r'$n = 1$',     fontsize=14)

ax3.text(200, 1e-5, r'$e_{0} = 0.99$', color='black', fontsize=12,
        bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=0.8'))
ax3.text(430, 2e-8, r'$n = 2,3,4$', fontsize=14)
ax3.text(430, 2e-1, r'$n = 1$',     fontsize=14)

#ax.legend(loc='upper right', fontsize=12)
ax1.set_xlim(0,2001)
ax2.set_xlim(0,2001)
ax3.set_xlim(0,2001)

ax1.set_ylim(10e-16, 10)
ax2.set_ylim(10e-16, 10)
ax3.set_ylim(10e-16, 10)
#ax.set_ylabel(r'$\Delta E$', fontsize=20)
#ax.set_ylabel(r'$\Delta a$', fontsize=20)
ax1.set_ylabel(r'$|\Delta E/E_{0}|$', fontsize=20)
ax2.set_ylabel(r'$|\Delta E/E_{0}|$', fontsize=20)
ax3.set_ylabel(r'$|\Delta E/E_{0}|$', fontsize=20)

ax_1.set_ylabel(r'$T_{\rm CPU} [s]$', fontsize=20)
ax_2.set_ylabel(r'$T_{\rm CPU} [s]$', fontsize=20)
ax_3.set_ylabel(r'$T_{\rm CPU} [s]$', fontsize=20)
#ax2.set_ylabel(r'$\Delta e$', fontsize=20)
#ax3.set_ylabel(r'$\Delta E$', fontsize=20)

ax1.set_title(r'$P(EC)^n$ ($a = 1$, $T = 2\pi$, $q=1/10$)', fontsize=20)
ax3.set_xlabel(r'$T(\Omega_{0})^{-1}$', fontsize=20)
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
fig.savefig('m1_10m2.pdf', format='pdf')
#fig.show()
