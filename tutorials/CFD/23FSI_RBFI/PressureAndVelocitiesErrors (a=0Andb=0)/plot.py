import numpy as np
import matplotlib.pyplot as plt
#plt.rcParams['text.usetex'] = True
plt.rcParams.update({ "text.usetex": True, "font.family": "Avant Garde" })
plt.figure(figsize=(16, 8))

import sys
sys.path.append("..")
from errL2U88_mat import errL2U as errL2U88
from errL2U816_mat import errL2U as errL2U816
from errL2U1015_mat import errL2U as errL2U1015
from errL2U2020_mat import errL2U as errL2U2020
from errL2U1010_mat import errL2U as errL2U1010
from errL2U3030_mat import errL2U as errL2U3030
#from errL2U2525_mat import errL2U as errL2U2525

from errL2P88_mat import errL2P as errL2P88
from errL2P816_mat import errL2P as errL2P816
from errL2P1015_mat import errL2P as errL2P1015
from errL2P1010_mat import errL2P as errL2P1010
from errL2P2020_mat import errL2P as errL2P2020
from errL2P3030_mat import errL2P as errL2P3030
#from errL2P2525_mat import errL2P as errL2P2525


times = np.arange(0, 30.1, 0.1)
#plt.title(r"\textbf{Velocity errors Vs Time}", fontsize=15)
plt.title(r"\textbf{Pressure errors Vs Time}", fontsize=15)
plt.xlabel(r"\textbf{Time [s]}", fontsize=15)
plt.tick_params(axis="x", labelsize=20)
plt.ylabel(r"$\epsilon_p [\%]$", fontsize=25)
#plt.ylabel(r"$\epsilon_u [\%]$", fontsize=25)
plt.tick_params(axis="y", labelsize=20)

plt.plot(times, 100*errL2P88,  'b*',  label=r"${N}_u =8, {N}_p =8$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='b', markersize=8, markerfacecolor='white')
plt.plot(times, 100*errL2P816, 'ro',  label=r"${N}_u =8, {N}_p =16$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='r', markersize=8, markerfacecolor='white')
plt.plot(times, 100*errL2P1015, 'kd',  label=r"${N}_u =10, {N}_p =15$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='k', markersize=8, markerfacecolor='white')
plt.plot(times, 100*errL2P1010, 'gv',  label=r"  ${N}_u =10, {N}_p =10$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='g', markersize=8, markerfacecolor='white')
plt.plot(times, 100*errL2P2020, 'ms',  label=r" ${N}_u =20, {N}_p =20$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='m', markersize=8, markerfacecolor='white')
plt.plot(times, 100*errL2P3030, 'c^', label=r"${N}_u =30, {N}_p =30$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='c', markersize=8, markerfacecolor='white')
'''
plt.plot(times, 100*errL2U88,  'b*',  label=r"${N}_u =8, {N}_p =8$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='b', markersize=8, markerfacecolor='white')
plt.plot(times, 100*errL2U816, 'ro',  label=r"${N}_u =8, {N}_p =16$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='r', markersize=8, markerfacecolor='white')
plt.plot(times, 100*errL2U1015, 'kd',  label=r"${N}_u =10, {N}_p =15$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='k', markersize=8, markerfacecolor='white')
plt.plot(times, 100*errL2U1010, 'gv',  label=r"  ${N}_u =10, {N}_p =10$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='g', markersize=8, markerfacecolor='white')
plt.plot(times, 100*errL2U2020, 'ms',  label=r" ${N}_u =20, {N}_p =20$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='m', markersize=8, markerfacecolor='white')
plt.plot(times, 100*errL2U3030, 'c^', label=r"${N}_u =30, {N}_p =30$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='c', markersize=8, markerfacecolor='white')
#plt.plot(times, 100*errL2U2525, 'x-', label=r"${N}_u =25, {N}_p =25$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='b', markersize=20, markerfacecolor='white')
'''
plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)
plt.legend(loc='best', fontsize=20)
plt.savefig('PressureErrors.pdf')
#plt.savefig('VelocityErrors.pdf')

plt.show()
