import numpy as np
import matplotlib.pyplot as plt
#plt.rcParams['text.usetex'] = True
plt.rcParams.update({ "text.usetex": True, "font.family": "Avant Garde" })
plt.figure(figsize=(18, 8))

import sys
sys.path.append("..")
from pdcoeffrbf_mat1015 import pdcoeffrbf as pdcoeffrbf1015
from pdcoeffrbf_mat2020 import pdcoeffrbf as pdcoeffrbf2020
from pdcoeffrbf_mat1010 import pdcoeffrbf as pdcoeffrbf1010
from pdcoeffrbf_mat3030 import pdcoeffrbf as pdcoeffrbf3030
from pdcoeffrbf_mat88   import pdcoeffrbf as pdcoeffrbf88
from pdcoeffrbf_mat816  import pdcoeffrbf as pdcoeffrbf816

from coeffpd_mat import coeffpd


times = np.arange(0, 30.1, 0.1)
plt.title(r"\textbf{PointDisplacement coefficients Vs Time}", fontsize=15)
plt.xlabel(r"\textbf{Time [s]}", fontsize=15)
plt.ylabel(r"\textit{PointDisplacement coefficients}", fontsize=15)

plt.tick_params(axis="x", labelsize=20)
plt.tick_params(axis="y", labelsize=20)

plt.plot(times, pdcoeffrbf88,  'b*',  label=r"${N}_u =8, {N}_p =8$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='b', markersize=8, markerfacecolor='white')
plt.plot(times, pdcoeffrbf816, 'ro',  label=r"${N}_u =8, {N}_p =16$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='r', markersize=8, markerfacecolor='white')
plt.plot(times, pdcoeffrbf1015, 'kd',  label=r"${N}_u =10, {N}_p =15$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='k', markersize=8, markerfacecolor='white')
plt.plot(times, pdcoeffrbf1010, 'gv',  label=r"  ${N}_u =10, {N}_p =10$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='g', markersize=8, markerfacecolor='white')
plt.plot(times, pdcoeffrbf2020, 'ms',  label=r" ${N}_u =20, {N}_p =20$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='m', markersize=8, markerfacecolor='white')
plt.plot(times, pdcoeffrbf3030, 'c^', label=r"${N}_u =30, {N}_p =30$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='c', markersize=8, markerfacecolor='white')

plt.plot(times, coeffpd[0], 'y-', label="offline")


plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)
plt.legend(fontsize=20)
plt.savefig('pointDisplacementcoeff.pdf')

plt.show()
