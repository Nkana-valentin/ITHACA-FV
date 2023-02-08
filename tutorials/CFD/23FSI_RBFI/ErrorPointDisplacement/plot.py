import numpy as np
import matplotlib.pyplot as plt
#plt.rcParams['text.usetex'] = True
plt.rcParams.update({ "text.usetex": True, "font.family": "Avant Garde" })
plt.figure(figsize=(15, 5))

import sys
sys.path.append("..")
from ErrorPdispl_mat1015 import ErrorPdispl as ErrorPdispl1015
from ErrorPdispl_mat2020 import ErrorPdispl as ErrorPdispl2020
from ErrorPdispl_mat1010 import ErrorPdispl as ErrorPdispl1010
from ErrorPdispl_mat3030 import ErrorPdispl as ErrorPdispl3030
from ErrorPdispl_mat88 import ErrorPdispl   as ErrorPdispl88
from ErrorPdispl_mat816 import ErrorPdispl  as ErrorPdispl816


times = np.arange(0, 30.1, 0.1)
#plt.subplot(1, 2, 1)
#plt.title(r"\textbf{Velocityerrors Vs Time}", fontsize=12)
plt.title(r"\textbf{PointDisplacement  errors Vs Time}", fontsize=15)
plt.xlabel(r"\textbf{Time [s]}", fontsize=15)
plt.ylabel(r"$\epsilon_{pointDisplacement} [\%]$", fontsize=20)
plt.tick_params(axis="x", labelsize=18)
plt.tick_params(axis="y", labelsize=18)
#plt.ylabel(r"$\epsilon_u [\%]$", fontsize=15)

plt.plot(times, 100*ErrorPdispl88,  'b*',  label=r"${N}_u =8, {N}_p =8$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='b', markersize=8, markerfacecolor='white')
plt.plot(times, 100*ErrorPdispl816, 'ro',  label=r"${N}_u =8, {N}_p =16$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='r', markersize=8, markerfacecolor='white')
plt.plot(times, 100*ErrorPdispl1015, 'kd',  label=r"${N}_u =10, {N}_p =15$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='k', markersize=8, markerfacecolor='white')
plt.plot(times, 100*ErrorPdispl1010, 'gv',  label=r"  ${N}_u =10, {N}_p =10$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='g', markersize=8, markerfacecolor='white')
plt.plot(times, 100*ErrorPdispl2020, 'ms',  label=r" ${N}_u =20, {N}_p =20$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='m', markersize=8, markerfacecolor='white')
plt.plot(times, 100*ErrorPdispl3030, 'c^', label=r"${N}_u =30, {N}_p =30$",  linewidth=2.0,  linestyle='dashed',  markevery=0.1, mec='c', markersize=8, markerfacecolor='white')


plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)
plt.legend(fontsize=15)
plt.savefig('ErrorDisplacement.pdf')
#plt.savefig('VelocityErrors.pdf')

plt.show()
