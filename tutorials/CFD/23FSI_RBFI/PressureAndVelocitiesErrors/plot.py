import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.figure(figsize=(10, 4))

import sys
sys.path.append("..")
from errL2U128_mat import errL2U as errL2U128
from errL2U1212_mat import errL2U as errL2U1212
from errL2U1015_mat import errL2U as errL2U1015
from errL2U1510_mat import errL2U as errL2U1510
from errL2U2010_mat import errL2U as errL2U2010
from errL2U1010_mat import errL2U as errL2U1010
from errL2U2525_mat import errL2U as errL2U2525

from errL2P128_mat import errL2P as errL2P128
from errL2P1212_mat import errL2P as errL2P1212
from errL2P1015_mat import errL2P as errL2P1015
from errL2P1510_mat import errL2P as errL2P1510
from errL2P2010_mat import errL2P as errL2P2010
from errL2P1010_mat import errL2P as errL2P1010
from errL2P2525_mat import errL2P as errL2P2525


times = np.arange(0, 30.1, 0.1)
#plt.subplot(1, 2, 1)
#plt.title(r"\textbf{Velocityerrors Vs Time}", fontsize=10)
plt.title(r"\textbf{Pressure errors Vs Time}", fontsize=12)
plt.xlabel(r"\textbf{Time (second)}")
plt.ylabel(r"$\epsilon_p [\%]$", fontsize=11)
#plt.yscale('log')


plt.plot(times, 100*errL2P128,  'o-',  label=r"${N}_u =12, {N}_p =8$", linewidth=2.0, markevery=0.1)
plt.plot(times, 100*errL2P2525, '*-',  label=r"${N}_u =25, {N}_p =25$", linewidth=2.0, markevery=0.1)
plt.plot(times, 100*errL2P1212, 'x-',  label=r"${N}_u =12, {N}_p =12$", linewidth=2.0, markevery=0.1)
plt.plot(times, 100*errL2P1510, 'd-',  label=r"${N}_u =15, {N}_p =10$", linewidth=2.0,  markevery=0.1)
plt.plot(times, 100*errL2P1010, 's-',  label=r"  ${N}_u =10, {N}_p =10$", linewidth=2.0, markevery=0.1)
plt.plot(times, 100*errL2P1015, 'v-',  label=r" ${N}_u =10, {N}_p =15$", linewidth=2.0, markevery=0.1)
plt.plot(times, 100*errL2P2010, '^-',  label=r"${N}_u =20, {N}_p =10$", linewidth=2.0, markevery=0.1)


'''
plt.plot(times, 100*errL2U128,  'o-',  label=r"${N}_u =12, {N}_p =8$",  linewidth=2.0,  ms=4, markevery=0.5)
plt.plot(times, 100*errL2U1212, 's-',  label=r"${N}_u =12, {N}_p =12$",  linewidth=2.0,  ms=4, markevery=0.5)
plt.plot(times, 100*errL2U1510, '*-',  label=r"${N}_u =15, {N}_p =10$",  linewidth=2.0,  ms=4, markevery=0.5)
plt.plot(times, 100*errL2U1010, '^-',  label=r"  ${N}_u =10, {N}_p =10$",  linewidth=2.0,  ms=4, markevery=0.5)
plt.plot(times, 100*errL2U1015, 'P-',  label=r" ${N}_u =10, {N}_p =15$",  linewidth=2.0,  ms=4, markevery=0.5)
plt.plot(times, 100*errL2U2010, 'd-', label=r"${N}_u =20, {N}_p =10$",  linewidth=2.0,  ms=4, markevery=0.5)
plt.plot(times, 100*errL2U2525, 'x-', label=r"${N}_u =25, {N}_p =25$",  linewidth=2.0,  ms=4, markevery=0.5)
'''
plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)
plt.legend()
plt.savefig('PressureErrors.pdf')
#plt.savefig('VelocityErrors.pdf')

plt.show()

'''
plt.plot(times, 100*errL2P128,  'o-',  label=r"${N}_u =12, {N}_p =8$")
plt.plot(times, 100*errL2P1212, 's-',  label=r"${N}_u =12, {N}_p =12$")
plt.plot(times, 100*errL2P1510, '*-',  label=r"${N}_u =15, {N}_p =10$")
plt.plot(times, 100*errL2P1010, '^-',  label=r"  ${N}_u =10, {N}_p =10$")
plt.plot(times, 100*errL2P1015, 'P-',  label=r" ${N}_u =10, {N}_p =15$")
plt.plot(times, 100*errL2P2010, 'd-', label=r"${N}_u =20, {N}_p =10$")
'''
