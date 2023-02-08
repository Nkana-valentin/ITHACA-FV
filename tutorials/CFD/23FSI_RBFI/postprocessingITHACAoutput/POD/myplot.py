import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True




datap = np.loadtxt("Eigenvalues_p", skiprows=2, dtype=float)
datau = np.loadtxt("Eigenvalues_U", skiprows=2, dtype=float)
datapd = np.loadtxt("Eigenvalues_pointDisplacement", skiprows=2, dtype=float)

#print(data[0])
plt.yscale('log')
plt.xlabel(r'$\textbf{N}$', fontsize=16)
plt.ylabel(r'$1-\displaystyle\sum_{n=1}^N' r'\lambda_n$', fontsize=16)
plt.plot(datap, 'k--', label='p')
plt.plot(datau, 'b--', label='U')
plt.plot(datapd, 'r--', label='pointDisplacement')
plt.legend()
plt.show()
