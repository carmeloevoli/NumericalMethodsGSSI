import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np

def set_axis(ax):
    ax.set_xlabel('t')
    ax.set_xlim([20, 80])
    ax.set_ylabel(r'$\left\| f(x) - f_a(x) \right\|_2$')
    ax.set_yscale('log')
    ax.set_ylim([1e-7, 1e-2])

fig = plt.figure(figsize=(10.5,8))
ax = fig.add_subplot(111)
set_axis(ax)
    
filename = 'diffusion_time_test.txt'
t, rms1, rms2, rms3 = np.loadtxt(filename, skiprows=0, usecols=(0,1,2,3), unpack=True)

ax.plot(t, rms1, 'o', color='tab:blue', label=r'$2^8$')
#ax.plot(t, rms2, 'o', color='tab:green', label=r'$2^9$')
#ax.plot(t, rms3, 'o', color='tab:orange', label=r'$2^{10}$')

ax.legend()

plt.savefig('DiffusionTimeEvolution.pdf')
