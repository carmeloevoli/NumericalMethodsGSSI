import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np

def set_axis(ax):
    ax.set_xlabel('log$_2$ N')
    ax.set_xlim([0, 30])
    ax.set_ylabel(r'log($\tilde I$ - $I$)')
    ax.set_ylim([-15, 0])
    
fig = plt.figure(figsize=(10.5,8))
ax = fig.add_subplot(111)
set_axis(ax)
    
filename = 'build/trapezium_v0.txt'
N, I, I_err = np.loadtxt(filename, skiprows=0, usecols=(0,1,2), unpack=True)
I_err = abs(I_err)

ax.plot(N, np.log10(I_err), 'o')

m = np.polyfit(N[0:10], np.log10(I_err[0:10]), 1)

print (m)

plt.savefig('Integration.pdf')
