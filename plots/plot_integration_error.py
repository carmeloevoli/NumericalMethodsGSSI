import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np

def set_axis(ax):
    ax.set_xlabel('log N')
    ax.set_xlim([5, 30])
    ax.set_ylabel(r'log($\tilde I$ - $I$)')
    ax.set_ylim([-15, -1])
    
fig = plt.figure(figsize=(10.5,8))
ax = fig.add_subplot(111)
set_axis(ax)
    
filename = 'build/trapezium.txt'
log2N, I, I_err = np.loadtxt(filename, skiprows=0, usecols=(0,1,2), unpack=True)

N = np.power(2., log2N)
logN = np.log10(N)
logErr = np.log10(abs(I_err))
ax.plot(log2N, logErr, 'o')

m = np.polyfit(logN[0:10], logErr[0:10], 1)

print (f"slope = {m[0]:6.3f}")
print (f"norm = {m[1]:6.3f}")

plt.savefig('Integration.pdf')
