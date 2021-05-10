import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np

def set_axis(ax):
    ax.set_xlabel('log N')
    ax.set_xlim([1, 5])
    ax.set_ylabel(r'log(e)')
    ax.set_ylim([-4,  1])

def plot_n_slope(ax, log2N, I_err, color, label):
    logN = log2N * np.log10(2.)
    logErr = np.log10(abs(I_err))
    ax.plot(logN, logErr, 'o', color=color, label=label)

    m = np.polyfit(logN[-4:], logErr[-4:], 1)
    print(log2N[-4:])
    print (f"slope = {m[0]:6.3f}")
    print (f"norm = {m[1]:6.3f}")
#    x = np.linspace(0, 10, 1000)
#    ax.plot(x, m[1] + m[0] * x, color=color, linestyle=':')

fig = plt.figure(figsize=(10.5,8))
ax = fig.add_subplot(111)
set_axis(ax)
    
filename = 'advection_convergence_test.txt'
log2N, I1_err, I2_err, I3_err = np.loadtxt(filename, skiprows=0, usecols=(0,1,2,3), unpack=True)

plot_n_slope(ax, log2N, I1_err, 'tab:red', 'Lax-Friedrichs')
plot_n_slope(ax, log2N, I2_err, 'tab:orange', 'Upwind')
plot_n_slope(ax, log2N, I3_err, 'tab:green', 'Lax-Wendroff')

ax.legend()

plt.savefig('AdvectionError.pdf')
