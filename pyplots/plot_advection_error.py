import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np

def set_axis(ax):
    ax.set_xlabel('log N')
    ax.set_xlim([1, 5])
    ax.set_ylabel(r'log(<E>)')
    ax.set_ylim([-4,  1])

def plot_n_slope(ax, log2N, I_err, color, label):
    N = np.power(2., log2N)
    logN = np.log10(N)
    logErr = np.log10(abs(I_err))
    ax.plot(logN, logErr, 'o', color=color, label=label)

    m = np.polyfit(logN[-4:], logErr[-4:], 1)
    print(logN[-4:])
    print (f"slope = {m[0]:6.3f}")
    print (f"norm = {m[1]:6.3f}")
    x = np.linspace(0, 10, 1000)
    #ax.plot(x, m[1] + m[0] * x, color=color, linestyle=':')

fig = plt.figure(figsize=(10.5,8))
ax = fig.add_subplot(111)
set_axis(ax)
    
filename = 'advection_convergence_test.txt'
log2N, I1_err, I2_err, I3_err, I4_err = np.loadtxt(filename, skiprows=0, usecols=(0,1,2,3,4), unpack=True)

plot_n_slope(ax, log2N, I1_err, 'tab:red', 'Upwind')
plot_n_slope(ax, log2N, I2_err, 'tab:orange', 'Lax-Friedrichs')
plot_n_slope(ax, log2N, I3_err, 'tab:green', 'Lax-Wendroff')
#plot_n_slope(ax, log2N, I4_err, 'tab:purple')

ax.legend()

plt.savefig('AdvectionError.pdf')
