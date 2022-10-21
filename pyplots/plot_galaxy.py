import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np
import matplotlib.gridspec as gridspec
from scipy import interpolate

def plot_analytical_solution(ax, filename):
    z, f = np.loadtxt(filename, skiprows=0, usecols=(0,1), unpack=True)
    ax.plot(z, f, '--', lw=1, color='k')

def plot_numerical_solution(ax, filename, color, label, linestyle='-'):
    z, f = np.loadtxt(filename, skiprows=0, usecols=(0,1), unpack=True)
    ax.plot(z, f, lw=2.5, color=color, label=label, linestyle=linestyle)

def plot_numerical_residual(ax, filename, filename_a, color, label, linestyle='-'):
    z_a, f_a =  np.loadtxt(filename_a, skiprows=0, usecols=(0,1), unpack=True)
    fa_new = interpolate.splrep(z_a, f_a, s=0)
    z, f = np.loadtxt(filename, skiprows=0, usecols=(0,1), unpack=True)
    y = []
    for i in range(len(z)):
        f_temp = interpolate.splev(z[i], fa_new, der=0)
        y.append((f[i] - f_temp) / f_temp)
    ax.plot(z, y, color=color, label=label)

def plot_numerical_derivative(ax, filename, color, label, linestyle='-'):
    z, f = np.loadtxt(filename, skiprows=0, usecols=(0,1), unpack=True)
    dz = z[1] - z[0]
    dfdz = np.gradient(f, dz)
    ax.plot(z, dfdz, lw=2.5, color=color, label=label, linestyle=linestyle)

def plot_analytical_derivative(ax, filename):
    z, f = np.loadtxt(filename, skiprows=0, usecols=(0,1), unpack=True)
    dz = z[1] - z[0]
    dfdz = np.gradient(f, dz)
    ax.plot(z, dfdz, '--', lw=1, color='k')
    
def plot_CN():
    fig = plt.figure(figsize=(10.5,12.))

    gs = gridspec.GridSpec(2, 1, height_ratios=(2, 1))
    gs.update(wspace=0.025, hspace=0.1)
    
    #ax1 = fig.add_subplot(211)
    ax1 = plt.subplot(gs[0])

    ax1.set_xlim([-10., 10.])
    ax1.set_ylim([0., 400.])
    ax1.set_ylabel('f(z) / f_0')
    ax1.set_title(r'$\beta = 5$')
    ax1.set_xticklabels([])
    
    plot_numerical_solution(ax1, 'output/galaxy_7_99.txt', 'tab:green', r'$N = 2^7$')
    plot_numerical_solution(ax1, 'output/galaxy_8_99.txt', 'tab:blue', r'$N = 2^8$')
    plot_numerical_solution(ax1, 'output/galaxy_9_99.txt', 'tab:olive', r'$N = 2^9$')
    plot_numerical_solution(ax1, 'output/galaxy_10_99.txt', 'tab:red', r'$N = 2^{10}$')
    plot_analytical_solution(ax1, 'output/solution_0.txt')
    
    ax1.legend()
    
    #ax2 = fig.add_subplot(212)
    ax2 = plt.subplot(gs[1])

    ax2.set_xlim([-10., 10.])
    ax2.set_ylim([-0.2, 0.2])
    ax2.set_ylabel('residual')
    ax2.set_xlabel('z')

    ax2.plot([-10.,10], [0.,0.], linestyle='--', color='tab:gray', lw=1.5)
    plot_numerical_residual(ax2, 'output/galaxy_7_99.txt', 'output/solution_0.txt', 'tab:green', r'$N = 2^7$')
    plot_numerical_residual(ax2, 'output/galaxy_8_99.txt', 'output/solution_0.txt', 'tab:blue', r'$N = 2^8$')
    plot_numerical_residual(ax2, 'output/galaxy_9_99.txt', 'output/solution_0.txt', 'tab:olive', r'$N = 2^9$')
    plot_numerical_residual(ax2, 'output/galaxy_10_99.txt', 'output/solution_0.txt', 'tab:red', r'$N = 2^{10}$')

    plt.savefig('galaxy.pdf',bbox_inches='tight')
    
def plot_der():
    fig, ax = make_plot()
    ax.set_ylabel('df/dz')
    ax.set_xlim([-1,1])
    ax.set_ylim([-2,2])
    ax.set_title(r'$\beta = 5$')

    plot_numerical_solution(ax, 'output/galaxy_7_99.txt', 'tab:green', r'$N = 2^7$')
    plot_numerical_solution(ax, 'output/galaxy_8_99.txt', 'tab:blue',  r'$N = 2^8$')
    plot_numerical_derivative(ax, 'output/galaxy_9_99.txt', 'tab:olive', r'$N = 2^9$')
    plot_numerical_derivative(ax, 'output/galaxy_10_99.txt', 'tab:red', r'$N = 2^{10}$')
    plot_analytical_derivative(ax)

    ax.legend()
    plt.savefig('galaxy_derivative.pdf')
    
if __name__== "__main__":
    plot_CN()
    #plot_der()
