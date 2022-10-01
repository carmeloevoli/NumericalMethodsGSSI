import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np

def plot_numerical_solution(ax, filename, color, linestyle):
    z, f = np.loadtxt(filename, skiprows=0, usecols=(0,1), unpack=True)
    ax.plot(z, f, lw=2.5, color=color, linestyle=linestyle)
    
def plot_analytical_solution(ax):
    z, f = np.loadtxt('output/solution_0.txt', skiprows=0, usecols=(0,1), unpack=True)
    ax.plot(z, f, '--', lw=1, color='k')

def make_plot():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    ax.set_xlabel('z')
    ax.set_ylabel('f(z)')
    #ax.set_xlim([-1., 1.])
    #ax.set_ylim([0.0, 1.1])
    return fig, ax
    
def plot_CN():
    fig, ax = make_plot()
    plot_numerical_solution(ax, 'output/galaxy_0.txt', 'tab:red', '-')
    plot_numerical_solution(ax, 'output/galaxy_3.txt', 'tab:green', '-')
    plot_numerical_solution(ax, 'output/galaxy_6.txt', 'tab:orange', '-')
    plot_numerical_solution(ax, 'output/galaxy_9.txt', 'tab:olive', '-')
    plot_analytical_solution(ax)

#    ax.set_title('Crank-Nicolson')
    plt.savefig('galaxy.pdf')
    
if __name__== "__main__":
    plot_CN()
