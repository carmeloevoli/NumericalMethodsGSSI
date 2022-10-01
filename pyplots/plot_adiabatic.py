import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np

def plot_numerical_solution(ax, filename, color, linestyle):
    t, f = np.loadtxt(filename, skiprows=0, usecols=(0,1), unpack=True)
    ax.plot(t, f, lw=2.5, color=color, linestyle=linestyle)

def plot_numerical_solution_log(ax, filename, color, linestyle):
    t, f = np.loadtxt(filename, skiprows=0, usecols=(0,1), unpack=True)
    ax.plot(np.exp(t), f, lw=2.5, color=color, linestyle=linestyle)
    
def plot_analytical_solution(ax):
    p, f = np.loadtxt('output/solution_0.txt', skiprows=0, usecols=(0,1), unpack=True)
    ax.plot(p, f, '--', lw=1, color='k')

def make_plot():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    ax.set_xlabel('p')
    ax.set_xscale('log')
    ax.set_ylabel('f(p)')
    #ax.set_xlim([-1., 1.])
    #ax.set_ylim([0.0, 1.1])
    return fig, ax
    
def plot_CN():
    fig, ax = make_plot()
    plot_analytical_solution(ax)
    plot_numerical_solution(ax, 'output/CN_0.txt', 'tab:red', '-')
    plot_numerical_solution(ax, 'output/CN_2.txt', 'tab:green', '-')
    plot_numerical_solution(ax, 'output/CN_4.txt', 'tab:orange', '-')
    plot_numerical_solution(ax, 'output/CN_6.txt', 'tab:blue', '-')

    plot_numerical_solution(ax, 'output/upwind_0.txt', 'tab:red', ':')
    plot_numerical_solution(ax, 'output/upwind_2.txt', 'tab:green', ':')
    plot_numerical_solution(ax, 'output/upwind_4.txt', 'tab:orange', ':')
    plot_numerical_solution(ax, 'output/upwind_6.txt', 'tab:blue', ':')

    plot_numerical_solution_log(ax, 'output/CNLog_0.txt', 'tab:red', '--')
    plot_numerical_solution_log(ax, 'output/CNLog_2.txt', 'tab:green', '--')
    plot_numerical_solution_log(ax, 'output/CNLog_4.txt', 'tab:orange', '--')
    plot_numerical_solution_log(ax, 'output/CNLog_6.txt', 'tab:blue', '--')

#    ax.set_title('Crank-Nicolson')
    plt.savefig('adiabatic.pdf')
    
if __name__== "__main__":
    plot_CN()
