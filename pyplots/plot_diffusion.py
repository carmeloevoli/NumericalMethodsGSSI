import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np

def plot_filename(ax, filename, color, linestyle):
    t, f = np.loadtxt(filename, skiprows=0, usecols=(0,1), unpack=True)
    ax.plot(t, f, lw=2.5, color=color, linestyle=linestyle)

def plot_analytical_solution(ax):
    x = np.linspace(-1, 1, 2000)
    y = np.cos(0.5 * np.pi * x)
    ax.plot(x, y, '--', lw=1, color='k')

def plot_numerical_solutions(ax, filename, linestyle, color):
#    plot_filename(ax, filename + '_0.txt', 'k', ':')
#    plot_filename(ax, filename + '_2.txt', 'tab:orange', linestyle)
#    plot_filename(ax, filename + '_4.txt', 'tab:olive', linestyle)
#    plot_filename(ax, filename + '_6.txt', 'tab:green', linestyle)
#    plot_filename(ax, filename + '_8.txt', 'tab:blue', linestyle)
    plot_filename(ax, filename + '_1.txt', color, linestyle)

def make_plot():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    ax.set_xlabel('x')
    ax.set_ylabel('f(x)')
    ax.set_xlim([-1., 1.])
    ax.set_ylim([0.0, 1.1])
    return fig, ax
    
def plot_CN():
    fig, ax = make_plot()
    plot_numerical_solutions(ax, 'output/diffusion_10_8', '-', 'tab:red')
    plot_numerical_solutions(ax, 'output/diffusion_20_8', '-', 'tab:orange')
#    plot_numerical_solutions(ax, 'output/CN_0.6_9', '-', 'tab:orange')
#    plot_numerical_solutions(ax, 'output/CN_0.6_10', '-', 'tab:green')
#    plot_numerical_solutions(ax, 'output/CN_0.6_11', '-', 'tab:blue')
    plot_analytical_solution(ax)
    ax.set_title('Crank-Nicolson')
#    ax.text(0.75, 0.90, r'$2^8$', color='tab:red')
#    ax.text(0.75, 0.80, r'$2^9$', color='tab:orange')
#    ax.text(0.75, 0.70, r'$2^{10}$', color='tab:green')
#    ax.text(0.75, 0.60, r'$2^{11}$', color='tab:blue')
    plt.savefig('DiffusionCN.pdf')
    
if __name__== "__main__":
    plot_CN()
