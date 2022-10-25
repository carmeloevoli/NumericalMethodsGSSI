import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np

def plot_filename(ax, filename, color, linestyle):
    t, f = np.loadtxt(filename, skiprows=0, usecols=(0,1), unpack=True)
    ax.plot(t, f, lw=2.5, color=color, linestyle=linestyle)

def gaussian(x, mu, sigma):
    return np.exp(-0.5 * np.power((x - mu) / sigma, 2.0))

def plot_analytical_solutions_gaussian(ax):
    x = np.linspace(0, 1, 2000)
    y = gaussian(x, 0.1, 0.01)
    v = 1.
    ax.plot(x + v * 0.6, y, '--', lw=1, color='k')

def wave(x):
    return np.sin(8.0 * x * np.pi);

def plot_analytical_solutions(ax):
    x = np.linspace(0, 1, 2000)
    y = wave(x)
    ax.plot(x, y, '--', lw=1, color='k')

def plot_numerical_solutions(ax, filename, linestyle, color):
    plot_filename(ax, filename + '_0.txt', 'k', ':')
    plot_filename(ax, filename + '_2.txt', 'tab:orange', linestyle)
    plot_filename(ax, filename + '_4.txt', 'tab:olive', linestyle)
    plot_filename(ax, filename + '_6.txt', 'tab:green', linestyle)
    plot_filename(ax, filename + '_8.txt', 'tab:blue', linestyle)
    plot_filename(ax, filename + '_1.txt', color, linestyle)

def make_plot():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    ax.set_xlabel('x')
    ax.set_ylabel('f(x)')
    ax.set_xlim([0., 1.0])
    ax.set_ylim([-1.1, 1.1])
    return fig, ax

def plot_Upwind():
    fig, ax = make_plot()
    plot_numerical_solutions(ax, 'output/upwind_1_7', '-', 'tab:red')
    plot_numerical_solutions(ax, 'output/upwind_1_8', '--', 'tab:green')
    plot_numerical_solutions(ax, 'output/upwind_1_9', ':', 'tab:blue')
    plot_analytical_solutions(ax)
    ax.set_title('Upwind')
    plt.savefig('AdvectionUpwind.pdf')

def plot_Laxf():
    fig, ax = make_plot()
    plot_numerical_solutions(ax, 'output/laxf_1_7', '-', 'tab:red')
    plot_numerical_solutions(ax, 'output/laxf_1_8', '--', 'tab:green')
    plot_numerical_solutions(ax, 'output/laxf_1_9', ':', 'tab:blue')
    plot_analytical_solutions(ax)
    ax.set_title('Laxf')
    plt.savefig('AdvectionLF.pdf')

def plot_Laxw():
    fig, ax = make_plot()
    plot_numerical_solutions(ax, 'output/laxw_1_7', '-', 'tab:red')
    plot_numerical_solutions(ax, 'output/laxw_1_8', '--', 'tab:green')
    plot_numerical_solutions(ax, 'output/laxw_1_9', ':', 'tab:blue')
    plot_analytical_solutions(ax)
    ax.set_title('Laxw')
    plt.savefig('AdvectionLW.pdf')
    
def plot_CN():
    fig, ax = make_plot()
    ax.set_xlim([0.5,0.8])
    ax.set_ylim([0, 1.1])
    plot_numerical_solutions(ax, 'output/CN_0.6_8', '-', 'tab:red')
    plot_numerical_solutions(ax, 'output/CN_0.6_9', '-', 'tab:orange')
    plot_numerical_solutions(ax, 'output/CN_0.6_10', '-', 'tab:green')
    plot_numerical_solutions(ax, 'output/CN_0.6_11', '-', 'tab:blue')
    plot_analytical_solutions_gaussian(ax)
    ax.set_title('Crank-Nicolson')
    ax.text(0.75, 0.90, r'$2^8$', color='tab:red')
    ax.text(0.75, 0.80, r'$2^9$', color='tab:orange')
    ax.text(0.75, 0.70, r'$2^{10}$', color='tab:green')
    ax.text(0.75, 0.60, r'$2^{11}$', color='tab:blue')
    plt.savefig('AdvectionCN.pdf')
    
def plot_box_solution(ax, filename, color, linestyle = '-'):
    z, f = np.loadtxt(filename + '.txt', skiprows=0, usecols=(0,1), unpack=True)
    ax.plot(z, f, color=color, linestyle=linestyle)
    
def plot_box():
    fig, ax = make_plot()
    ax.set_xlim([0,1])
    ax.set_ylim([0, 1.1])
    plot_box_solution(ax, 'output/box_5_0', 'tab:blue')
    plot_box_solution(ax, 'output/box_5_2', 'tab:green', linestyle='--')
    plot_box_solution(ax, 'output/box_5_4', 'tab:red', linestyle='--')
    plot_box_solution(ax, 'output/box_5_6', 'tab:orange', linestyle='--')
    plot_box_solution(ax, 'output/box_5_8', 'tab:olive', linestyle='--')
    plt.savefig('AdvectionBox4.pdf')

if __name__== "__main__":
#    plot_Upwind()
#    plot_Laxw()
#    plot_Laxf()
#    plot_CN()
    plot_box()
