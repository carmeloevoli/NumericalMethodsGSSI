import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np

def plot_y_vs_t():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    ax.set_xlabel('t')
    ax.set_ylabel('y')
    ax.set_xlim([0, 50])
    ax.set_ylim([-30, 30])
    filename = 'build/lorenz.txt'
    t, y = np.loadtxt(filename, skiprows=0, usecols=(0,2), unpack=True)
    ax.plot(t, y, lw=2, color='tab:red')
    plt.savefig('LorenzY.pdf')

def plot_z_vs_x():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    ax.set_xlabel('x')
    ax.set_xlim([-20, 20])
    ax.set_ylabel('z')
    ax.set_ylim([0, 50])
    filename = 'build/lorenz.txt'
    x, z = np.loadtxt(filename, skiprows=0, usecols=(1,3), unpack=True)
    ax.plot(x, z, lw=2, color='tab:green')
    plt.savefig('LorenzXZ.pdf')

if __name__== "__main__":
    plot_y_vs_t()
    plot_z_vs_x()
