import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np
import cmath

def plot_montecarlo_variance(ax, filename, color, label):
    draw, I = np.loadtxt(filename, usecols=(0,1), skiprows=0, unpack=True)
    ax.plot(draw, I, 'o', color=color, label=label)
    mean = np.mean(I)
    sdev = np.std(I)
    ax.plot([0, 100], [mean, mean], ':', color=color)
    ax.fill_between([0, 100], [mean - sdev, mean - sdev], [mean + sdev, mean + sdev], color=color, alpha=0.2)

def plot_montecarlo_comparison():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'draw')
    ax.set_xlim([0, 100])
    ax.set_ylabel(r'I')
    ax.set_ylim([0.7, 1.0])

    plot_montecarlo_variance(ax, 'output/mcmc_0.txt', 'tab:red', 'mean value')
    plot_montecarlo_variance(ax, 'output/mcmc_1.txt', 'tab:blue', 'importance sampling')

    ax.text(10, 0.97, 'N = 10$^3$', fontsize=22)

    ax.legend(loc='upper right')
    plt.savefig('IntegrationMontecarloComparison.pdf')

def plot_integration_results():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    
    ax.set_xlabel(r'log N')
    ax.set_xlim([1, 6])
    ax.set_ylabel(r'log($\tilde I$ - $I$)')
    ax.set_ylim([-6, 0])
    
    log2N, I = np.loadtxt('build/montecarlo.txt', usecols=(0,1), unpack=True)
    N = np.power(2., log2N)
    
    ax.plot(np.log10(N), np.log10(abs(I - 0.5) + 1e-40), 'o')
    
    x = np.linspace(1,6,100)
    ax.plot(x, -0.5 - 0.5 * x)
    plt.savefig('IntegrationMontecarloResults.pdf')

if __name__ == "__main__":
    plot_montecarlo_comparison()
    #plot_integration_results()
