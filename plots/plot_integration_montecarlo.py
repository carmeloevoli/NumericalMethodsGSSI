import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np
import cmath

def plot_integration_function():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'x')
    ax.set_xlim([0, 1])
    ax.set_ylabel(r'f(x)')
    ax.set_ylim([0, 1])
 
    x = np.linspace(0, 1, 1000)
    f = 0.5 * np.sin(x * 8.0 * cmath.pi) + 0.5

    ax.plot(x, f, color='tab:red')
    ax.fill_between(x, f, 0, color='tab:red', alpha=0.3)

    N = 1000

    x = np.random.rand(N)
    f = 0.5 * np.sin(x * 8.0 * cmath.pi) + 0.5
    y = np.random.rand(N)

    ax.plot(x[np.where(f > y)], y[np.where(f > y)], 'o', color='tab:red')
    ax.plot(x[np.where(f < y)], y[np.where(f < y)], 'o', color='tab:blue')

    plt.savefig('IntegrationMontecarlo.pdf')

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

plot_integration_results()
