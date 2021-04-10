import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import numpy as np
import cmath

def set_axis(ax):
    ax.set_xlabel(r'x')
    ax.set_xlim([0, 1])
    ax.set_ylabel(r'f(x)')
    ax.set_ylim([0, 1])
    
fig = plt.figure(figsize=(10.5,8))
ax = fig.add_subplot(111)
set_axis(ax)
 
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
