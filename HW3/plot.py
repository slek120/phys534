from scipy import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def E(x, k, a, E0):
  return 0.5*k*(x-a)*(x-a)+E0

dat = loadtxt('total_E.dat')
x   = dat[:,0]
y   = dat[:,1]
k   = 5.
a0  = 2.46
E0  = -22.79119210

popt, pcov = curve_fit(E, x, y, p0=(k, a0, E0))

x0  = linspace(x[0],x[-1],20)
y0  = E(x0, *popt)
lab = '$k=%.4f$\n$a_0=%.4f$\n$E_0=%.4f$'%(popt[0], popt[1], popt[2])

plt.plot(x,y,label='Data')
plt.plot(x0,y0,label=lab)
plt.xlabel('$a$')
plt.ylabel('$E$')
plt.legend(loc='best')
plt.savefig('total_E.png')

print(popt)
