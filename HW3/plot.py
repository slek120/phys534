from scipy import loadtxt, linspace
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Function to fit
def E(x, k, a, E0):
  return 0.5*k*(x-a)*(x-a)+E0

# Data to fit
dat = loadtxt('total_E.dat')
# 1st column as x
x   = dat[:,0]
# 2nd column as y
y   = dat[:,1]

# Initial guess
k0  = 5.
a0  = 2.46
E0  = -22.79119210

# Get fit and covariance
popt, pcov = curve_fit(E, x, y, p0=(k0, a0, E0))

# Make x and y array from fit parameters
xFit  = linspace(x[0],x[-1],50)
yFit  = E(xFit, *popt)
# String for label
lab = '$k=%.4f$\n$a_0=%.4f$\n$E_0=%.4f$'%(popt[0], popt[1], popt[2])

# Plot
plt.plot(x,y,label='Data')
plt.plot(xFit,yFit,label=lab)
plt.xlabel('$a$')
plt.ylabel('$E$')
plt.legend(loc='best')
plt.savefig('total_E.eps')
plt.savefig('total_E.png')
plt.show()
