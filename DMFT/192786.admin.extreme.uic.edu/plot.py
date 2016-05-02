from scipy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

UvsEfile = loadtxt('UvsE.dat')
U = UvsEfile[:,0]
E = UvsEfile[:,3]
plt.plot(U,E)
plt.savefig('UvsE.eps')
