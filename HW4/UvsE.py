from scipy import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

f = loadtxt('UvsE_ins.dat')
U = f[:,0]
Ekin = f[:,1]
Epot = f[:,2]
Etot = f[:,3]

plt.plot(U, Ekin, label='$E_{kin}^{Insulator}$')
plt.plot(U, Epot, label='$E_{pot}^{Insulator}$')
plt.plot(U, Etot, label='$E_{tot}^{Insulator}$')

f = loadtxt('UvsE_met.dat')
U = f[:,0]
Ekin = f[:,1]
Epot = f[:,2]
Etot = f[:,3]

plt.plot(U, Ekin, label='$E_{kin}^{Metal}$')
plt.plot(U, Epot, label='$E_{pot}^{Metal}$')
plt.plot(U, Etot, label='$E_{tot}^{Metal}$')

plt.legend(loc='best')
plt.ylabel('$E$')
plt.xlabel('$U$')
plt.savefig('UvsE.png')