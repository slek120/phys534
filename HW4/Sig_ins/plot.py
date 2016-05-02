from scipy import *
import matplotlib.pyplot as plt

f = loadtxt('UvsA.dat')
x = f[:,0]
y = f[:,1]

plt.plot(x,y)
plt.ylabel('$A(E_f)$')
plt.xlabel('$U$')
plt.savefig('UvsA.png')