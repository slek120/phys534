import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

Gf = numpy.loadtxt('Gf_HL.out',skiprows=1)
Sig= numpy.loadtxt('Sig_HL.out',skiprows=1)

plt.plot(Gf[:,0], Gf[:,1], 'o-', label='real')
plt.plot(Gf[:,0], Gf[:,2], 'o-', label='imag')
plt.xlabel('$i\omega_n$')
plt.ylabel('$G(i\omega_n)$')
plt.legend(loc='best')
plt.xlim(0,4)
plt.savefig('Gf_HL.eps')
plt.show()

plt.clf()

plt.plot(Sig[:,0], Sig[:,1], 'o-', label='real')
plt.plot(Sig[:,0], Sig[:,2], 'o-', label='imag')
plt.xlabel('$i\omega_n$')
plt.ylabel('$\Sigma(i\omega_n)$')
plt.legend(loc='best')
plt.xlim(0,4)
plt.savefig('Sig_HL.eps')
plt.show()