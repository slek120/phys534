import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

GfHL = numpy.loadtxt('Gf_HL.out',skiprows=1)
SigHL= numpy.loadtxt('Sig_HL.out',skiprows=1)

GfSL = numpy.loadtxt('Gf_SL.out',skiprows=1)
SigSL= numpy.loadtxt('Sig_SL.out',skiprows=1)

plt.plot(GfHL[:,0], GfHL[:,1], 'o-', label='real')
plt.plot(GfHL[:,0], GfHL[:,2], 'o-', label='imag')
plt.xlabel('$i\omega_n$')
plt.ylabel('$G(i\omega_n)$')
plt.legend(loc='best')
plt.xlim(0,4)
plt.savefig('Gf_HL.png')
plt.show()

plt.clf()

plt.plot(SigHL[:,0], SigHL[:,1], 'o-', label='real')
plt.plot(SigHL[:,0], SigHL[:,2], 'o-', label='imag')
plt.xlabel('$i\omega_n$')
plt.ylabel('$\Sigma(i\omega_n)$')
plt.legend(loc='best')
plt.xlim(0,4)
plt.savefig('Sig_HL.png')
plt.show()

plt.clf()

plt.plot(GfSL[:,0], GfSL[:,1], 'o-', label='real')
plt.plot(GfSL[:,0], GfSL[:,2], 'o-', label='imag')
plt.xlabel('$i\omega_n$')
plt.ylabel('$G(i\omega_n)$')
plt.legend(loc='best')
plt.xlim(0,4)
plt.savefig('Gf_SL.png')
plt.show()

plt.clf()

plt.plot(SigSL[:,0], SigSL[:,1], 'o-', label='real')
plt.plot(SigSL[:,0], SigSL[:,2], 'o-', label='imag')
plt.xlabel('$i\omega_n$')
plt.ylabel('$\Sigma(i\omega_n)$')
plt.legend(loc='best')
plt.xlim(0,4)
plt.savefig('Sig_SL.png')
plt.show()