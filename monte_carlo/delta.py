from scipy import *
from scipy import integrate
import matplotlib.pyplot as plt

mu=0.
W=2.0
rom=201
nom=6000
T=0.01

ommesh=linspace(-W/2,W/2,rom)
iom=(2*arange(nom)+1)*pi*T
DOS=ones(rom)*(1./W)

Delta=zeros(nom,dtype=complex)
for i in range(nom):
  denom=iom[i]**2+(mu-ommesh)**2
  Del_re=DOS*(mu-ommesh)/denom
  Del_im=DOS*iom[i]/denom
  Delta[i]=integrate.simps(Del_re,ommesh)-1j*integrate.simps(Del_im,ommesh)

with open('Delta.inp', 'w') as f:
  for i in range(nom):
    f.write('%.8f %.8f %.8f\n'%(iom[i], Delta[i].real,Delta[i].imag))

plt.plot(iom, Delta.real, 'o-', label='real')
plt.plot(iom, Delta.imag, 'o-', label='imag')
plt.xlabel('$i\omega_n$')
plt.ylabel('$\Delta(i\omega_n$')
plt.legend(loc='best')
plt.xlim(0,4)
plt.savefig('Delta.eps')
plt.show()