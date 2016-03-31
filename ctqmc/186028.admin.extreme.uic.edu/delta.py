import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scipy import *
from scipy import integrate
import TB

W=2.0
rom=201
nom=6000
T=0.01

t  = 1.
tp = 0.
e0 = 0.
mu = 0.
nk = 300
ndos = 101
SL = TB.square_lattice(nk,mu,hopping=(t,tp,e0))
SL.Create_irr_klist()
SL.Fourier_Transform()
SL.Compute_eigvals()
ommesh, DOS = TB.Compute_DOS(ndos,SL.eigvals,SL.wklist)

iom=(2*arange(nom)+1)*pi*T

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
plt.ylabel('$\Delta(i\omega_n)$')
plt.legend(loc='best')
plt.xlim(0,4)
plt.savefig('Delta.eps')
plt.show()