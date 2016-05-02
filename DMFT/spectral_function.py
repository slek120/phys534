from scipy import *
from scipy.integrate import simps

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def Compute_Aw(om, DOS, Sig, ommesh, delta=0.1j):
  Aw = zeros(len(om))
  for i in range(len(om)):
    DMFTW = DOS / (om[i] - Sig[i] - ommesh + delta)
    Aw[i] = simps(-1.*DMFTW.imag/pi, ommesh)
  return Aw

if __name__ == '__main__':
  # Load DOS
  DOSfile = loadtxt('2D_SL_DOS')
  # 1st column as energies
  ommesh = DOSfile[:,0]
  # 2nd column as DOS
  DOS = DOSfile[:,1]
  # Normalize
  DOS = DOS / simps(DOS, ommesh)

  # Load Sig
  Sigfile = loadtxt('Sig.out')
  # 1st column as frequencies
  om = Sigfile[:,0]
  # 2nd, 3rd column as self energy
  Sig = Sigfile[:,1] + 1j * Sigfile[:,2]

  Aw = Compute_Aw(om, DOS, Sig, ommesh)
  plt.plot(om, Aw)
  plt.legend(['U=6'], loc='best')
  plt.ylabel('$A(\omega)$')
  plt.xlabel('$\omega$')
  plt.savefig('Aw_U6.eps')