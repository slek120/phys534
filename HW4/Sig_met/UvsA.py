from scipy import *
from scipy.integrate import simps
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def Compute_Aw(om, DOS, Sig, ommesh, delta=0.1j):
  DMFTW = DOS / (om - Sig - ommesh + delta)
  return simps(-1.*DMFTW.imag/pi, ommesh)

# Load DOS
DOSfile = loadtxt('2D_SL_DOS')
# 1st column as energies
ommesh = DOSfile[:,0]
# 2nd column as DOS
DOS = DOSfile[:,1]
# Normalize
DOS = DOS / simps(DOS, ommesh)

for f in os.listdir(os.getcwd()):
  if f.startswith('Sig.out'):
    # Load Sig
    Sigfile = loadtxt(f)
    # 1st column as frequencies
    om = Sigfile[0,0]
    # 2nd, 3rd column as self energy
    Sig = Sigfile[0,1] + 1j * Sigfile[0,2]
    Aw = Compute_Aw(om, DOS, Sig, ommesh, 0.01j)

    with open('UvsA.dat','a') as dat:
      dat.write('%s\t%f\n'%(f,Aw))