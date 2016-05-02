import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scipy import *
from scipy import integrate
import time

import os
import subprocess

def now():
  return time.strftime("%H:%M:%S\t")

def Compute_Delta(niom, T, mu, Sig, iom, DOS, ommesh):
  # Local Green function
  Gloc = zeros(niom, dtype=complex)
  for i in range(niom):
    Re = mu - ommesh - Sig[i].real
    Im = iom[i] - Sig[i].imag
    denom = 1/(Re**2 + Im**2)
    ReInt = DOS*Re*denom
    ImInt = DOS*Im*denom
    Gloc[i] = integrate.simps(ReInt, ommesh) - 1j*integrate.simps(ImInt, ommesh)

  Delta = 1j*iom+mu-Sig-1./Gloc
  with open('Delta.inp', 'w') as f:
    for i in range(niom):
      f.write('%.8f\t%.8f\t%.8f\n'%(iom[i], Delta[i].real, Delta[i].imag))
  return Gloc

def Create_PARAMS(params):
  with open('PARAMS', 'w') as f:
    for p in params:
      f.write(p + '\t' + str(params[p][0]) + '\t' + params[p][1] + '\n')

def fermi(eps):
  if eps > 100:
    return 0.
  elif eps < -100:
    return 1.
  else:
    return 1./(1.+exp(eps))

def Compute_Nlatt(niom, T, mu, Sig, iom, DOS, ommesh):
  SigInf = Sig[-1].real
  
  Nlatt = 0
  for i in range(niom):
    DMFTW = DOS * (1./(1j*iom[i] + mu - ommesh - Sig[i]) - 1./(1j*iom[i] + mu - ommesh - SigInf))
    Nlatt += 2 * T * integrate.simps(DMFTW.real,ommesh)
  
  Dist = array([fermi((eps+SigInf-mu)/T) for eps in ommesh])
  Nlatt += integrate.simps(DOS * Dist, ommesh)
  return Nlatt

def Compute_Ekin(niom, T, mu, Sig, iom, DOS, ommesh):
  SigInf = Sig[-1].real
  
  Ekin = 0.
  for i in range(niom):
    DMFTW = ommesh * DOS * (1./(1j*iom[i] + mu - ommesh - Sig[i]) - 1./(1j*iom[i] + mu - ommesh - SigInf))
    Ekin += 2 * T * integrate.simps(DMFTW.real,ommesh)
  
  Dist = array([fermi((eps+SigInf-mu)/T) for eps in ommesh])
  Ekin += integrate.simps(ommesh * DOS * Dist, ommesh)
  return Ekin

def Compute_Epot(niom, T, Sig, iom, Gloc, Nlatt):
  SigInf = Sig[-1].real
  Sig1   = Sig[-1].imag * iom[-1]

  Epot = 0.
  for i in range(niom):
    Epot += ((Sig[i]-SigInf)*Gloc[i]).real - Sig1/iom[i]**2
  Epot *= T

  Epot += Sig1/(8.*T) + 0.5*SigInf*Nlatt
  return Epot

if __name__=='__main__':
  params = {
    'Delta' : ['Delta.inp',    '# Input hybridizations'],
    'cix'   : ['impurity.cix', '# Input atomic state'],
    'U'     : [ 6.,            '# Coulomb interaction'],
    'mu'    : [ 3.,            '# Chemical potential'],
    'beta'  : [100.,           '# Inverse tempurature'],
    'M'     : [2e6,            '# Number of Monte Carlo steps'],
    'nom'   : [80,             '# Number of sampled points'],
    'sampleGtau' : [1000,      '# Sampling step of Gtau'],
    'warmup'     : [250000,    '# Number of warm-up steps']
  }

  niom = 6000
  Niter= 40
  U    = params['U'][0]
  T    = 1./params['beta'][0]
  mu   = params['mu'][0]
  Sig  = zeros(niom, dtype=complex)
  mix_Sig = 0.4

  # Matsubara Frequencies
  iom = pi*T*(2*arange(niom)+1)
  # Load DOS
  DOSfile = loadtxt('2D_SL_DOS')
  # 1st column as energies
  ommesh = DOSfile[:,0]
  # 2nd column as DOS
  DOS = DOSfile[:,1]
  # Normalize
  DOS = DOS / integrate.simps(DOS, ommesh)

  out = open('OUTPUT', 'w')

  # Start DMFT loop
  for it in range(Niter):
    Sig_old = array(Sig) # Copy Sig
    if os.path.exists('Sig.out'):
      with open('Sig.out', 'r') as f:
        lines = f.readlines()[1:]
        for i,l in enumerate(lines):
          Sig[i] = float(l.split()[1]) + 1j*float(l.split()[2])
    else:
      print('Sig.out not found. Using default U/2.')
      Sig[:] = 0.5*U

    if it>0:
      Sig = mix_Sig*Sig + (1-mix_Sig)*Sig_old

    Gloc  = Compute_Delta(niom, T, mu, Sig, iom, DOS, ommesh)
    Nlatt = Compute_Nlatt(niom, T, mu, Sig, iom, DOS, ommesh)
    Ekin  =  Compute_Ekin(niom, T, mu, Sig, iom, DOS, ommesh)
    Epot  =  Compute_Epot(niom, T, Sig, iom, Gloc, Nlatt)
    out.write(now()+'Loop #%d\tmu=%.2f\tN=%.2f\tEtot=%.4f\tEkin=%.4f\tEpot=%.4f\n'\
      %(it+1, mu, 2*Nlatt, 2*(Ekin+Epot), 2*Ekin, 2*Epot))
    out.flush()

    Create_PARAMS(params)

    with open('para_com.dat', 'r') as f:
      para_com = str(f.readline())[:-1]

    out.write(now()+'Starting CTQMC\n')
    out.flush()

    cmd = para_com + ' ./ctqmc > ctqmc.log 2> ctqmc.err'
    subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    cmd = 'cp Gf.out Gf.out.' + str(it)
    print os.popen(cmd).read()
    cmd = 'cp Sig.out Sig.out.' + str(it)
    print os.popen(cmd).read()
    cmd = 'cp Delta.inp Delta.inp.' + str(it)
    print os.popen(cmd).read()

    out.write(now()+"Finish CTQMC\t"+"\n")
    out.flush()

  out.close()
