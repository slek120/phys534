import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mpi4py import MPI
from scipy import *
import time

def elapsed(start):
    return "%.2f sec"%(time.time()-start)

def ComputeMCSteps(IS, NStep, NWarm, NMeas):
  MC_E = IS.E
  MC_M = IS.M
  L = IS.L
  Navg = 0.
  Eavg = 0.
  Mavg = 0.
  E2avg = 0.
  M2avg = 0.

  for istep in range(NStep):
    # Pick random site (i,j)
    Site = int(rand()*L*L)
    i,j = Site//L, Site%L
    # Calculate S_i.S_j for n.n.
    SI = IS.spins[i,j]
    SJ = IS.spins[(i-1)%L,j]\
       + IS.spins[(i+1)%L,j]\
       + IS.spins[i,(j-1)%L]\
       + IS.spins[i,(j+1)%L]
    # Get probability to flip spin
    P = IS.FP[4+SI*SJ]
    if rand() < P:
      # Flip spin at site (i,j)
      IS.spins[i,j] = -SI
      MC_E += 2*SI*SJ
      MC_M -= 2*SI
    if istep > NWarm and istep%NMeas == 0:
      Navg += 1.
      Eavg += MC_E
      Mavg += MC_M
      E2avg += MC_E**2
      M2avg += MC_M**2
  E = Eavg/Navg/IS.L2
  M = Mavg/Navg/IS.L2
  E2 = E2avg/Navg/IS.L2
  M2 = M2avg/Navg/IS.L2
  return E, M, E2, M2

class Ising:
  def __init__(self,L):
    """Create Ising spin class with LxL sites"""
    self.L = L
    self.L2 = L**2
    # self.spins = array(list(map(int,sign(rand(L**2)-0.5)))).reshape(L,L)
    self.spins = ones(L**2).reshape(L,L)

  def ComputeE(self):
    self.E = 0.
    for i in range(self.L):
      for j in range(self.L):
        SI = self.spins[i,j]
        SJ = self.spins[(i-1)%self.L,j]\
           + self.spins[(i+1)%self.L,j]\
           + self.spins[i,(j-1)%self.L]\
           + self.spins[i,(j+1)%self.L]
        self.E -= SI*SJ
    self.E *= 0.5

  def ComputeM(self):
    self.M = sum(self.spins)

  def UpdateT(self, T):
    self.FP = zeros(9, dtype=float)
    self.FP[0] = 1.
    self.FP[2] = 1.
    self.FP[4] = 1.
    self.FP[6] = exp(-4./T)
    self.FP[8] = exp(-8./T)
    self.T = T


def main():
  comm = MPI.COMM_WORLD
  size = comm.Get_size()
  rank = comm.Get_rank()

  L = 100
  Nstep,Nwarm,Nmeasure = (1000000, 1000, 100)
  IS = Ising(L)

  Tlist = linspace(4., 0.5, 60)
  Elist = []
  Mlist = []
  Cvlist = []
  chilist = []

  t1 = time.time()
  for T in Tlist:
    t2 = time.time()
    IS.UpdateT(T)
    IS.ComputeE()
    IS.ComputeM()
    E,M,E2,M2 = ComputeMCSteps(IS,Nstep,Nwarm,Nmeasure)
    E = comm.allreduce(E, op=MPI.SUM)
    M = comm.allreduce(abs(M), op=MPI.SUM)
    E2 = comm.allreduce(abs(E2), op=MPI.SUM)
    M2 = comm.allreduce(abs(M2), op=MPI.SUM)
    Eavg = E/size
    Mavg = M/size
    E2avg= E2/size
    M2avg= M2/size
    Cv = (E2avg-Eavg**2)/T**2
    chi= (M2avg-Mavg**2)/T

    Elist.append(Eavg)
    Mlist.append(Mavg)
    Cvlist.append(Cv)
    chilist.append(chi)
    print('T: %f, E: %f, M: %f, Cv: %f, chi: %f'%(T, Eavg, Mavg, Cv, chi))
    print(elapsed(t2))
  print(elapsed(t1))

  if rank==0:
    plt.plot(Tlist, Elist, label='$E(T)$')
    plt.plot(Tlist, Mlist, label='$M(T)$')
    plt.plot(Tlist, Cvlist, label='$C_v(T)$')
    plt.xlabel('$T$')
    plt.legend(loc='best')
    plt.savefig('EMCV.eps')
    plt.close()
    plt.plot(Tlist, chilist, label='$\chi(T)$')
    plt.xlabel('T')
    plt.legend(loc='best')
    plt.savefig('chi.eps')

if __name__ == '__main__':
  main()