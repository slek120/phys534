from scipy import *
import time

def elapsed(start):
    return "%.2f sec"%(time.time()-start)

def ComputeMCSteps(IS, NStep, NWarm, NMeas):
  MC_E = IS.E
  MC_M = IS.M
  L = IS.L
  L2 = IS.L2
  Navg = 0.
  Eavg = 0.
  Mavg = 0.
  E2avg = 0.
  M2avg = 0.

  for istep in range(NStep):
    # Pick random site (i,j)
    Site = int(rand()*L2)
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
  E = Eavg/Navg/L2
  M = Mavg/Navg/L2
  Cv = (E2avg/Navg-(Eavg/Navg)**2)/L2/IS.T**2
  chi = (M2avg/Navg-(Mavg/Navg)**2)/L2/IS.T
  return E, M, Cv, chi

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

  def UpdateFP(self, T):
    self.FP = zeros(9, dtype=float)
    self.FP[0] = 1.
    self.FP[2] = 1.
    self.FP[4] = 1.
    self.FP[6] = exp(-4./T)
    self.FP[8] = exp(-8./T)
    self.T = T


def main():
  L = 100
  Nstep,Nwarm,Nmeasure = (1000000, 1000, 100)

  IS = Ising(L)
  IS.ComputeE()
  IS.ComputeM()
  print('E: %f, M: %f'%(IS.E, IS.M))

  IS.UpdateFP(4.)
  for i in range(10):
    print('E: %f, M: %f, Cv: %f, chi: %f'%ComputeMCSteps(IS,Nstep,Nwarm,Nmeasure))


  # for temp in [0.01,0.1,1.,10.]:
  #   start = time.time()
  #   IS.UpdateFP(temp)
  #   print('E: %f, M: %f, Cv: %f, chi: %f'%ComputeMCSteps(IS,Nstep,Nwarm,Nmeasure))
  #   print(elapsed(start))

if __name__ == '__main__':
  main()