from scipy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



if __name__=='__main__':
   EF=-0.4494
   dos=loadtxt('graphene.dos')
   plt.plot(dos.T[0]-EF, dos.T[1], '-', label='dos')
   plt.xlabel('Energy [eV]')
   plt.ylabel('DOS')
   plt.legend(loc='best')
   plt.xlim(-5,5)
   plt.savefig('dos.eps')
