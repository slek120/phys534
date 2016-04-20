from scipy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



if __name__=='__main__':

   band=loadtxt('graphene.bands.gnu')
   plt.plot(band.T[0], band.T[1], '-', label='band')
   plt.xlabel('k')
   plt.ylabel('Energy')
   plt.legend(loc='best')
   plt.xticks([0,0.6667,1,1.5773],['$\Gamma$','K','X','$\Gamma$'])
   plt.savefig('band.eps')
