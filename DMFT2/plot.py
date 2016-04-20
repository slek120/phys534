from scipy import *
from scipy import integrate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



if __name__=='__main__':

   file=['Sig.out.36','Sig.out.37','Sig.out.38','Sig.out.39']
   lt_list=['o-','s-','^-','>-']
   cl_list=['red','blue','purple','black']
   for i,fi in enumerate(file):
      Sig=loadtxt(fi)
      plt.plot(Sig.T[0], Sig.T[1], lt_list[i], color=cl_list[i], label=fi)
   plt.xlabel('i$\omega_n$')
   plt.ylabel('Re$\Sigma(i\omega_n)$')
   plt.legend(loc='best')
   plt.xlim(0,4)
   plt.ylim(1.3,1.7)
   plt.savefig('ReSigma.eps')
   plt.close()
   for i,fi in enumerate(file):
      Sig=loadtxt(fi)
      plt.plot(Sig.T[0], Sig.T[2], lt_list[i], color=cl_list[i], label=fi)
   plt.xlabel('i$\omega_n$')
   plt.ylabel('Im$\Sigma(i\omega_n)$')
   plt.legend(loc='best')
   plt.xlim(0,4)
   plt.ylim(-0.4,0)
   plt.savefig('ImSigma.eps')
