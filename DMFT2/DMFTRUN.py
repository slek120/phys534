from scipy import *
from scipy import integrate
import os

def Compute_Delta(niom,T,mu,Sig):

   iom=(arange(niom)*2+1)*pi*T # Matsubara frequency
   DOS_f=loadtxt('2D_SL_DOS')
   ommesh=DOS_f.T[0]; DOS=DOS_f.T[1] 
   DOS=DOS/integrate.simps(DOS,ommesh) # Normalization
   Gloc=zeros(niom,dtype=complex)
   for i in range(niom):
      denom_re=mu-ommesh-Sig[i].real
      denom_im=iom[i]-Sig[i].imag
      denom2=denom_re**2+denom_im**2
      Gloc_re=DOS*denom_re/denom2
      Gloc_im=DOS*denom_im/denom2
      Gloc[i]=integrate.simps(Gloc_re,ommesh)-1j*integrate.simps(Gloc_im,ommesh)

   Delta=1j*iom+mu-Sig-1./Gloc
   fi=open('Delta.inp','w')
   for i in range(niom):
      print >>fi, '%.8f %.8f %.8f' %(iom[i], Delta[i].real, Delta[i].imag)
   fi.close()


if __name__=='__main__':

   niom=6000; T=0.01; mu=0.0
   Sig=zeros(niom,dtype=complex) 
   Compute_Delta(niom,T,mu,Sig) 

   if os.path.exists('Sig.out'):
      with open('Sig.out', 'r') as f:
         lines = f.readlines()[1:]
         for i,l in enumerate(lines):
            print(i,l)
   else:
      print('Sig.out was not found')

