from scipy import *
from scipy import integrate
import os,subprocess,time

def now():
   return time.strftime("at %H:%M:%S")

def Compute_Delta(niom,iom,DOS,ommesh,T,mu,Sig):

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

def Create_PARAMS(params):
   f=open('PARAMS','w')
   print >>f, "# Input file for CTQMC"
   for p in params:
      print >>f, p, params[p][0], params[p][1]
   f.close()

def Fermi(epsilon):
   if epsilon>100: return 0.
   elif epsilon<-100: return 1.
   else: return 1./(1.+exp(epsilon))

def Compute_Nlatt(Sig,niom_cut,DOS,iom,mu,ommesh,T):
   Sigoo=Sig[-1].real
   Nlatt=0.
   for i in range(niom_cut):
      DMFTW=DOS*(1./(1j*iom[i]+mu-ommesh-Sig[i])-1./(1j*iom[i]+mu-ommesh-Sigoo))
      Nlatt+=2*T*integrate.simps(DMFTW.real,ommesh)
   Nlatt+=integrate.simps(DOS*array([Fermi((eps+Sigoo-mu)/T) for eps in ommesh]),ommesh)
   return Nlatt 	

if __name__=='__main__':

   params={"Delta" : ['Delta.inp', '# Input hybridization'],
           "cix" :['impurity.cix', '# Input atomic state'],
           "U" :  [ 12.0, '# Coulomb interaction'],
           "mu" : [ 3.0, '# Chemical potential'],
           "beta":[100.0, '# Inverse temperature'],
           "M"   :[2e6, '# Number of Monte Carlo steps'],
           "nom" :[ 80, '# Number of sampled points'],
         "warmup":[250000, '# Number of warm-up steps']
           }
   niom=6000; Niter=40; mix_Sig=0.8 
   T=1./params["beta"][0]; mu=params["mu"][0]; U=params["U"][0]
   Sig=zeros(niom,dtype=complex) 

   iom=(arange(niom)*2+1)*pi*T # Matsubara frequency
   DOS_f=loadtxt('2D_SL_DOS')
   ommesh=DOS_f.T[0]; DOS=DOS_f.T[1] 
   DOS=DOS/integrate.simps(DOS,ommesh) # Normalization

   main_out=open("OUTPUT","w")
   iter_out=open("ITER_INFO","w")
   iter_out.write("it  mu  Nlatt"+'\n')
   iter_out.flush()
   ############  Starting DMFT Loop ###############
   for it in range(Niter):
      main_out.write("Starting DMFT loop "+str(it+1)+" at "+now()+"\n")
      main_out.flush()
    
      Sig_old=array(Sig)
      if os.path.exists('Sig.out'):
         lines=open('Sig.out','r').readlines()[1:]
         for i,line in enumerate(lines):
            Sig[i]=float(line.split()[1])+1j*float(line.split()[2])
      else: Sig[:]=U/2.
        
      if it>0: Sig=Sig_old+mix_Sig*(Sig-Sig_old)

      Compute_Delta(niom,iom,DOS,ommesh,T,mu,Sig) 

      Nlatt=Compute_Nlatt(Sig,niom,DOS,iom,mu,ommesh,T)
   
      iter_out.write("%2d %.2f %.2f \n" %(it+1, mu, 2*Nlatt))
      iter_out.flush()
      Create_PARAMS(params)

      fipa=open('para_com.dat','r')
      para_com=str(fipa.readline())[:-1]
      fipa.close()

      main_out.write("Starting CTQMC at "+now()+"\n") 
      main_out.flush()

      cmd=para_com+" ./ctqmc > ctqmc.log 2> ctqmc.err" 
      subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

      cmd = "cp Gf.out Gf.out."+str(it)
      print os.popen(cmd).read() # copying Gf 
      cmd = "cp Sig.out Sig.out."+str(it)
      print os.popen(cmd).read() # copying Sig 
      cmd = "cp Delta.inp Delta.inp."+str(it) 
      print os.popen(cmd).read() # copying Delta





