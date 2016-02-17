from scipy import *
from scipy import linalg
import matplotlib.pyplot as plt
import time

def elapsed(start):
    return "%.2f sec"%(time.time()-start)

def Compute_DOS(ndos, eklist, wklist, delta=1e-2):
    tmp = sort(eklist.flatten())
    emin = tmp[0]
    emax = tmp[-1]
    ommesh = linspace(emin,emax,ndos)
    dos = zeros(ndos,dtype=float)
    totnk=len(eklist)
    for i,om in enumerate(ommesh):
        for ik, ekl in enumerate(eklist):
            for ek in ekl:
                dos[i] += wklist[ik]*delta/((om-ek)**2+delta**2)
    dos=dos/pi
    return ommesh, dos


def Compute_Tot_E(mu, eklist, wklist):
    Tot_E=0.
    Tot_N=0.
    for ik,ekl in enumerate(eklist):
        for ek in ekl:
            if ek<=mu:
                Tot_E += wklist[ik]*ek
                Tot_N += wklist[ik]
    return Tot_E, Tot_N


def Create_kpath(KPoints,nk_band):
    def dist(a,b): return sqrt(sum((array(a)-array(b))**2))
    # returns the distance of given a and b points
    KPoints=array(KPoints)
    path_len=[]
    for i in range(len(KPoints)-1):
        path_len.append(dist(KPoints[i+1],KPoints[i]))
    path_nk=map(int,nk_band*array(path_len)/sum(path_len))
    klist=[]; dist_K=[0.]; dist_SK=[0.]
    for i,nkk in enumerate(path_nk):
        for n in range(nkk):
            klist.append(KPoints[i]+(KPoints[i+1]-KPoints[i])*n/nkk)
            if len(klist)>1: dist_K.append(dist_K[-1]+dist(klist[-1],klist[-2]))
        dist_SK.append(dist_SK[-1]+path_len[i])
    # Add the ending point
    klist.append(KPoints[-1])
    dist_K.append(dist_K[-1]+dist(klist[-1],klist[-2]))
    return array(klist), array(dist_K), array(dist_SK)

    
class square_lattice:
    def __init__(self,nk,mu,hopping):
        """ nk: number of k-points in kx axis, mu: chemical potential, hopping: hopping matrix elements """
        self.nk=nk; self.mu=mu
        self.t, self.tp, self.e0 = hopping
        self.Hr={
            (1,0)   : [[-self.t]],
            (0,1)   : [[-self.t]],
            (-1,0)  : [[-self.t]],
            (0,-1)  : [[-self.t]],
            (1,1)   : [[self.tp]],
            (1,-1)  : [[self.tp]],
            (-1,1)  : [[self.tp]],
            (-1,-1) : [[self.tp]],
            (0,0)   : [[self.e0]]
        }

    def Create_irr_klist(self):
        kxlist=linspace(0,2*pi,self.nk+1)[:-1]
        self.klist=[]
        self.wklist=[]
        if self.nk%2==0:
            for i in range(self.nk/2+1):
                for j in range(i+1):
                    self.klist.append([kxlist[i],kxlist[j]])
                    if i==0 and j==0:
                        self.wklist.append(1.)
                    elif i==self.nk/2:
                        if j==i:
                            self.wklist.append(1.)
                        elif j==0:
                            self.wklist.append(2.)
                        else:
                            self.wklist.append(4.)
                    elif j==0 or i==j:
                        self.wklist.append(4.)
                    else:
                        self.wklist.append(8.)
        else:
            for i in range((self.nk+1)/2):
                for j in range(i+1):
                    self.klist.append([kxlist[i],kxlist[j]])
                    if i==0 and j==0:
                        self.wklist.append(1.)
                    elif j==0 or i==j:
                        self.wklist.append(4.)
                    else:
                        self.wklist.append(8.)
        self.totnk=len(self.klist)
        print("Total number of irr. k poins=", self.totnk)
        print("Total sum of multiplicity=", sum(self.wklist))
        self.klist=array(self.klist)
        self.wklist=array(self.wklist)/sum(self.wklist) # Normalize

    def Create_full_klist(self):
        kxlist=linspace(0,2*pi,self.nk+1)[:-1]
        self.klist=[[kx,ky] for kx in kxlist for ky in kxlist] 
        self.klist=array(self.klist)
        self.totnk=len(self.klist) # Total number of k-points
        self.wklist=ones(self.totnk)/self.totnk
        print("Total number of k-points=", self.totnk)

    def Fourier_Transform(self):
        self.ndim=len(self.Hr[(0,0)])
        self.Hk=zeros((self.totnk,self.ndim,self.ndim),dtype=complex)
        for ik,kp in enumerate(self.klist):
            for rvec in self.Hr.keys():
                self.Hk[ik,:,:]+=array(self.Hr[rvec])*exp(1j*dot(kp,rvec))

    def Compute_eigvals(self):
        self.eigvals=[]
        for H in self.Hk:
            self.eigvals.append(linalg.eigvalsh(H))
        self.eigvals=array(self.eigvals)

    def Update_klist(self,klist):
        self.klist=array(klist)
        self.totnk=len(self.klist)

class hexagonal_lattice:
    def __init__(self,nk,mu,hopping):
        """ nk: number of k-points in kx axis, mu: chemical potential, hopping: hopping matrix elements """
        self.nk=nk
        self.mu=mu
        self.t, self.tp, self.e0 = hopping
        self.Hr={
            ( 0.5*sqrt(3), 0.5) : [[self.tp, -self.t], [0, self.tp]],       # a1
            ( 0.5*sqrt(3),-0.5) : [[self.tp, -self.t], [0, self.tp]],       # a2
            (-0.5*sqrt(3),-0.5) : [[self.tp, 0], [-self.t, self.tp]],       #-a1
            (-0.5*sqrt(3), 0.5) : [[self.tp, 0], [-self.t, self.tp]],       #-a2
            ( 0, 1)             : [[self.tp, 0], [0, self.tp]],             # a1-a2
            ( 0,-1)             : [[self.tp, 0], [0, self.tp]],             # a2-a1
            ( 0, 0)             : [[self.e0, -self.t], [-self.t, self.e0]]  # 0
        }

    def Create_irr_klist(self):
        # Initial list of k: range=[0,2*pi)
        kxlist=linspace(0,2*pi,self.nk+1)[:-1]
        self.klist=[] # List of (kx,ky) coordinates
        self.wklist=[] # List of multiplicity of corresponding coordinate
        for i,kx in enumerate(kxlist):
            for j,ky in enumerate(kxlist):
                # Fractional coordinate
                k=[(kx+ky)/sqrt(3),kx-ky]
                # Only add to list if inside Gamma-X-K-Gamma triangle
                if k[0]<=2*pi/sqrt(3) and 0 <= k[1] and k[0] <= 3*k[1]:
                    self.klist.append([(kx+ky)/sqrt(3),kx-ky])
                    # Gamma point
                    if k==[0,0]:
                        self.wklist.append(1.)
                    # Along X-K
                    elif k[0]==2*pi/sqrt(3):
                        # K point
                        if k[0] == 3*k[1]:
                            self.wklist.append(3.)
                        # X point 
                        elif 0 == k[1]:
                            self.wklist.append(3.)
                        # Along X-K
                        else:
                            self.wklist.append(6.)
                    # Along Gamma-K
                    elif k[0] == 3*k[1]:
                        self.wklist.append(6.)
                    # Along Gamma-X
                    elif 0 == k[1]:
                        self.wklist.append(6.)
                    # Inside Gamma X K triangle
                    else:
                        self.wklist.append(12.)
        self.totnk=len(self.klist)
        print("Total number of irr. k-points=", self.totnk)
        self.klist=array(self.klist)
        self.wklist=array(self.wklist)/sum(self.wklist) # Normalize

    def Create_full_klist(self):
        kxlist=linspace(0,2*pi,self.nk+1)[:-1]
        self.klist=[[(kx+ky)/sqrt(3),kx-ky] for kx in kxlist for ky in kxlist] 
        self.klist=array(self.klist)
        self.totnk=len(self.klist) # Total number of k-points
        self.wklist=ones(self.totnk)/self.totnk
        print("Total number of k-points=", self.totnk)

    def Fourier_Transform(self):
        self.ndim=len(self.Hr[(0,0)])
        self.Hk=zeros((self.totnk,self.ndim,self.ndim),dtype=complex)
        for ik,kp in enumerate(self.klist):
            for rvec in self.Hr.keys():
                self.Hk[ik,:,:]+=array(self.Hr[rvec])*exp(1j*dot(kp,rvec))

    def Compute_eigvals(self):
        self.eigvals=[]
        for H in self.Hk:
            self.eigvals.append(linalg.eigvalsh(H))
        self.eigvals=array(self.eigvals)

    def Update_klist(self,klist):
        self.klist=array(klist)
        self.totnk=len(self.klist)


if __name__=='__main__':
    t  = 1.
    tp = 0.
    e0 = 0.
    mu = 0.
    nk = 300
    ndos = 101
    HL = hexagonal_lattice(nk,mu,hopping=(t,tp,e0))

    # Calculate Energy dispersion along gamma-X-M-gamma path
    start = time.time()
    KPoints=[[0,0],[2*pi/sqrt(3),0],[2*pi/sqrt(3),2*pi/3],[0,0]]
    nk_band=500
    klist, dist_K, dist_SK = Create_kpath(KPoints,nk_band)
    HL.Update_klist(klist)
    HL.Fourier_Transform()
    HL.Compute_eigvals()
    print("Dispersion: ",elapsed(start))

    # Plot energy dispersion
    label=['$\Gamma$','X','K','$\Gamma$']
    plt.plot(dist_K,HL.eigvals[:,:])
    plt.grid(True)
    plt.xticks(dist_SK,label)
    plt.xlabel('kpath')
    plt.ylabel('E(k)')
    plt.title('Band Structure')
    plt.show()


    # Calculate DOS
    start = time.time()
    HL.Create_irr_klist()
    # HL.Create_full_klist()
    HL.Fourier_Transform()
    HL.Compute_eigvals()
    ommesh, dos = Compute_DOS(ndos,HL.eigvals,HL.wklist)
    print("DOS: ",elapsed(start))

    # Plot DOS
    plt.plot(ommesh,dos,'o-')
    plt.xlabel('$\omega$')
    plt.ylabel('$\\rho(\omega)$')
    plt.title('DOS')
    plt.show()


    # Calculate E vs n_k
    start = time.time()
    nklist = arange(10,200)
    Elist = []
    for nk in nklist:
        HL=hexagonal_lattice(nk,mu,hopping=(t,tp,e0))
        HL.Create_irr_klist()
        # HL.Create_full_klist()
        HL.Fourier_Transform()
        HL.Compute_eigvals()
        Tot_E, Tot_N = Compute_Tot_E(HL.mu, HL.eigvals, HL.wklist)
        Elist.append(Tot_E)
    print("E vs nk: ",elapsed(start))

    # Plot E vs n_k
    plt.plot(nklist,Elist,'o-')
    plt.xlabel('$n_k$')
    plt.ylabel('$E$')
    plt.title('$E \; vs \; n_k$')
    plt.show()
