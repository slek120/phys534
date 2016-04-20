import os
import subprocess

def Create_scf(i):
  with open('graphene.scf%d.in'%i, 'w') as f:
    a = 2.4623 + (i-5)*0.01
    s = "&CONTROL\n\
      calculation = 'scf',\n\
     restart_mode = 'from_scratch',\n\
       pseudo_dir = './',\n\
           outdir = './',\n\
           prefix = 'graphene',\n\
/\n\
&SYSTEM\n\
            ibrav = 4,\n\
                a = %.4f,\n\
                c = 10\n\
              nat = 2,\n\
             ntyp = 1,\n\
      occupations = 'smearing',\n\
         smearing = 'methfessel-paxton',\n\
          degauss = 0.02,\n\
          ecutwfc = 60,\n\
          ecutrho = 720,\n\
             nbnd = 8,\n\
/\n\
&ELECTRONS\n\
         conv_thr = 1.0d-10,\n\
      mixing_mode = 'plain',\n\
      mixing_beta = 0.7,\n\
  diagonalization = 'cg',\n\
/\n\
ATOMIC_SPECIES\n\
C 12.0107 C.pbe-rrkjus.UPF\n\
ATOMIC_POSITIONS {crystal}\n\
C  0.333333333  0.666666666  0.500000000\n\
C  0.666666666  0.333333333  0.500000000\n\
K_POINTS {automatic}\n\
 42 42 1   0 0 0\n\
"%a
    f.write(s)

if __name__=='__main__':
  with open('para_com.dat', 'r') as f:
    para_com = str(f.readline())[:-1]

  for i in range(10):
    Create_scf(i)
    cmd = para_com + ' pw.x -i graphene.scf%d.in > espresso.log 2> espresso.err'%(it)
    subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    cmd = 'cp out out.' + str(it)
    print os.popen(cmd).read()
