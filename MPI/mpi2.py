from mpi4py import MPI
from scipy import *

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()
print "Hello, World! I am process %d of %d on %s\n" %(rank, size, name),

if rank==0:
  data = arange(10,dtype=int)
else:
  data = empty(10,dtype=int)

comm.Bcast(data,root=0)
for i in range(10):
  assert data[i]==i