from mpi4py import MPI
from scipy import *
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

print("Hello, World! I am process %d of %d on %s."%(rank, size, name))

sendbuf=None
if rank==0:
  sendbuf=empty([size,100],dtype=int)
  sendbuf.T[:,:]=range(size)

recvbuf=empty(100,dtype=int)
comm.Scatter(sendbuf,recvbuf,root=0) # Use Captial letter

assert allclose(recvbuf, rank)
print('data=', recvbuf, 'rank=', rank)