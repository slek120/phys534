from mpi4py import MPI
from scipy import *
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

print("Hello, World! I am process %d of %d on %s."%(rank, size, name))

sendbuf = zeros(100, dtype=int)+rank
recvbuf = None

if rank == 0:
  recvbuf=empty([size, 100], dtype=int)

comm.Gather(sendbuf, recvbuf, root=0)

if rank == 0:
  for i in range(size):
    assert allclose(recvbuf[i,:], i)
else:
  assert recvbuf == None

print('data=', recvbuf, 'rank=', rank)