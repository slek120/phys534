from mpi4py import MPI
from scipy import *
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

print("Hello, World! I am process %d of %d on %s."%(rank, size, name))

data = ones(10,dtype=int)*(rank+1)**2

data_sum=zeros(10,dtype=int)

print 'data=', data, ' at processor', rank
comm.Reduce(data, data_sum, op=MPI.SUM, root=0)

if rank==0:
  print 'summed data=', data_sum
