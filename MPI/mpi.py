from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()
print "Hello, World! I am process %d of %d on %s\n" %(rank, size, name),

if rank==0:
  data = {'a':7,'b':3.14}
  # req = comm.isend(data,dest=1,tag=11)
  # req.wait()
else:
  data = None
  # req = comm.irecv(source=0,tag=11)
  # data = req.wait()

data=comm.bcast(data,root=0)
print "data a is %d, b is %f, in proc %d" %(data['a'], data['b'], rank)