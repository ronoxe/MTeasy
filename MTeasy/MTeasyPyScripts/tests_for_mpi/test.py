from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
import hello 
print rank 
