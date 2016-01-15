from MAST.structopt_stem import Optimizer
from mpi4py import MPI

rank = MPI.COMM_WORLD.Get_rank()
if rank==0:
    print 'Running Serial...'
    A = Optimizer('input-clus.txt')
    A.run()
    done = True
else:
    done = False
done = MPI.COMM_WORLD.bcast(done,root=0)
print done

print 'Running parallel...'
A = Optimizer('input-clus-p.txt')
A.run()

print 'Running Island_Method'
A = Optimizer('input-clus-im.txt')
A.run()
