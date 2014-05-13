from MAST.structopt import Optimizer
from mpi4py import MPI

rank = MPI.COMM_WORLD.Get_rank()
parameters = {'structure' : 'Cluster',
            'optimizer_type' : 'GA',
            'atomlist' : [('Au',1,196.9665,0.0)],
            'natoms' : 10,
            'nindiv' : 10,
            'maxgen' : 10,
            'reqrep' : 50,
            'size' : 25.0**0.333*2.28,
            'CX_SCHEME' : 'rotct_rand',
            'FIT_SCHEME' : 'enpafit',
            'SELECTION_SCHEME' : 'Tournament2',
            'NAT_SELECT' : 'FUSS',
            'tournsize' : 3,
            'fusslimit' : 5,
            'convergence_scheme': 'Gen-Rep-Min',
            'predator' : 'Mutation_Dups',
            'mutation_options' : ['Lattice_Alteration','Lattice_Alteration_Group','Rotation_geo','Rotation'],
            'pair_style' : 'eam',
            'pot_file' : 'Au_u3.eam',
            'keep_Lammps_files' : True,
            'LammpsMin' : '1e-8 1e-8 5000 10000',
            'Lmin_style' : 'cg\nmin_modify line quadratic',
            'genealogy' : True,
            'allenergyfile' : True,
            'BestIndsList' : True}
if rank==0:
    print 'Running Serial...'
    parameters['filename'] = 'Test-Serial'
    A = Optimizer(parameters)
    A.run()
    done = True
else:
    done = False
done = MPI.COMM_WORLD.bcast(done,root=0)
print done


print 'Running parallel...'
parameters['filename'] = 'Test-Parallel'
parameters['parallel'] = True
A = Optimizer(parameters)
A.algorithm_parallel()

print 'Running Island_Method'
parameters['filename'] = 'Test-IM'
parameters['algorithm_type'] = 'Island_Method'
parameters['migration_percent'] = 0.1
A = Optimizer(parameters)
#else:
#    A = None
#A = MPI.COMM_WORLD.bcast(A,root=0)
A.island_algorithm()
