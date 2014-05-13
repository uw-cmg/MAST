from MAST.structopt import Optimizer
from ase.lattice.cubic import BodyCenteredCubic as BCC
from ase import Atoms,Atom
from MAST.structopt.io.write_xyz import write_xyz
import random
import numpy
import os
from mpi4py import MPI

def main():
	mincomd = '1e-8 1e-8 5000 10000'
    parameters={'structure':'Defect',
            'optimizer_type': 'GA',
            'atomlist':[('Fe',1,55.845,-4.120933)],
            'natoms': 3,
            'filename': 'Fe3',
            'nindiv': 64,
            'maxgen': 500,
            'reqrep': 100,
            'tolerance':0.001,
            'size': 5.0,
            'sf' : 6.0,
            'supercell': (5,5,5),
            'solidfile': 'Bulk.xyz',
            'solidcell': None,
            'evalsolid': False,
            'cxpb': 0.8,
            'mutpb': 0.2,
            'cx_scheme': 'rotct_rand',
            'fitness_scheme': 'totalenfit',
            'natural_selection_scheme': 'fuss',
            'selection_scheme': 'tournament2',
            'predator':'mutation_dups',
            'convergence_shceme': 'gen_rep_min',
            'indiv_defect_write':True,
            'demin': 0.001,
            'number_of_bests':100,
            'mutation_options':['lattice_alteration_rdrd','lattice_alteration_group','rotation_geo','rotation','zp_rotation','quench']
            'quench_max_temp':800,
            'forcing': 'Concentration',
            'calc_method': 'LAMMPS',
            'pair_style':'eam/fs',
            'pot_file':'Fe_mm.eam.fs',
            'lammps_thermo_steps':100,
            'keep_Lammps_files': True,
            'ase_min': False,
            'lammps_min': mincomd,
            'lammps_min_style': 'cg\nmin_modify line quadratic',
            'genealogy': True,
            'output_format': 'FormationEnergyPerInt',
            'allenergyfile':True,
            'bestindslist': True,
            'finddefects':True,
            'parallel' : True,
            'algorithm_type' : 'lambda+mu',
            'debug':['None']}
    bsfe = BCC('Fe', size=parameters['supercell'])
	rank = MPI.COMM_WORLD.Get_rank()
	if rank==0:
    	write(parameters['SolidFile'],bsfe)
	cell=numpy.maximum.reduce(bsfe.get_cell())
	parameters['SolidCell']=cell
    A = Optimizer(parameters)
    A.algorithm_parallel()
                                
                                
if __name__ == "__main__":
        main()
                   
