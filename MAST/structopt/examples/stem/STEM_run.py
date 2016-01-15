from MAST.structopt import Optimizer
from MAST.structopt.tools.StemCalc import ConvStem
from MAST.structopt import inp_out
from mpi4py import MPI
import numpy as np
import math
import time
import scipy.interpolate

rank = MPI.COMM_WORLD.Get_rank()
if rank==0:
    # experimental device information 
    aber=[[0,0],[0,0],[22.56,-20.1],[22.08,-7.5],[0.1198,0],[0.9018,-170.1],[0.04964,20.9],[28.43,-120.6],[11.84,153.8],[8.456,76.1],[0.622,0],[2.811,-125.5]]
    autostemparameters={'Electron energy': 200,'Spherical aberration': 1.4,'Defocus':0,'Aperture semiangle': 24.5,'Source size': 0.882,'Slice size': 25.0,'Pixels':976,'Chromatic aberration Coefficient':1.4,'Delta E':0.73,'aber':aber,'Scale Factor':0.00570113}

    # based on precalculated PSF("PSF.txt") and structure information("STEM_ref"), calculate the reference STEM image
    A = ConvStem(parameters=autostemparameters,calc_exp=False)
    atoms_ref=inp_out.read_xyz('STEM_ref',0)

    nk = autostemparameters['Pixels']
    A.psf = np.empty([nk,nk],dtype=float)
    fileobj = open('PSF.txt', 'r')
    lines = fileobj.readlines()
    for x in range(0,nk):
       A.psf[x] = lines[x].split()
    fileobj.close()

    STEM_ref = A.get_image(A.psf, atoms_ref, autostemparameters['Slice size'], autostemparameters['Pixels'])
    autostemparameters['Exp_Image'] = STEM_ref
    
    #setup parameters
    autostemparameters['Grid_sim2exp'] = 1
    autostemparameters['Pixelshift'] = False

    parameters={'structure':'Cluster',
     'optimizer_type':'GA',
     'atomlist':[('Au',1,196.9665,-3.6)],'natoms': 309,
     'nindiv': 20, 'maxgen': 4000, 'reqrep': 1200,
     'restart': False,
     'generate_flag': 'sphere',
     'tolerance': 0.001,
     'r_ab': 2.5, 'size': 23.0,
     'cxpb': 0.8, 'mutpb': 0.2, 'cx_scheme': 'cxtp', #rotct_rand', 
     'mutation_options': ['lattice_alteration','lattice_alteration_group','rotation','stem_mutation','zp_rotation'],
     'selection_scheme': 'tournament2', 'natural_selection_scheme': 'fuss', 'tournsize': 3, 'fusslimit': 10,
     'convergence_scheme': 'gen_rep_min', 'fingerprinting': False, 'predator': 'fitpred', 'demin': 0.05,
     'pair_style': 'eam', 'pot_file': 'Au_u3.eam',
     'filename': 'Output',
     'lammps_keep_files': True, 'lammps_min': '1e-8 1e-8 5000 10000', 'lammps_min_style': 'cg\nmin_modify line quadratic', 'lammps_thermo_steps': 1000,
     'genealogy': True, 'allenergyfile': True, 'best_inds_list': True, 'postprocessing': True,
     'parallel': True,
     'energy_cutoff_factor': 4.0,
     'stem_keep_files': False, 'stem_coeff': 10.0
    }
    parameters['stem_parameters']=autostemparameters
    parameters['fitness_scheme']='stem_cost'
else:
    parameters=None
parameters = MPI.COMM_WORLD.bcast(parameters,root=0)
C = Optimizer(parameters)
C.run()
