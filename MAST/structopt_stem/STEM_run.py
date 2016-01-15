from MAST.structopt_stem import Optimizer
from MAST.structopt_stem.tools.StemCalc import ConvStem
from MAST.structopt_stem import inp_out
from mpi4py import MPI
import numpy as np
import math
import time
import scipy.interpolate

rank = MPI.COMM_WORLD.Get_rank()
if rank==0:

    import ast
    fp = open('structoptinput.txt','r')
    data = fp.readlines()
    for i in range(len(data)):
        data[i] = data[i].strip()
    data=' '.join(data)
    data='{'+data+'}'
    parameters = ast.literal_eval(data)
    parameters['fitness_scheme']='stem_cost'
    del parameters['calc_method']
    del parameters['vaspcalc']

    aber=[[0,0],[0,0],[22.56,-20.1],[22.08,-7.5],[0.1198,0],[0.9018,-170.1],[0.04964,20.9],[28.43,-120.6],[11.84,153.8],[8.456,76.1],[0.622,0],[2.811,-125.5]]
    autostemparameters={'Electron energy': 200,'Spherical aberration': 1.4,'Defocus':0,'Aperture semiangle': 24.5,'Source size': 0.882,'Slice size': 25.0,'Pixels':976,'Chromatic aberration Coefficient':1.4,'Delta E':0.73,'aber':aber,'Scale Factor':0.00570113}

    # calculate ref STEM image
    '''
    fp = open('../archive_input_options.txt')
    data = fp.readlines()
    for line in range(len(data)):
        if "'psf':" in data[line]:
            parameters['psf'] = ast.literal_eval(data[line].strip())['psf']
        if "'stem_ref':" in data[line]:
            parameters['stem_ref'] = ast.literal_eval(data[line].strip())['stem_ref']
    '''
    A = ConvStem(parameters=autostemparameters,calc_exp=False)
    try:
        atoms_ref=inp_out.read_xyz(parameters['stem_ref'],0)
    except KeyError:
        atoms_ref=inp_out.read_xyz('/home/usitguest/USIT/dropbox_app/STEM_ref',0)

    nk = autostemparameters['Pixels']
    A.psf = np.empty([nk,nk],dtype=float)
    try:
        fileobj = open(parameters['psf'], 'r')
    except KeyError:
        fileobj = open('/home/usitguest/USIT/dropbox_app/PSF.txt', 'r')
    lines = fileobj.readlines()
    for x in range(0,nk):
       A.psf[x] = lines[x].split()
    fileobj.close()

    STEM_ref = A.get_image(A.psf, atoms_ref, autostemparameters['Slice size'], autostemparameters['Pixels'])
    
    parameters['stem_parameters']=autostemparameters
    autostemparameters['Exp_Image'] = STEM_ref
    autostemparameters['Grid_sim2exp'] = 1
    autostemparameters['Pixelshift'] = False
else:
    parameters=None
parameters = MPI.COMM_WORLD.bcast(parameters,root=0)
C = Optimizer(parameters)
C.run()
