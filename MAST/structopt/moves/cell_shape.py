import random
import os
import copy
import numpy
from MAST.structopt.tools.lammps import LAMMPS

def cell_shape(indiv, Optimizer):
    """Move function to forcefully alter the unit cell shape
    Inputs:
        indiv = Individual class object to be altered
        Optimizer = Optimizer class object with needed parameters
    Outputs:
        indiv = Altered Individual class object
    """
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False
    Optimizer.output.write('Cell Shape change mutation\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    structure=random.choice(Optimizer.cell_shape_options)
    cello=indiv[0].get_cell()
    if structure=='cubic':
        #Set to cubic shape
        an,bn,cn = [numpy.linalg.norm(v) for v in cello]
        a=(an+bn+cn)/3.0
        celln=numpy.array([[a,0,0],[0,a,0],[0,0,a]])
        Optimizer.output.write('Mutating to cubic\n')
        muttype='CSC'
    elif structure=='hexagonal':
        #Set to hexagonal shape
        an,bn,cn = [numpy.linalg.norm(v) for v in cello]
        a=(an+bn)/2.0
        c=cn
        if c<=a:
            c=random.uniform(a+1,10)
        trans=numpy.array([[1,0,0],[-0.5,(3.0**0.5)/2.0,0],[0,0,1]])
        trans[0]=[a*i for i in trans[0]]
        trans[1]=[a*i for i in trans[1]]
        trans[2]=[c*i for i in trans[2]]
        celln=trans
        Optimizer.output.write('Mutating to Hexagonal\n')
        muttype='CSH'
    elif structure=='tetragonal':
        #Set to tetragonal shape
        an,bn,cn = [numpy.linalg.norm(v) for v in cello]
        a=(an+bn)/2.0
        c=cn
        if c==a:
            c=random.uniform(1,10)
        celln=numpy.array([[a,0,0],[0,a,0],[0,0,c]])
        Optimizer.output.write('Mutating to tetragonal\n')
        muttype='CSTe'
    elif structure=='orthorhombic':
        #Set to orthorhombic
        a=random.uniform(2,10)
        b=random.uniform(2,10)
        c=random.uniform(2,10)
        celln=numpy.array([[a,0,0],[0,b,0],[0,0,c]])
        Optimizer.output.write('Mutating to orthorhombic\n')
        muttype='CSO'
    elif structure=='monoclinic':
        #Set to monoclinic
        a,b,c = [numpy.linalg.norm(v) for v in cello]
        if a==b:
            b=random.uniform(1,10)
        trans=numpy.array([(1+random.random())*c, 0, (1+random.random())*c])
        celln=numpy.array([[a,0,0],[0,b,0],[0,0,0]])
        celln[2]=trans
        Optimizer.output.write('Mutating to monoclinic\n')
        muttype='CSM'
    elif structure=='triclinic':
        #Set to triclinic
        a,b,c = [numpy.linalg.norm(v) for v in cello]
        celln=numpy.array([[a,0,0],[(1+random.random())*b,(1+random.random())*b,0],[(1+random.random())*c,0,(1+random.random())*c]])
        Optimizer.output.write('Mutating to triclinic\n')
        muttype='CSTr'
    indiv[0].set_cell(celln)
    #Relax new box shape with Lammps box/reax for cell
    Optimizer.output.write('LAMMPS fix box/relax performed\n')
    cwd1=os.getcwd()
    parameters=copy.deepcopy(Optimizer.calc.parameters)
    try:
        parameters['mass'][len(parameters['mass'])-1] += '\nfix 1 all box/relax iso 0.0 vmax 0.001'
    except KeyError:
        parameters['pair_coeff'][0] += '\nfix 1 all box/relax iso 0.0 vmax 0.001'
    filesL = [ Optimizer.pot_file ]
    if Optimizer.lammps_keep_files:
        calc2 = LAMMPS(parameters=parameters, files=filesL,keep_tmp_files=True, tmp_dir=os.getcwd()+'/'+Optimizer.filename+'/Lammpsrun2/')
    else:
        calc2 = LAMMPS(parameters=parameters, files=filesL)
    atmsdup=indiv[0].copy()
    atmsdup.set_calculator(calc2)
    try:
        OUT=atmsdup.calc.calculate(atmsdup)
        indiv[0].set_cell(OUT['atoms'].get_cell())
        indiv[0].set_positions(OUT['atoms'].get_positions())
        Optimizer.output.write('Energy = '+repr(OUT['thermo'][-1]['pe']/indiv[0].get_number_of_atoms())+'\n')
    except:
        Optimizer.output.write('WARNING: Minimizer quit - Box relaxation unsuccessful\n')
    #pdb.set_trace()
    os.chdir(cwd1)
    calc2.clean()
    Optimizer.output.write(repr(indiv[0].get_cell())+'\n')
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv