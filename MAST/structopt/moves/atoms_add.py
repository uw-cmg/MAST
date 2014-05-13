import random
import numpy
from ase import Atom, Atoms

def atoms_add(indiv,Optimizer):
    """Move function to randomly add atoms or groups of atoms to a structure
    Input:
        indiv = individual class object to be altered
        Optimizer = Optimizer class object with needed parameters
    Output:
        indiv = Altered Individual class object
    """
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False
    cell=indiv[0].get_cell()
    nato=indiv[0].get_number_of_atoms()
    xlen=cell[0][0]-cell[0][1]*cell[1][0]/cell[1][1]
    ylen=cell[1][1]-cell[1][0]*cell[0][1]/cell[0][0]
    zlen=cell[2][2]-cell[2][1]*cell[1][2]/cell[1][1]
    cell_min=[0,0,0]
    cell_max=[xlen,ylen,zlen]
    for i in range(random.randint(1,4)):
        atomsym=random.choice(Optimizer.atomlist)[0]
        pos=[random.uniform(cell_min[0],cell_max[0]), random.uniform(cell_min[1],cell_max[1]), random.uniform(cell_min[2],cell_max[2])]
        natom=Atom(atomsym,pos)
        indiv[0].append(natom)
    Optimizer.output.write('Add atom mutation\n')
    if Optimizer.structure=='Crystal':
        a, b, c = cell
        #normalized cell
        an, bn, cn = [numpy.linalg.norm(v) for v in cell]
        cell[0]=a/an
        cell[1]=b/bn
        cell[2]=c/cn
        indiv[0].set_cell(cell*[(indiv[0].get_number_of_atoms())**0.33*Optimizer.r_ab])
    muttype='AA'+repr(indiv[0].get_number_of_atoms()-nato)
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    Optimizer.output.write(repr(indiv[0].get_positions())+'\n')
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv
