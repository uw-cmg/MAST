import numpy
import random
from MAST.structopt.tools.find_defects import find_defects

def lattice_alteration(indiv, Optimizer):
    """Move function to perform Lattice Alteration of atoms
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
    cell_max=numpy.maximum.reduce(indiv[0].get_positions())
    cell_min=numpy.minimum.reduce(indiv[0].get_positions())
    if Optimizer.structure=='Defect':
        if Optimizer.isolate_mutation:
            indc,indb,vacant,swap,stro = find_defects(indiv[0],Optimizer.solidbulk,0)
            positions=indc.get_positions()
        else:
            positions=indiv[0].get_positions()
    else:
        positions=indiv[0].get_positions()
    if len(positions) != 0:
        try:
            natomsmove=random.randint(1,len(positions)/5)
        except ValueError:
            natomsmove=1
        ratmlocnew=[0]*natomsmove
        for i in range(natomsmove):
            ratmloc=random.randint(0,len(positions)-1)
            ratmlocnew[i]=(random.uniform(cell_min[0],cell_max[0]), random.uniform(cell_min[1],cell_max[1]), random.uniform(cell_min[2],cell_max[2]))
            positions[ratmloc]=ratmlocnew[i]
        if Optimizer.structure=='Defect':
            if Optimizer.isolate_mutation:
                indc.set_positions(positions)
                indiv[0]=indc.copy()
                indiv[0].extend(indb)
            else:
                indiv[0].set_positions(positions)
        else:
            indiv[0].set_positions(positions)
    else:
        natomsmove=0
        ratmlocnew=0
    Optimizer.output.write('Lattice Alteration Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write('Number of atoms moved = '+repr(natomsmove)+'\n')
    Optimizer.output.write(repr(ratmlocnew)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='LA'+repr(natomsmove)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv