import random
import numpy
from MAST.structopt.tools.find_defects import find_defects

def move_la(indiv, Optimizer):
    """Move function to move atoms in structure by lattice constant.  Intended for use in Defect optimization.
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
    #Move structure by lattice constant
    Optimizer.output.write('Lattice Constant Move Mutation performed on individual\n')
    if Optimizer.structure=='Defect':
        if Optimizer.isolate_mutation:
            indc,indb,vacant,swaps,stro = find_defects(indiv[0],Optimizer.solidbulk,0)
            ind = indc.copy()
        else:
            ind = indiv[0]
        la = numpy.maximum.reduce(numpy.maximum.reduce(ind.get_cell())/Optimizer.solidcell)
        ax = [[1,0,0],[0,1,0],[0,0,1],[-1,0,0],[0,-1,0],[0,0,-1]]
        selax = [la*i for i in random.choice(ax)]
        ind.translate(selax)
        indiv[0]=ind.copy()
        if Optimizer.isolate_mutation:
            indiv[0].extend(indb)
    else:
        Optimizer.output.write('WARNING: Move Mutation performed on non-Defect structure\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='LM'
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv
    