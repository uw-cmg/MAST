import random
import numpy
from MAST.structopt.tools.find_defects import find_defects

def permutation(indiv, Optimizer):
    """Move function to perform Permutation of one atom based on atomlist
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
    a1=indiv[0][random.randint(0,indiv[0].get_number_of_atoms()-1)]
    opts=[inds for inds in indiv[0] if inds.symbol !=a1.symbol]
    try:
        a2=opts[random.randint(0,len(opts)-1)]
    except ValueError:
        a2=indiv[0][random.randint(0,indiv[0].get_number_of_atoms()-1)]
        Optimizer.output.write('WARNING: Permutation Mutation attempted on single atom structure system\n')
    indiv[0][a1.index].symbol, indiv[0][a2.index].symbol=a2.symbol, a1.symbol
    Optimizer.output.write('Permutation Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write(repr(a1.symbol) + ' swapped with ' + repr(a2.symbol)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='P'
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv