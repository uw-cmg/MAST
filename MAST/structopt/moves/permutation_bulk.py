import random
from MAST.structopt.tools.find_defects import find_defects

def permutation_bulk(indiv, Optimizer):
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
    if Optimizer.structure != 'Defect':
        Optimizer.output.write('WARNING: Crystal extended permutation attempted on non-Defect system. Skipping Mutation\n')
    else:
        a1=indiv.bulki[random.randint(0,indiv.bulki.get_number_of_atoms()-1)]
        opts=[inds for inds in indiv.bulki if inds.symbol !=a1.symbol]
        try:
            a2=opts[random.randint(0,len(opts)-1)]
        except ValueError:
            a2=indiv.bulki[random.randint(0,indiv.bulki.get_number_of_atoms()-1)]
            Optimizer.output.write('WARNING: Permutation Mutation attempted on single atom structure system\n')
        indiv.bulki[a1.index].symbol, indiv.bulki[a2.index].symbol=a2.symbol, a1.symbol
        Optimizer.output.write('Bulk Permutation Mutation performed on individuals bulk\n')
        Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
        Optimizer.output.write(repr(a1.symbol) + ' swapped with ' + repr(a2.symbol)+'\n')
        Optimizer.output.write(repr(indiv.bulki)+'\n')
        muttype='PB'
        if indiv.energy==0:
            indiv.history_index=indiv.history_index+'m'+muttype
        else:
            indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv