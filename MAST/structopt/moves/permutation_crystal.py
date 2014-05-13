import random
import numpy
from ase import Atom, Atoms
from MAST.structopt.tools.find_defects import find_defects

def permutation_crystal(indiv, Optimizer):
    """Move function to perform Permutation of one atom based on atomlist with entire crystal in Defect
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
        indexbulk=len(indiv[0])
        solid=Atoms()
        solid.extend(indiv[0])
        solid.extend(indiv.bulki)
        solid.set_pbc(True)
        solid.set_cell(indiv.bulko.get_cell())
        a1=solid[random.randint(0,indiv[0].get_number_of_atoms()-1)]
        opts=[inds for inds in solid if inds.symbol != a1.symbol]
        try:
            a2=opts[random.randint(0,len(opts)-1)]
        except ValueError:
            a2=indiv[0][random.randint(0,indiv[0].get_number_of_atoms()-1)]
            Optimizer.output.write('WARNING: Permutation Mutation attempted on single atom structure system\n')
        solid[a1.index].symbol, solid[a2.index].symbol=a2.symbol, a1.symbol
        indiv.bulki=solid[indexbulk::]
        indiv[0]=solid[0:indexbulk]
        Optimizer.output.write('Crystal-wide Permutation Mutation performed on individual\n')
        Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
        Optimizer.output.write(repr(a1.symbol) + ' swapped with ' + repr(a2.symbol)+'\n')
        Optimizer.output.write(repr(indiv[0])+'\n')
        for sym,c,m,u in Optimizer.atomlist:
            nc=len([atm for atm in solid if atm.symbol==sym])
            Optimizer.output.write('Defect configuration contains '+repr(nc)+' '+repr(sym)+' atoms\n')
        muttype='PC'
        if indiv.energy==0:
            indiv.history_index=indiv.history_index+'m'+muttype
        else:
            indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv