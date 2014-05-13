import random
import numpy
from ase import Atom, Atoms
from MAST.structopt.tools.find_defects import find_defects

def permutation_crystal_multi(indiv, Optimizer):
    """Move function to perform Permutation of multiple atom based on atomlist with entire crystal in Defect
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
        solid.set_cell(indiv[0].get_cell())
        a1=solid[random.randint(0,indiv[0].get_number_of_atoms()-1)]
        opts=[inds for inds in solid if inds.symbol != a1.symbol]
        inpts=[inds for inds in solid if inds.symbol == a1.symbol]
        maxperm=len(inpts)
        if maxperm > len(opts):
            maxperm=len(opts)
        Optimizer.output.write('Multiatom Crystal-wide Permutation Mutation performed on individual\n')
        Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
        if maxperm==0:
            Optimizer.output.write('WARNING: Permutation Mutation attempted on single atom structure system\n')
            natperm=0
        else:
            natperm=random.randint(1,maxperm/2)
        Optimizer.output.write('Number of atoms permutated = '+repr(natperm)+'\n')
        for i in range(natperm):
            a1=inpts[random.randint(0,len(inpts)-1)]
            a2=opts[random.randint(0,len(opts)-1)]
            solid[a1.index].symbol, solid[a2.index].symbol=a2.symbol, a1.symbol
            Optimizer.output.write(repr(a1.symbol) + ' swapped with ' + repr(a2.symbol)+'\n')
        indiv.bulki=solid[indexbulk::]
        indiv[0]=solid[0:indexbulk]
        Optimizer.output.write(repr(indiv[0])+'\n')
        for sym,c,m,u in Optimizer.atomlist:
            nc=len([atm for atm in solid if atm.symbol==sym])
            Optimizer.output.write('Defect configuration contains '+repr(nc)+' '+repr(sym)+' atoms\n')
        muttype='PCMA'+repr(natperm)
        if indiv.energy==0:
            indiv.history_index=indiv.history_index+'m'+muttype
        else:
            indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv