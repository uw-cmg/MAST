import copy
import random
import numpy
from MAST.structopt.tools import fitness_switch
from MAST.structopt.moves.atoms_add import atoms_add
from MAST.structopt.moves.atoms_remove import atoms_remove

def basin_hop_ra_atoms(indiv, Optimizer):
    """Move function to perform mini-basin hopping run to add/remove atoms
    Inputs:
        indiv = Individual class object to be altered
        Optimizer = Optimizer class object with needed parameters
    Outputs:
        indiv = Altered Individual class object
    *** needs work ***
    """
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False
    Optimizer.output.write('Performing Basin Hopping Add/Remove Atoms Mutation on indiv '+repr(indiv.index)+'\n')
    #Save copy of starting positions
    startindiv=indiv[0].copy()
    alloysetting=Optimizer.alloy
    fingerprintsetting=Optimizer.fingerprinting
    if Optimizer.structure=='Defect':
        finddefectsetting=Optimizer.finddefects
        startbulk=indiv.bulki
        indivmark=len(startindiv)
    mutation_options = copy.deepcopy(Optimizer.mutation_options)
    startindiv = indiv.duplicate()
    
    #Get starting fitness
    fitmin=indiv.fitness
    if indiv.fitness==0:
        out=fitness_switch([Optimizer,indiv])
        fitmin=out[0].fitness
    options = ['atoms_add', 'atoms_remove']
    
    totalsteps=Optimizer.bh_steps
    kt=Optimizer.bh_temp
    flag=False
    for step in range(totalsteps):
        prevind = indiv.duplicate()
        scheme = random.choice(options)
        indiv = eval(scheme+'(indiv, Optimizer)')
        out=fitness_switch([Optimizer,indiv])
        fitnew=out[0].fitness
        if fitnew < fitmin:
            flag=True
            break
        else:
            accept=numpy.exp((fitmin - fitnew) / kt) > random.random()
            if not accept:
                indiv = prevind

    Optimizer.output.write('Evaluated '+repr(step)+' steps\n')
    if flag==False:
        Optimizer.output.write('Failed to find lower energy add/removed atoms structure\n')
        #indiv = startindiv
    else:
        Optimizer.output.write('Found lower energy add/removed atoms structure\n')
    Optimizer.alloy=alloysetting
    Optimizer.fingerprinting=fingerprintsetting
    if Optimizer.structure=='Defect':
        Optimizer.finddefects=finddefectsetting
    Optimizer.mutation_options = mutation_options
    muttype='BHAR'+repr(step)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv
