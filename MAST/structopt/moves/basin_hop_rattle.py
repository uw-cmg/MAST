import numpy
import copy
import random
from MAST.structopt.switches import fitness_switch
from MAST.structopt.moves.rattle import rattle

def basin_hop_rattle(indiv, Optimizer):
    """Move function to perform mini-basin hopping run to Rattle atoms
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
    Optimizer.output.write('Performing Basin Hopping Rattle Mutation on indiv '+repr(indiv.index)+'\n')
    #Save copy of starting positions
    startindiv=indiv[0].copy()
    startbulk=indiv.bulki
    indivmark=len(startindiv)
    alloysetting=Optimizer.alloy
    fingerprintsetting=Optimizer.fingerprinting
    if Optimizer.structure=='Defect':
        finddefectsetting=Optimizer.finddefects
    mutation_options = copy.deepcopy(Optimizer.mutation_options)
    startindiv = indiv.duplicate()
    
    #Get starting fitness
    fitmin=indiv.fitness
    if indiv.fitness==0:
        out=fitness_switch([Optimizer,indiv])
        fitmin=out[0].fitness
    options = ['rattle']
    
    totalsteps=Optimizer.bh_steps
    kt=Optimizer.bh_temp
    flag=False
    rattlemin = None
    for step in range(totalsteps):
        prevind = indiv.duplicate()
        scheme = random.choice(options)
        indiv = eval(scheme+'(indiv, Optimizer)')
        out = fitness_switch([Optimizer,indiv])
        fitnew=out[0].fitness
        indiv = out[0]
        if fitnew < fitmin:
            flag=True
            break
        elif fitnew > (fitmin+Optimizer.tolerance*2):
            accept=numpy.exp((fitmin - fitnew) / kt) > random.random()
            if not rattlemin:
                rattlemin = indiv.duplicate()
            else:
                if fitnew < rattlemin.fitness:
                    rattlemin = indiv.duplicate()
            if not accept:
                indiv = prevind
    if rattlemin:
        indiv = rattlemin.duplicate()
    Optimizer.output.write('Evaluated '+repr(step)+' steps\n')
    if flag==False:
        Optimizer.output.write('Failed to find lower energy rattled structure\n')
        #indiv = startindiv
    else:
        Optimizer.output.write('Found lower energy rattled structure\n')
    Optimizer.alloy=alloysetting
    Optimizer.fingerprinting=fingerprintsetting
    if Optimizer.structure=='Defect':
        Optimizer.finddefects=finddefectsetting
    Optimizer.mutation_options = mutation_options
    muttype='BHRat'+repr(step)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv