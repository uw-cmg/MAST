from MAST.structopt.tools import get_best
from MAST.structopt.switches import selection_switch

def lambdacommamu(pop, Optimizer):
    """Selection function to employ a lambda,mu GA scheme
    """
    fits = [ind.fitness for ind in pop]
    minindex = [i for i in range(len(pop)) if pop[i].fitness==min(fits)]
    try:
        Optimizer.mark
    except:
        Optimizer.mark = len(pop)/2
    parents = pop[0:Optimizer.mark]
    offspring = pop[Optimizer.mark::]
    if minindex < Optimizer.mark:
        offspring.append(pop[minindex])
    if len(offspring) < Optimizer.nindiv:
        diff = Optimizer.nindiv-len(offspring)
        if Optimizer.natural_selection_scheme=='elitism':
            addins = get_best(parents,diff)
            STR = 'Adding in '+repr(diff)+' lowest fitness parents\n'
        else:
            addins = selection_switch(parents, diff, Optimizer.natural_selection_scheme, Optimizer)
            STR = 'Adding in '+repr(diff)+' parents based on natural selection\n'
        for one in addins:
            offspring.append(one)
    elif len(offspring) > Optimizer.nindiv:
        diff = len(offspring)-Optimizer.nindiv
        if Optimizer.natural_selection_scheme=='elitism':
            offspring = get_best(offspring,Optimizer.nindiv)
            STR = 'Removing lowest '+repr(diff)+' fitness offspring\n'
        else:
            offspring = selection_switch(offspring, Optimizer.nindiv, Optimizer.natural_selection_scheme, Optimizer)
            STR = 'Removing '+repr(diff)+' offspring by natural selection\n'
    else:
        STR = 'Number of offspring = {0}\n'.format(len(offspring))
    return offspring, STR