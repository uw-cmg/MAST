import random

def tournament2(pop, nkeep, Optimizer):
    """Selection function to order and choose individuals based on Tournament schemes
    Function will allow repeats but not back to back
    """
    tournsize = Optimizer.tournsize
    newpop = []
    prevind = None
    while len(newpop) < nkeep:
        subgroup = []
        for i in range(tournsize):
            subgroup.append(random.choice(pop))
        fitnesses = [ind.fitness for ind in subgroup]
        mine = min(fitnesses)
        select = [ind for ind in subgroup if ind.fitness==mine]
        if select[0].index != prevind:
            newpop.append(select[0])
            prevind = select[0].index
    return newpop