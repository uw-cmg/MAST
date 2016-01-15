import random

def tournament(pop, nkeep, Optimizer):
    """Selection function to order and choose individuals based on Tournament schemes
    Function will allow for repeats
    """
    tournsize = Optimizer.tournsize
    newpop = []
    for j in range(nkeep):
        subgroup = []
        for i in range(tournsize):
            subgroup.append(random.choice(pop))
        fitnesses = [ind.fitness for ind in subgroup]
        mine = min(fitnesses)
        select = [ind for ind in subgroup if ind.fitness==mine]
        newpop.append(select[0])
    return newpop