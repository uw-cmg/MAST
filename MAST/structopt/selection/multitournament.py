import random

def multitournament(pop, nkeep, Optimizer):
    """Multiple child tournament scheme for selection
    """
    tournsize = Optimizer.tournsize
    #Will allow repeats but not back to back
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
    for i in range(nkeep/4):
        newpop.append(pop[0])
        newpop.append(pop[random.randint(1,len(pop)-1)])
    return newpop