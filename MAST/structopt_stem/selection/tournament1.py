import random

def tournament1(pop, nkeep, Optimizer):
    """Selection function that orders and chooses individuals based on tournament schemes
    Function will attempt to prevent repeats
    """
    tournsize = Optimizer.tournsize
    fits = [ind.fitness for ind in pop]
    minindex = min(xrange(len(fits)), key=fits.__getitem__)
    newpop = []
    indices = []
    indexo = [ind.index for ind in pop]
    for i in range(len(pop)):
        pop[i].index = i
    counter=0
    while len(newpop)< nkeep:
        subgroup = random.sample(pop,tournsize)
        mine = min([ind.fitness for ind in subgroup])
        select = [ind for ind in subgroup if ind.fitness==float(mine)][0]
        if select.index not in indices:
            newpop.append(select)
            indices.append(select.index)
        elif counter > 100000:
            try:
                nonpop = [ind for ind in pop if ind.index not in indices]
                select = random.choice(nonpop)
                newpop.append(select)
                indices.append(select.index)
            except:
                select = random.choice(pop)
                newpop.apend(select)
                indices.append(select.index)
        counter+=1
    for i in range(len(newpop)):
        newpop[i].index = indexo[indices[i]]
    if minindex not in indices:
        rm = random.choice(indices)
        newpop = [pop[inx] for inx in indices if inx != rm]
        newpop.append(pop[minindex])
    return newpop