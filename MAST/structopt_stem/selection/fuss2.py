import random

def fuss2(pop, nkeep, Optimizer):
    """Selection function to employ basic Fixed Uniform Selection Scheme to select and 
	order individuals in a population.
	fuss point must reside within fusslimit of minimum fitness
	"""
    newpop = []
    fits = [ind.fitness for ind in pop]
    minindex = min(xrange(len(fits)), key=fits.__getitem__)
    minf = min(fits)
    maxf = max(fits)
    if abs(maxf - minf) > fusslimit:
        maxf = minf + fusslimit
    indices = []
    loopcount = 0
    while len(newpop) < nkeep:
        pt = random.uniform(minf,maxf)
        i = 0
        distances = []
        for fit in fits:
            distances.append((abs(pt - fit),i))
            i += 1
        dlist = sorted(distances, key=lambda one: one[0])
        index = dlist[0][1]
        if index not in indices:
            newpop.append(pop[index])
            indices.append(index)
        else:
            loopcount += 1
        if loopcount > 100:
            indexl = [l for l in range(len(pop)) if l not in indices]
            if len(indexl)==0:
                randi = random.choice(range(len(pop)))
                newpop.append(pop[randi])
                indices.append(randi)
            else:
                newpop.append(pop[min(indexl)])
                indices.append(min(indexl))
    if minindex not in indices:
        rm = random.choice(indices)
        newpop = [pop[inx] for inx in indices if inx != rm]
        newpop.append(pop[minindex])
    return newpop