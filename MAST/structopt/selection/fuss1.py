import random

def fuss1(pop, nkeep, Optimizer):
	"""Selection function to employ basic Fixed Uniform Selection Scheme to select and 
	order individuals in a population.
	No constraint on fuss -> population will diverge
	"""
	newpop = []
	fits = [ind.fitness for ind in pop]
	minindex = min(xrange(len(fits)), key=fits.__getitem__)
	minf = min(fits)
	maxf = max(fits)
	indices = []
	while len(newpop) < nkeep:
		pt = random.uniform(minf,maxf)
		distances = [abs(pt-fit) for fit in fits]
		index = min(xrange(len(distances)), key=distances.__getitem__)
		if index not in indices:
			newpop.append(pop[index])
			indices.append(index)
	if minindex not in indices:
		rm = random.choice(indices)
		newpop = [pop[inx] for inx in indices if inx != rm]
		newpop.append(pop[minindex])
	return newpop