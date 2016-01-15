import random

def fuss1r(pop, nkeep, Optimizer):
	"""Selection function to employ basic Fixed Uniform Selection Scheme to select and 
	order individuals in a population.
	No constraint on fuss -> population will diverge
	"""
	newpop = []
	#Collect fitnesses
	fits = [ind.fitness for ind in pop]
	#Identify the minimum and maximum fitness and index
	minindex = min(xrange(len(fits)), key=fits.__getitem__)
	minf = min(fits)
	maxf = max(fits)
	#Scale fitnesses to range 0-1
	fits = [(f-minf)/(maxf-minf) for f in fits]
	indices = []
	while len(newpop) < nkeep:
		pt = random.random()
		if debug: print pt
		distances = [abs(pt-fit) for fit in fits]
		index = min(xrange(len(distances)), key=distances.__getitem__)
		if index not in indices:
			newpop.append(pop[index])
			indices.append(index)
	return newpop