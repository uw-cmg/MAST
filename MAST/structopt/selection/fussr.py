import random

def fussr(pop,nkeep,Optimizer):
	"""Selection function to employ basic Fixed Uniform Selection Scheme to select and 
	order individuals in a population.
	Scales fuss values
	Picks one fuss point for the selection process
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
	#Pick FUSS point and calculate distances for each fitness
	pt = random.random()
	distances = [[abs(pt-fits[i]),i] for i in range(len(fits))]
	distances.sort()
	indices = []
	for i in range(nkeep):
		newpop.append(pop[distances[i][1]])
		indices.append(distances[i][1])
	if minindex not in indices:
		rm = random.choice(indices)
		newpop = [pop[inx] for inx in indices if inx != rm]
		newpop.append(pop[minindex])
	return newpop