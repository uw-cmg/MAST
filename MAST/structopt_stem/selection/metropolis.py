import random
import math

def metropolis(pop, nkeep, Optimizer):
	"""Selection function that applies metropolis algorithm to determine which individuals to keep
	Inputs:
		pop = initial list of Individual class structures
		nkeep = number of structure to keep in population
			nkeep must be less than or equal to pop
		Optimizer = Optimizer class object with metropolis temperature
	Outputs:
		newpop = new list of Individual class structures of length nkeep
	"""
	newpop = []
	#Identify temperature for algorithm
	T = Optimizer.metropolis_temp
	#Collect fitnesses
	fits = [ind.fitness for ind in pop]
	#Identify the minimum fitness
	minf = min(fits)
	print fits
	while len(newpop)<nkeep:
		ind = random.choice(pop)
		r=random.random()
		if r <= math.exp((minf-ind.fitness)/T):
			newpop.append(ind)
	return newpop