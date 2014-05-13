from MAST.structopt.tools import get_best
import random

def cost(pop, nkeep, Optimizer):
	"""Selection function to order and select structures based on simple relative cost
	"""
	pop = get_best(pop,len(pop))
	newpop = []
	fitnesses = [ind.fitness for ind in pop]
	prob = []
	cumprob = []
	try:
		norms = [fit - pop[nkeep+1].fitness for fit in fitnesses[0:nkeep]]
	except:
		norms = [fit - pop[nkeep-1].fitness for fit in fitnesses[0:nkeep]]
	sumn = sum(norms)
	for i in range(nkeep):
		prob.append(norms[i] / sumn)
		cumprob.append(sum(prob))
	prevcounter = []
	for i in range(nkeep):
		rand = random.random()
		counter = 0
		while cumprob[counter] < rand:
			counter += 1
		if counter in prevcounter:
			a = random.choice(pop)
			newpop.append(a)
			prevcounter.append(a.index)
		else:
			newpop.append(pop[counter])
			prevcounter.append(counter)
	return newpop