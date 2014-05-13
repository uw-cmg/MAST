import random
from MAST.structopt.tools import get_best

def random_pick(pop, nkeep, Optimizer):
	"""Selection function that randomly choose structures from a population to survive.
	"""
	pop = get_best(pop,len(pop))
	newpop = []
	for i in range(nkeep):
		sel = random.choice(pop)
		newpop.append(sel)
	return newpop