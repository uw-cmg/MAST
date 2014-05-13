from MAST.structopt.fingerprinting import fingerprint_dist
import random

def fussf(pop, nkeep, Optimizer):
	"""Selection function to employ basic Fixed Uniform Selection Scheme to select and 
	order individuals in a population based on fingerprint distance values.
	Picks one fuss point for the selection process
	FUSS point must reside within fusslimit
	"""
	fusslimit = Optimizer.fusslimit
	#FUSS for fingerprint
	newpop = []
	# Collect fitnesses
	fits = [ind.fitness for ind in pop]
	minindex = min(xrange(len(fits)), key=fits.__getitem__)
	fits = []
	STR = 'Distances: \n'
	for i in range(len(pop)):
		dist = fingerprint_dist(pop[minindex].fingerprint,pop[i].fingerprint)
		fits.append(dist)
		STR += repr(dist)+'\n'
	
	#Find min and max fitness
	minf = min(fits)
	maxf = max(fits)
	if abs(maxf-minf) > fusslimit:
		maxf = minf + fusslimit
	#Select random point on fitness line
	pt = random.uniform(minf,maxf)
	#Calculate the distance of each individual's fitness from that point
	i = 0
	distances = []
	for fit in fits:
		distances.append((abs(pt-fit),i))
		i += 1
	#Sort distances from min to max
	dlist = sorted(distances, key=lambda one: one[0])
	#Select individuals with lowest distance
	indices = []
	for i in range(nkeep):
		newpop.append(pop[dlist[i][1]])
		indices.append(dlist[i][1])
	#Always keep lowest energy individual
	if minindex not in indices:
		rm = random.choice(indices)
		newpop = [pop[inx] for inx in indices if inx != rm]
		newpop.append(pop[minindex])
	return newpop
