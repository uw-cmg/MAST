import random
from MAST.structopt.tools import get_best

def rank(pop, nkeep, Optimizer):
    """Selection function that chooses structures to survive and orders them based on their relative ranking
    """
    pop = get_best(pop,len(pop))
    newpop = []
    prob = []
    cumprob = []
    for i in range(len(pop)):
        rankprob = float((nkeep-i)) / float(sum(range(nkeep+1)))
        prob.append(rankprob)
        cumprob.append(sum(prob))
    prevcounter = None
    for i in range(nkeep):
        rand = random.random()
        counter = 0
        while cumprob[counter] < rand:
            counter += 1
        if prevcounter==counter:
            while True:
                a = random.choice(pop)
                if a.index != prevcounter:
                    break
            newpop.append(random.choice(pop))
            prevcounter = a.index
        else:
            newpop.append(pop[counter])
            prevcounter = counter
    return newpop