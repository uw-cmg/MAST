from operator import attrgetter

def get_best(pop, nkeep):
    """Function for sorting population by fitness attributes"""
    pop = sorted(pop, key=attrgetter('fitness')) 
    pop = pop[0:nkeep]
    return pop
