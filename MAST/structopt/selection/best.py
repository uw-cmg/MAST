from MAST.structopt.tools import get_best

def best(pop, nkeep, Optimizer=None):
    """Selection function to select the best individuals in a population
    Inputs:
        pop = list of Individual class structures
        nkeep = number of Individuals to keep in the population
            must be less than or equal to the number of individuals currently in the population
        Optimizer = dummy placeholder
    Outputs:
        npop = population with desired number of individuals
    """
    if nkeep <= len(pop):
        npop = get_best(pop,nkeep)
    else:
        print 'WARNING: selection.best nkeep > pop. Returning pop'
        npop = pop
    return npop