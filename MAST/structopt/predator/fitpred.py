from MAST.structopt.tools import get_best

def fitpred(pop,Optimizer):
    """Predator function to select best structures
    """
    pop = get_best(pop,Optimizer.nindiv)
    return pop, ''
