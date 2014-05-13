#Test the BestInds function

from MAST.structopt.tools.BestInds import BestInds
from MAST.structopt.generate.Individual import Individual
from ase import Atom, Atoms
from operator import attrgetter

class Optifake():
    def __init__(self, Number_of_Bests, tolerance, filename, structure):
        self.Number_of_Bests = Number_of_Bests
        self.demin = tolerance
        self.algorithm_type = 'lambda+mu'
        self.filename = filename
        self.structure = structure

def test_BestInds():
    tolerance = 0.99
    A = Optifake(10,tolerance,'Test','Cluster')
    bests = []
    for j in range(10):
        pop = []
        for i in range(5):
            pop.append(Individual(Atoms(),fitness=float(i)*0.9-j*1.3))
        bests = BestInds(pop,bests,A)
        bfits = [ind.fitness for ind in bests]
        pfits = [ind.fitness for ind in pop]

    diff = [abs(bfits[i]-bfits[i+1]) for i in range(len(bfits)-1)]
    if min(diff) < tolerance:
        print 'Error in BestInds test'