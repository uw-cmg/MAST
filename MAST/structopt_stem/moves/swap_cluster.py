import random
from ase import Atoms
import math
from math import fabs, sqrt, acos, pi, atan
import numpy 
from ase.calculators.neighborlist import NeighborList
import logging

def swap_cluster(indiv, Optimizer):
    """swap two clusters of atoms whose central atom are within low potential energy list
    Inputs:
        indiv = Individual class object to be altered; including indiv.lpealist giving low PE atom index
        Optimizer = Optimizer class object with needed parameters
    Outputs:
        indiv = Altered Individual class object
    """
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False

    logger = logging.getLogger(Optimizer.loggername)
    #logger.info('M:start swap-cluster')
    Optimizer.output.write('swap cluster Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')

    syms = [sym for sym,c,m,u in Optimizer.atomlist]
    numatom = [c for sym,c,m,u in Optimizer.atomlist]

    #totalsol = indiv[0] 
    #cutoffs=[1.5 for one in totalsol]
    #nl=NeighborList(cutoffs,bothways=True,self_interaction=False)
    #nl.update(totalsol)
    #nat=len(totalsol)
    #logger.info('M:NeighborList')
    positions=indiv[0].get_positions()

    lpealist = indiv.lpealist

    atmlist = [] 

    atom1 = random.choice(lpealist)
    lpealist.remove(atom1)
    atom2 = random.choice(lpealist)

    #indices1, offsets=nl.get_neighbors(atom1)
    #indices2, offsets=nl.get_neighbors(atom2)

    #print 'atom1',atom1,'indices1',indices1
    #print 'atom2',atom2,'indices2',indices2    

    indices1 = []
    for i in range(len(positions)):
        R2 = (positions[i][0]-positions[atom1][0])**2 + (positions[i][1]-positions[atom1][1])**2 + (positions[i][2]-positions[atom1][2])**2
        if R2 < 3.8**2 and i != atom1 :
          indices1.append(i)

    indices2 = []
    for i in range(len(positions)):
        R2 = (positions[i][0]-positions[atom2][0])**2 + (positions[i][1]-positions[atom2][1])**2 + (positions[i][2]-positions[atom2][2])**2
        if R2 < 3.8**2 and i != atom2 :
          indices2.append(i)
    #print 'indices1',indices1
    #print 'indices2',indices2
    #logger.info('M:Neighboratom')

    for ind in indices1:
       positions[ind][0] += positions[atom2][0]-positions[atom1][0]
       positions[ind][1] += positions[atom2][0]-positions[atom1][0]
       positions[ind][2] += positions[atom2][0]-positions[atom1][0]
       
    for ind in indices2:
       positions[ind][0] += positions[atom1][0]-positions[atom2][0]
       positions[ind][1] += positions[atom1][0]-positions[atom2][0]
       positions[ind][2] += positions[atom1][0]-positions[atom2][0]

    position_tmp = positions[atom1]
    positions[atom1] = positions[atom2]
    positions[atom2] = positions[atom1] 
       
    indiv[0].set_positions(positions) 
    muttype='swap-cluster'
    Optimizer.output.write('moved atoms='+repr(len(indices1)+1)+'+'+repr(len(indices2)+1)+'\n')
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    #logger.info('M:end swapcluster') 
    return indiv
