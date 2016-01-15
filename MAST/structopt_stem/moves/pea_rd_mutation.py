import random
from ase import Atoms
import math
from math import fabs, sqrt, acos, pi, atan
import numpy 

def pea_rd_mutation(indiv, Optimizer):
    """Move function to move atoms with high PE in random direction for random distance
    Inputs:
        indiv = Individual class object to be altered; including indiv.hpealist giving high PE atom index
        Optimizer = Optimizer class object with needed parameters
    Outputs:
        indiv = Altered Individual class object
    """
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False

    Optimizer.output.write('pea_rd Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')

    syms = [sym for sym,c,m,u in Optimizer.atomlist]
    numatom = [c for sym,c,m,u in Optimizer.atomlist]
    positions=indiv[0].get_positions()

    movelist = indiv.hpealist
    try:
        natomsmove=random.randint(1,len(movelist)/len(syms))
    except ValueError:
        natomsmove=1

    atmlist = [] 
    for i in range(natomsmove): 
       iatom = random.choice(movelist)
       atmlist.append(iatom)
    list(set(atmlist))

    d_max=numpy.minimum.reduce(numpy.maximum.reduce(indiv[0].get_positions())-numpy.minimum.reduce(indiv[0].get_positions()))
    r=random.uniform(0.3,d_max)

    for i in range(len(atmlist)):
       ind = atmlist[i]
       theta=math.radians(random.uniform(0,360))
       phi=math.radians(random.uniform(0,180))
       direction=[r*math.sin(theta)*math.cos(phi),r*math.sin(theta)*math.sin(phi),r*math.cos(theta)]    
       positions[ind][0] += direction[0]
       positions[ind][1] += direction[1]
       positions[ind][2] += direction[2]
       
    indiv[0].set_positions(positions) 
    muttype='pea-rd'+repr(natomsmove)
    Optimizer.output.write('moved atoms='+repr(atmlist)+'\n')
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype

    return indiv
