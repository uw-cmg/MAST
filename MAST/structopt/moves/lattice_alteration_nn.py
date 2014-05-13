import random
import numpy
from MAST.structopt.tools.find_defects import find_defects
from ase import Atom, Atoms
from ase.calculators.neighborlist import NeighborList

def lattice_alteration_nn(indiv, Optimizer):
    """Move function to perform random move along random axis for nearest neighbor distance to random atoms
    Inputs:
        indiv = Individual class object to be altered
        Optimizer = Optimizer class object with needed parameters
    Outputs:
        indiv = Altered Individual class object
    """
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False
    if Optimizer.structure=='Defect':
        if Optimizer.isolate_mutation:
            indc,indb,vacant,swaps,stro = find_defects(indiv[0],Optimizer.solidbulk,0)
            ind = indc.copy()
            ind.extend(indb)
        else:
            ind=indiv[0].copy()
            indc=indiv[0].copy()
    else:
        ind=indiv[0].copy()
        indc=indiv[0].copy()
    if len(indc) != 0:
        ctoff1 = [1.0 for one in ind]
        nl = NeighborList(ctoff1, bothways=True, self_interaction=False)
        nl.update(ind)
        try:
            natomsmove=random.randint(1,len(indc)/2)
        except ValueError:
            natomsmove=1
        passn=0
        for count in range(natomsmove):
            try:
                indexmv = random.choice([i for i in range(len(indc))])
                indices, offsets = nl.get_neighbors(indexmv)
                nns = Atoms()
                nns.append(ind[indexmv])
                for index, d in zip(indices,offsets):
                    index = int(index)
                    pos = ind[index].position + numpy.dot(d,ind.get_cell())
                    nns.append(Atom(symbol=ind[index].symbol, position=pos))
                dist = [nns.get_distance(0,i) for i in range(1, len(nns))]
                r = sum(dist)/len(dist)
                dir = random.choice([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
                ind[indexmv].position += [i*r for i in dir]
            except:
                passn+=1
        indiv[0]=ind.copy()
    else:
        natomsmove=0
        passn=0
    Optimizer.output.write('Lattice Alteration NN Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    natomsmove-=passn
    Optimizer.output.write('Number of atoms moved = '+repr(natomsmove)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='LANN'+repr(natomsmove)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv