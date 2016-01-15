import random
import numpy
import copy
from ase import Atom, Atoms
from ase.calculators.neighborlist import NeighborList

def swap_int_local(indiv, Optimizer):
    """Move function to perform n atom structure swap based on swaplist for cluster in local region
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
    Optimizer.output.write('Local Cluster Swap Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    # Identify the local atoms available for swapping
    atmsst = indiv[0].copy()
    nat = len(atmsst)
    totalsol = atmsst.copy()
    if Optimizer.structure=='Defect':
        totalsol.extend(indiv.bulki.copy())
        ctoff = [Optimizer.sf for one in totalsol]
    else:
        ctoff = [1.5 for one in totalsol]
    nl = NeighborList(ctoff, bothways=True, self_interaction=False)
    nl.update(totalsol)
    cluster = Atoms()
    indices = [] #store indices of the atoms for later
    for i in range(nat):
        indices1, offsets = nl.get_neighbors(i)
        for index, d in zip(indices1,offsets):
            index = int(index)
            #Make sure to prevent repeats of the same atom
            if index not in indices:
                pos = totalsol[index].position + numpy.dot(d,totalsol.get_cell())
                cluster.append(Atom(symbol=totalsol[index].symbol,position=pos))
                indices.append(index)
    # Identify Number of atoms to swap
    swapmax=sum([c for sym,c in indiv.swaplist])
    if swapmax >1:
        natomsswap=random.randint(1,swapmax)
    elif swapmax==1:
        natomsswap=1
    else:
        Optimizer.output.write('Maximum number of swaps is '+repr(swapmax)+'. Unable to perform mutation')
        natomsswap=0
    if len(indiv.swaplist)==1:
        Optimizer.output.write('WARNING: Swap Mutation attempted on single atom structure system\n')
        natomsswap=0
    Optimizer.output.write('Number of swaps = '+repr(natomsswap)+'\n')
    # Identify symbols that can be swapped
    syms=[sym for sym,c in indiv.swaplist]
    slist = copy.deepcopy(indiv.swaplist)
    if debug: print 'Starting swaplist = '+repr(indiv.swaplist)
    #Sanity check
    sanch = [[sym,0] for sym in syms]
    for i in range(len(sanch)):
        nc = len([atm for atm in cluster if atm.symbol==sanch[i][0]])
        nc += [c for sym,c in slist if sym==sanch[i][0]][0]
        sanch[i][1]=nc
    # Swap Atoms
    for i in range(natomsswap):
        while True:
            at1 = random.choice(cluster)
            osym = at1.symbol
            nsyml=[sym for sym,c in slist if sym != osym and c > 0]
            if len(nsyml) > 0:
                break
        nsym = random.choice(nsyml)
        cluster[at1.index].symbol = nsym
        Optimizer.output.write('Swapped '+nsym+' atom with '+osym+'\n')
        for i in range(len(slist)):
            if slist[i][0]==nsym:
                slist[i][1]-=1
            elif slist[i][0]==osym:
                slist[i][1]+=1
    # Apply the new swap to the actual atoms in the structure
    for i in range(len(indices)):
        totalsol[indices[i]].symbol = cluster[i].symbol
    # Separate out bulk from defects
    indiv[0] = totalsol[0:nat]
    if Optimizer.structure=='Defect':
        indiv.bulki = totalsol[nat::]
    # Update new swaplist
    indiv.swaplist = slist
    #Sanity check
    sanchn = [[sym,0] for sym in syms]
    for i in range(len(sanchn)):
        nc = len([atm for atm in cluster if atm.symbol==sanchn[i][0]])
        nc += [c for sym,c in slist if sym==sanchn[i][0]][0]
        sanchn[i][1]=nc
    if debug:
        for i in range(len(sanch)):
            print 'INTSWAP: Starting Atom structure '+repr(sanch[i][0])+' = '+repr(sanch[i][1])
            print 'INTSWAP: After Mutation Atom structure '+repr(sanchn[i][0])+' = '+repr(sanchn[i][1])
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='SCL'+repr(natomsswap)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv