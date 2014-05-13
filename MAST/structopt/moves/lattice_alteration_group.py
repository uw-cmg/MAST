import random
import numpy
from ase import Atom, Atoms
from MAST.structopt.tools.find_defects import find_defects

def lattice_alteration_group(indiv, Optimizer):
    """Move function to perform Lattice Alteration of group of atoms based on location
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
            indc,indb,vacant,swap,stro = find_defects(indiv[0],Optimizer.solidbulk,0)
            atms=indc.copy()
        else:
            atms=indiv[0].copy()
    else:
        atms=indiv[0].copy()
    if len(atms) != 0:
        try:
            natomsmove=random.randint(1,len(atms)/5)
        except ValueError:
            natomsmove=1
        #Select random position in cluster
        cellx=numpy.maximum.reduce(atms.get_positions())
        cellm=numpy.minimum.reduce(atms.get_positions())
        pt=[random.uniform(cellm[0],cellx[0]),random.uniform(cellm[1],cellx[1]),random.uniform(cellm[2],cellx[2])]
        #Get distance of atoms from random point
        atpt=Atom(position=pt)
        atms.append(atpt)
        dist=[]
        for i in range(len(atms)-1):
            dist.append(atms.get_distance(i,len(atms)-1))
        atms.pop()
        dlist=zip(dist,atms)
        dlist=sorted(dlist, key=lambda one: one[0], reverse=True)
        # Select atoms closest to random point
        atmst=Atoms()
        indexlist=[]
        for i in range(natomsmove):
            atmst.append(dlist[i][1])
            indexlist.append(dlist[i][1].index)
        trans=(random.uniform(0,cellx[0]-cellm[0]), random.uniform(0,cellx[1]-cellm[1]), random.uniform(0,cellx[2]-cellm[2]))
        atmst.translate(trans)
        for i in range(len(indexlist)):
            index=indexlist[i]
            atms[index].position=atmst[i].position
        if Optimizer.structure=='Defect':
            if Optimizer.isolate_mutation:
                indiv[0]=atms.copy()
                indiv[0].extend(indb)
            else:
                indiv[0] = atms.copy()
        else:
            indiv[0]=atms.copy()
    else:
        natomsmove=0
        trans=0
    Optimizer.output.write('Group Lattice Alteration Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write('Number of atoms moved = '+repr(natomsmove)+'\n')
    Optimizer.output.write(repr(trans)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='LAGC'+repr(natomsmove)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv