import random
import numpy
from ase import Atom, Atoms
from MAST.structopt.tools.find_defects import find_defects
from MAST.structopt.generate import gen_pop_box

def random_replacement(indiv, Optimizer):
    """Move function to replace selection of atoms with randomly generated group
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
            atms,indb,vacant,swap,stro = find_defects(indiv[0],Optimizer.solidbulk,0)
        else:
            atms = indiv[0].copy()
    else:
        atms=indiv[0].copy()
    nat=len(atms)
    if nat != 0:
        #Select number of atoms to replace
        if nat<=1:
            natrep2=1
            natrep1=0
        elif nat<=5:
            natrep2=2
            natrep1=0
        else: 
            natrep1=random.randint(1,nat/2)
            while True:
                natrep2=random.randint(2,nat/2)
                if natrep2 != natrep1:
                    break
        natrep=abs(natrep2 - natrep1)
        #Select random position in cluster
        maxcell = numpy.maximum.reduce(atms.get_positions())
        mincell = numpy.minimum.reduce(atms.get_positions())
        pt=[random.uniform(mincell[0],maxcell[0]),random.uniform(mincell[1],maxcell[1]),random.uniform(mincell[2],maxcell[2])]
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
        atmsr=Atoms()
        indexlist=[]
        for i in range(natrep):
            atmsr.append(dlist[i][1])
            indexlist.append(dlist[i][1].index)
        natomlist=[0]*len(Optimizer.atomlist)
        for i in range(len(Optimizer.atomlist)):
            atms1=[inds for inds in atmsr if inds.symbol==Optimizer.atomlist[i][0]]
            natomlist[i]=(Optimizer.atomlist[i][0], len(atms1),Optimizer.atomlist[i][2],Optimizer.atomlist[i][3])
        nsize = max(numpy.maximum.reduce(atmsr.get_positions())-numpy.minimum.reduce(atmsr.get_positions()))
        repcenter = atmsr.get_center_of_mass()
        atmsn = gen_pop_box(natomlist,nsize)
        atmsn.translate(repcenter)
        #Update individual with new atom positions
        for i in range(len(indexlist)):
            index=indexlist[i]
            atms[index].position=atmsn[i].position
        if Optimizer.structure=='Defect':
            if Optimizer.isolate_mutation:
                atms.extend(indb)
        indiv[0]=atms.copy()
    else:
        natrep=0
        pt=0
    Optimizer.output.write('Random Group Replacement Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write('Number of atoms replaced = '+repr(natrep)+'\n')
    Optimizer.output.write('Geometry point = '+repr(pt)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='RGR'+repr(natrep)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv
