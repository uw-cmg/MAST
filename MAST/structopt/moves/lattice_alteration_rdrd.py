import random
import math
import numpy
from MAST.structopt.tools.find_defects import find_defects

def lattice_alteration_rdrd(indiv, Optimizer):
    """Move function to move random atoms in random direction for random distance
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
    d_max=numpy.minimum.reduce(numpy.maximum.reduce(indiv[0].get_positions())-numpy.minimum.reduce(indiv[0].get_positions()))
    if Optimizer.structure=='Defect':
        if Optimizer.isolate_mutation:
            indc,indb,vacant,swap,stro = find_defects(indiv[0],Optimizer.solidbulk,0)
            positions=indc.get_positions()
        else:
            positions = indiv[0].get_positions()
    else:
        positions=indiv[0].get_positions()
    if len(positions) != 0:
        try:
            natomsmove=random.randint(1,len(positions)/5)
        except ValueError:
            natomsmove=1
        r=random.uniform(0.5,d_max)
        theta=math.radians(random.uniform(0,360))
        phi=math.radians(random.uniform(0,180))
        direction=[r*math.sin(theta)*math.cos(phi),r*math.sin(theta)*math.sin(phi),r*math.cos(theta)]
        ratmlocnew=[0]*natomsmove
        for i in range(natomsmove):
            ratmloc=random.randint(0,len(positions)-1)
            ratmlocnew[i]=(direction[0]+positions[ratmloc][0],direction[1]+positions[ratmloc][1],direction[2]+positions[ratmloc][2])
            positions[ratmloc]=ratmlocnew[i]
        if Optimizer.structure=='Defect':
            if Optimizer.isolate_mutation:
                indc.set_positions(positions)
                indiv[0]=indc.copy()
                indiv[0].extend(indb)
            else:
                indiv[0].set_positions(positions)
        else:
            indiv[0].set_positions(positions)
    else:
        natomsmove=0
        ratmlocnew=0
    Optimizer.output.write('Lattice Alteration Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write('Number of atoms moved = '+repr(natomsmove)+'\n')
    Optimizer.output.write(repr(ratmlocnew)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='LARD'+repr(natomsmove)	
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv