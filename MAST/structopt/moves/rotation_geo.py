import random
import numpy
from ase import Atom, Atoms
from MAST.structopt.tools.find_defects import find_defects

def rotation_geo(indiv, Optimizer):
    '''Function to handle a geometry based rotation mutation
    Rotates a small group of atoms based on location within structure
    Input:
        indiv = structopt Individual class object to be mutated
        Optimizer = structopt Optimizer class object
    Output:
        indiv = structopt Individual class object that has been mutated
    '''
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False
    #Rotate group of atoms based on geometry
    if Optimizer.structure=='Defect':
        if Optimizer.isolate_mutation:
            atms,indb,vacant,swap,stro = find_defects(indiv[0],Optimizer.solidbulk,0)
        else:
            atms = indiv[0].copy()
    else:
        atms=indiv[0].copy()
    nat=len(atms)
    if nat > 1:
        #Select number of atoms to rotate
        if nat==2:
            natrot=2
        else:
            natrot = random.randint(2,nat)
        #Select random position in cluster
        cell=numpy.maximum.reduce(atms.get_positions())-numpy.minimum.reduce(atms.get_positions())
        pt=[random.uniform(0,cell[0]),random.uniform(0,cell[1]),random.uniform(0,cell[2])]
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
        for i in range(natrot):
            atmsr.append(dlist[i][1])
            indexlist.append(dlist[i][1].index)
        #Identify random rotation and rotate selected atoms
        ax=['x', '-x','y','-y','z','-z']
        rax=ax[random.randint(0,len(ax)-1)]
        rang=random.uniform(30,180)
        #rang=random.random()*90
        atmsr.rotate(rax,a=rang,center='COM',rotate_cell=False)
        #Update individual with new atom positions
        for i in range(len(indexlist)):
            index=indexlist[i]
            atms[index].position=atmsr[i].position
        if Optimizer.structure=='Defect':
            if Optimizer.isolate_mutation:
                atms.extend(indb)
        indiv[0]=atms.copy()
        #Optimizer.output.write('New positions = '+repr(indiv[0].get_positions())+'\n')
    else:
        natrot=0
        pt=0
        rax=0
        rang=0
    Optimizer.output.write('Geometry Rotation Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    Optimizer.output.write('Number of atoms rotated = '+repr(natrot)+'\n')
    Optimizer.output.write('Geometry point = '+repr(pt)+'\n')
    Optimizer.output.write('Rotation vector = '+repr(rax)+'\n')
    Optimizer.output.write('Rotation angle = '+repr(rang)+'\n')
    Optimizer.output.write(repr(indiv[0])+'\n')
    muttype='GR'+repr(natrot)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv