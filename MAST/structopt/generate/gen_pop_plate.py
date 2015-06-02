import random
import numpy
import math
from ase import Atom, Atoms
#from MAST.structopt.generate import __init__ 
from MAST.structopt.generate.rot_vec import * 

#def gen_pop_plate(atomlist, size, crystal=False):
def gen_pop_plate(atomlist, sizeinput, dir_par=[1,-1,0], dir_per=[1,1,1], crystal=False):
    '''
     Last edited by Hyunseok Ko on 10-15-2014
     Function to generate defects in a plate. The center of plate is equal to the center of mass of the bulk structure. 
     -- Plate orientation can be specified by following parameters::
     		dir_par and dir_per is a directional vector (do not necessarily be a unit vector) which sets the axis of the plate
     -- Plate dimension can be chosen (Needs development. Manual setup for now) by setting 'size' in this module. 
    '''
# Setup axis vectors for plate 
    size=[sizeinput, sizeinput, 3]
    #print 'HKK :: plate vector', size
    # Set unit vectors for plate
    dir_crs=numpy.cross(dir_par, dir_per)
    uvec_par=[x/((dir_par[0]**2+dir_par[1]**2+dir_par[2]**2)**0.5) for x in dir_par]
    uvec_crs=[x/((dir_crs[0]**2+dir_crs[1]**2+dir_crs[2]**2)**0.5) for x in dir_crs]
    uvec_per=[x/((dir_per[0]**2+dir_per[1]**2+dir_per[2]**2)**0.5) for x in dir_per]
    #print 'Unit vectors of new xyz axis'
    #print uvec_par, uvec_crs, uvec_per
    
    # Set vectors for plate
    vec_par=[x*size[0] for x in uvec_par]
    vec_crs=[x*size[1] for x in uvec_crs]
    vec_per=[x*size[2] for x in uvec_per]
    # Compute rotation angles from primary axis to axis for plate (If loop to avoid floating error)
    temp=[-uvec_per[1]/((1-uvec_per[2]**2)**0.5), uvec_per[2],uvec_crs[2]/((1-uvec_per[2]**2)**0.5)]
    if [temp[x] for x in range(3)]  >=1:
        temp[x]=1
    alpha=numpy.arccos(temp[0])
    beta=numpy.arccos(temp[1])
    gamma=numpy.arccos(temp[2])
    
# Get list of atom types for all atoms in cluster	
    indiv=Atoms()
    indiv.set_cell(size)

    for s,c,m,u in atomlist:
        if c > 0:
            for i in range(c):
                temppos = [random.uniform(0,size[j]) for j in range(3)]
                # Rotating the defect position to the plate 
                pos = rot_vec(temppos,alpha,beta,gamma)
                at=Atom(symbol=s,position=pos)
                indiv.append(at)
    if crystal:
        stro=''
        natoms=sum([c for s,c,m,u in atomlist])
        pos=indiv.get_scaled_positions()
        structure=random.choice(crystal)
        cello=indiv.get_cell()
        if structure=='cubic':
            #Set to cubic shape
            an,bn,cn = [numpy.linalg.norm(v) for v in cello]
            a=(an+bn+cn)/3.0
            celln=numpy.array([[a,0,0],[0,a,0],[0,0,a]])
            stro+='Setting cell to cubic\n'
        elif structure=='orthorhombic':
            #Set to orthorhombic
            a=random.uniform(2,natoms**0.3333*size)
            b=random.uniform(2,natoms**0.3333*size)
            c=random.uniform(2,natoms**0.3333*size)
            celln=numpy.array([[a,0,0],[0,b,0],[0,0,c]])
            stro+='Setting cell to orthorhombic\n'
        elif structure=='tetragonal':
            #Set to tetragonal shape
            an,bn,cn = [numpy.linalg.norm(v) for v in cello]
            a=(an+bn)/2.0
            c=cn
            if c==a:
                c=random.uniform(1,natoms**0.3333*size)
            celln=numpy.array([[a,0,0],[0,a,0],[0,0,c]])
            stro+='Setting cell to tetragonal\n'
        elif structure=='hexagonal':
            #Set to hexagonal shape
            an,bn,cn = [numpy.linalg.norm(v) for v in cello]
            a=(an+bn)/2.0
            c=cn
            if c<=a:
                c=random.uniform(a+1,natoms**0.3333*size)
            trans=numpy.array([[1,0,0],[-0.5,(3.0**0.5)/2.0,0],[0,0,1]])
            trans[0]=[a*i for i in trans[0]]
            trans[1]=[a*i for i in trans[1]]
            trans[2]=[c*i for i in trans[2]]
            celln=trans
            stro+='Setting cell to Hexagonal\n'
        elif structure=='monoclinic':
            #Set to monoclinic
            a,b,c = [numpy.linalg.norm(v) for v in cello]
            if a==b:
                b=random.uniform(1,natoms**0.3333*size)
            trans=numpy.array([(1+random.random())*c, 0, (1+random.random())*c])
            celln=numpy.array([[a,0,0],[0,b,0],[0,0,0]])
            celln[2]=trans
            stro+='Setting cell to monoclinic\n'
        elif structure=='triclinic':
            #Set to triclinic
            a,b,c = [numpy.linalg.norm(v) for v in cello]
            celln=numpy.array([[a,0,0],[(1+random.random())*b,(1+random.random())*b,0],[(1+random.random())*c,0,(1+random.random())*c]])
            stro+='Setting cell to triclinic\n'
        indiv.set_cell(celln)
        indiv.set_scaled_positions(pos)
        stro+=repr(indiv.get_cell())+'\n'
        return indiv, stro
    return indiv

