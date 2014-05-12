import random
import numpy
from ase import Atom, Atoms

def gen_pop_sphere(atomlist,size,crystal=False):
    """Function to generate a random structure of atoms within a sphere of given size.
    Inputs:
        atomlist = List of tuples with structure of atoms and quantity
            [('Sym1',int(concentration1), float(mass1),float(chempotential1)),
            ('Sym2',int(concentration2), float(mass2),float(chempotential2)),...]
        size = Float of length of side of cube within which to generate atoms
        crystal = False/List of crystal cell shape options
            list('cubic','orthorhombic','tetragonal','hexagonal','monoclinic','triclinic')
            cell shape will be adjusted accordingly
    Outputs:
        Returns individual of class Atoms (see ase manual for info on Atoms class)
        and if crystal list provided also outputs combined string with output information
    """
    size = float(size)
    indiv=Atoms()
    indiv.set_cell([size,size,size])
    # Get list of atom types for all atoms in cluster
    for s,c,m,u in atomlist:
        if c > 0:
            for i in range(c):
                r = random.random()*size/2.0
                d = [random.uniform(-1.0,1.0) for j in range(3)]
                u = [d[j]/numpy.linalg.norm(d) for j in range(3)]
                pos = [r*u[j]+size/2.0 for j in range(3)]
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

