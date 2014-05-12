from ase import Atom, Atoms
from ase.io import read, write
import math
import numpy

def finite_size_scale(standard, ssize, primordial, fsize, psize=[1,1,1],writefile=True):
    """Function to perform finite size scaling for defect structure relaxation
    **Assumes cell information can be retrieved by reading file**
    Inputs:
        standard = File of structure containing defect
        ssize = Supercell size of structure with defect (list of size 3)
        primordial = File of structure for basic unit of perfect cell
        psize = Supercell size of structure for padding. Default = [1,1,1] (list of size 3)
        fsize = Desired supercell size of final structure (list of size 3)
    Outputs:
        POSCAR file of structure containing defect and padding"""
    
    # Check if the input sizes work out with the desired final size
    padding = [0,0,0]
    srcon = [0,0,0]
    for i in range(3):
        diff = fsize[i] - ssize[i]
        if diff < 0:
            raise RuntimeError('Desired final size of the structure must be larger than existing defect structure size. Defect Size = '+repr(ssize)+' Final Size = '+repr(fsize))
        elif diff >= 0:
             if math.fmod(diff,psize[i]):
                raise RuntimeError('Primordial structure and defect structure sizes cannot be used to form desired final size.  Reduce size of primordial structure. Defect Size = '+repr(ssize)+' Final Size = '+repr(fsize)+' Primordial size = '+repr(psize))
             else:
                padding[i] = diff/psize[i]
    
    # Load the defect structure and primordial structure
    defst = read(standard)
    pst = read(primordial)
    
    # Pad the structure
    positions = pst.get_positions()
    syms = pst.get_chemical_symbols()
    final = defst.copy()
    lv = [one/ssize for one in defst.get_cell()]
    vect = []
    for m0 in range(padding[0]):
        for m1 in numpy.arange(0,fsize[1],psize[1]):
            for m2 in numpy.arange(0,fsize[2],psize[2]):
                vect.append([ssize[0]+m0*psize[0],m1,m2])

    for m1 in range(padding[1]):
        for m0 in numpy.arange(0,ssize[0],psize[0]):
            for m2 in numpy.arange(0,fsize[2],psize[2]):
                vect.append([m0,ssize[1]+m1*psize[1],m2])

    for m2 in range(padding[2]):
        for m0 in numpy.arange(0,ssize[0],psize[0]):
            for m1 in numpy.arange(0,ssize[1],psize[1]):
                vect.append([m0,m1,ssize[2]+m2*psize[2]])

    for m0,m1,m2 in vect:
        npos = positions + numpy.dot((m0, m1, m2), lv)
        for i in range(len(npos)):
            final.append(Atom(symbol=syms[i],position=npos[i]))
    
    final.set_cell(numpy.array([fsize[c] * lv[c] for c in range(3)]))
    
    # Write output as POSCAR
    if writefile:
        write('POSCAR_Final', final)
    
    return final

# Potential issues:
#    --Doesn't account for internal periodicity of relaxed defect structure.  
#        -Could add this by checking distance
#    --Currently uses ASE but should be able to switch to pymatgen without too much issue
#    --Assumes defect structure is aligned from [0,0,0]

if __name__ == "__main__":
    #Sanity checks    
    fil=open('Test_FSS.xyz','a')
    
    #Test 1 - Final size cannot be constructed directly from primordial structure
    from ase.lattice.cubic import BodyCenteredCubic as BCC
    ssize=[3,3,3]
    fsize=[7,3,3]
    psize=[2,3,3]
    defectstructure = BCC('Fe', size = ssize)
    primordialstructure = BCC('Fe', size = psize)
    f = BCC('Fe', size = fsize)
    write('POSCAR_defect1',defectstructure)
    write('POSCAR_prim1',primordialstructure)
    final = finite_size_scale('POSCAR_defect1', ssize, 'POSCAR_prim1', fsize, psize)
    write(fil,f,'xyz')
    write(fil,final,'xyz')
    fil.flush()
    
    #Test 2 - Test in all 3 dimensions
    ssize=[2,2,2]
    fsize=[4,4,4]
    psize=[1,1,1]
    defectstructure = BCC('Fe', size = ssize)
    primordialstructure = BCC('Fe', size = psize)
    f = BCC('Fe', size = fsize)
    write('POSCAR_defect2',defectstructure)
    write('POSCAR_prim2',primordialstructure)
    final = finite_size_scale('POSCAR_defect2', ssize, 'POSCAR_prim2', fsize, psize)
    write(fil,f,'xyz')
    write(fil,final,'xyz')
    fil.flush()
    
    #Test 3 - Test for not cubic lattice
    from ase.lattice.spacegroup import crystal
    a = 3.21
    c = 5.21
    mg = crystal('Mg', [(1./3., 2./3., 3./4.)], spacegroup=194,cellpar=[a, a, c, 90, 90, 120])
    defectstructure = mg.repeat(ssize)
    primordialstructure = mg.repeat(psize)
    f = mg.repeat(fsize)
    write('POSCAR_defect3',defectstructure)
    write('POSCAR_prim3',primordialstructure)
    final = finite_size_scale('POSCAR_defect3', ssize, 'POSCAR_prim3', fsize, psize)
    write(fil,f,'xyz')
    write(fil,final,'xyz')
    fil.flush()
    fil.close()
    
    