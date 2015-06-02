#!/usr/bin/env python
import sys, getopt
from numpy import *
import pymatgen as mg
from pymatgen.core.structure import Structure


def get_shift_vector(center,size,Size):
    shift_vector = array([0.0, 0.0, 0.0])
    for i in range(3):
        if abs(center[0]-0.5*size/Size)>0.25*size/Size:
            shift_vector[i] = 0.5*size/Size - center[i]
    return shift_vector

def shift_coords(struct,shift_vector):
    for i in range(len(struct)):
        coords = struct[i].frac_coords
        for j in range(3):
            if coords[j]+shift_vector[j]>=1:
                coords[j] = coords[j]+shift_vector[j]-1
            elif coords[j]+shift_vector[j]<0:
                coords[j] = coords[j]+shift_vector[j]+1
            else: coords[j] = coords[j]+shift_vector[j]
        struct.replace(i,struct[i].specie,coords)

def find_atom(atom,struct):
    index = None
    struct.append('X',atom)
    for i in range(len(struct)-1):
        if struct.get_distance(i,-1)<1.0: 
            index = i
            break
    struct.__delitem__(-1)
    return index

def compare_coords(relax,defect,supercell,size,Size): # Correct the positions of the relaxed structure according to the defect cell.
    for i in range(len(relax)): # relax is the CONTCAR and defect is the POSCAR.
        r_coords = relax[i].frac_coords
        d_coords = defect[i].frac_coords
        for j in range(3):
            if r_coords[j]-d_coords[j]>0.5:
                r_coords[j] = r_coords[j] - 1
            elif d_coords[j]-r_coords[j]>0.5:
                r_coords[j] = r_coords[j] + 1
        defect.append('X',r_coords)
        if defect.get_distance(i,-1)>1.0:    
            print "Warning: The %sth atom in CONTCAR deviates too much from the POSCAR"%i
            index = find_atom(d_coords*size/Size,supercell)
            supercell.replace(index,supercell[index].specie,r_coords*size/Size)
        defect.__delitem__(-1)
        #relax.replace(i,relax[i].specie,r_coords)

def rescale(struct,size,Size): # Necessary for fractional coordinates.
    for i in range(len(struct)):
        coords = struct[i].frac_coords
        for j in range(3):
            coords[j] = coords[j]*size/Size
        struct.replace(i,struct[i].specie,coords)

def replace(supercell,relax):
    for i in range(len(relax)):
        r_coords = relax[i].frac_coords
        index = find_atom(r_coords,supercell)
        supercell.replace(index,supercell[index].specie,r_coords)

if __name__=='__main__':
    center = [0.0, 0.0, 0.0]
    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv,"s")
    Size = int(args[0])
    size = int(args[1])

    shift_vector = get_shift_vector(center,size,size)
    Shift_Vector = get_shift_vector(center,size,Size)

    supercell = mg.read_structure('POSCAR_super')
    shift_coords(supercell,Shift_Vector)

    relax = mg.read_structure('POSCAR_relax')
    shift_coords(relax,shift_vector)

    defect = mg.read_structure('POSCAR_defect')
    shift_coords(defect,shift_vector)

    compare_coords(relax,defect,supercell,size,Size)
    rescale(relax,size,Size)
    
    replace(supercell,relax)
    shift_coords(supercell,-Shift_Vector)

    mg.write_structure(supercell,'POSCAR_final')

