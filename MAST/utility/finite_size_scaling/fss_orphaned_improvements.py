#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Wei Xie
# Last updated: 2014-02-07
##############################################################
import sys, getopt, os
import shutil
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

from pymatgen import Lattice, Structure
from pymatgen.io.smartio import read_structure, write_structure

stdX = read_structure("POSCAR_std")
prmX = read_structure("POSCAR_prm")
stdL = Structure.lattice(stdX)
prmL = Structure.lattice(stdX)


  
def finite_size_scale(standard, ssize, primordial, fsize, psize=[1,1,1]):
    """Function to perform finite size scaling for defect structure relaxation
    Inputs:
        standard = POSCAR file of structure containing defect
        ssize = Supercell size of structure with defect (list of size 3)
        primordial = POSCAR file of structure for basic unit of perfect cell
        psize = Supercell size of structure for padding. Default = [1,1,1] (list of size 3)
        fsize = Desired supercell size of final structure (list of size 3)
    Outputs:
        POSCAR file of structure containing defect and padding"""
    
    # Check if the input sizes work out with the desired final size
    padding = [0,0,0]
    for i in range(3):
        diff = fsize[i] - ssize[i]
        if diff < 0:
            raise RuntimeError('Desired final size of the structure must be larger than \
existing defect structure size. Defect Size = '+repr(ssize)+' Final Size = '+repr(fsize))
        elif diff >= 0:
             if math.fmod(diff,psize[i]):
                raise RuntimeError('Primordial structure and defect structure sizes cannot \
be used to form desired final size.  Reduce size of primordial structure. Defect Size = '+
    repr(ssize)+' Final Size = '+repr(fsize)+' Primordial size = '+repr(psize))
             else:
                padding[i] = diff/psize[i]
    
    # Load the defect structure and primordial structure
    try:
        defst = read_structure(standard)
    except:
        raise RuntimeError('Error: Unable to read standard structure.  Please check file. Filename: '+\
            standard)
    try:
        pst = read_structure(primordial)
    except:
        raise RuntimeError('Error: Unable to read primordial structure.  Please check file. Filename: '+\
            primordial)
    
    # Pad the structure
    positions = [site.coords for site in pst]
    syms = [str(site.specie.symbol) for site in pst]
    lv = [one/ssize for one in defst.lattice.matrix]
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
    
    #Construct a new structure with desired size
    new_lat = Lattice(numpy.array([fsize[c] * lv[c] for c in range(3)]))
    final = Structure(new_lat, defst.species_and_occu,defst.cart_coords,
            coords_are_cartesian=True)
    for m0,m1,m2 in vect:
        npos = positions + numpy.dot((m0, m1, m2), lv)
        for i in range(len(npos)):
            final.append(syms[i],npos[i],coords_are_cartesian=True)
    
    #Check for periodic issues in final structure
    final = check_periodic(final,defst)
    
    # Write output as POSCAR
    write_structure(final, 'POSCAR_Final')
    
    return final

def check_periodic(final,defect):
    # Identify nearest neighbor atoms for each atom in perfect structure
    #Identify trouble atoms
    natoms = len(defect)
    issueatoms = []
    for i in range(natoms):
        dist = [j for j in range(natoms,len(final)) 
            if i != j and final[i].distance(final[j])<=0.6]
        if len(dist) !=0:
            issueatoms.append(dist[0])
    #Idenfity nearest boundary and new coordinates
    coords = final.cart_coords
    for i in issueatoms:
        xflag = False
        yflag = False
        zflag = False
        scaledpos = numpy.linalg.solve(defect.lattice.matrix.T,final.cart_coords[i].T)
        relative = 0.1
        #Try reflection across x axis
        if 1-relative*0.1 < scaledpos[0] < 1+relative:
            npos = final.cart_coords[i]-numpy.dot(numpy.array([1,0,0]),defect.lattice.matrix)
            xflag = True
            coords[i][0] = npos[0]
        #Try reflection across y axis
        if 1-relative*0.1 < scaledpos[1] < 1+relative:
            npos = final.cart_coords[i]-numpy.dot(numpy.array([0,1,0]),defect.lattice.matrix)
            yflag = True
            coords[i][1] = npos[1]
        #Try reflection across z axis
        if 1-relative*0.1 < scaledpos[2] < 1+relative:
            npos = final.cart_coords[i]-numpy.dot(numpy.array([0,0,1]),defect.lattice.matrix)
            zflag = True
            coords[i][2] = npos[2]
        #print xflag, yflag, zflag
    #Generate new structure based on coordinates
    nfinal = Structure(final.lattice.matrix, final.species_and_occu, coords,
                coords_are_cartesian=True)
    return nfinal

# Potential issues:
#	--Poor accounting of internal periodicity  
#		-definition of relative term can be ambiguous depending on structure
#	--Assumes defect structure and primordial structures are aligned

if __name__ == "__main__":
    #Sanity checks	
    #Test 1 - Final size cannot be constructed directly from primordial structure
    ssize=[5,5,5]
    fsize=[9,5,5]
    psize=[2,5,5]
    felat = Lattice.cubic(2.87)
    primf = Structure(felat,["Fe","Fe"],[[0,0,0],[0.5,0.5,0.5]])
    primordialstructure = primf.copy()
    primordialstructure.make_supercell(psize)
    f = primf.copy()
    f.make_supercell(fsize)
    write_structure(primordialstructure,'POSCAR_prim1')
    final = finite_size_scale('POSCAR_Fdefect', ssize, 'POSCAR_prim1', fsize, psize)
    write_structure(f,'POSCAR_F1expected')
    write_structure(final,'POSCAR_F1')
    
    #Test 2 - Test in all 3 dimensions
    ssize=[5,5,5]
    fsize=[8,8,8]
    psize=[1,1,1]
    primordialstructure = primf.copy()
    primordialstructure.make_supercell(psize)
    f = primf.copy()
    f.make_supercell(fsize)
    write_structure(primordialstructure, 'POSCAR_prim2')
    final = finite_size_scale('POSCAR_Fdefect', ssize, 'POSCAR_prim2', fsize, psize)
    write_structure(f,'POSCAR_F2expected')
    write_structure(final,'POSCAR_F2')
    
    #Test 3 - Test for non-cubic lattice
    ssize = [4,4,4]
    fsize=[8,8,8]
    psize=[1,1,1]
    mglat = Lattice.hexagonal(3.21,5.21)
    primmg = Structure(mglat,["Mg","Mg"],[[1./3.,2./3.,3./4.],[2./3.,1./3.,1./4.]])
    primordialstructure = primmg.copy()
    primordialstructure.make_supercell(psize)
    f = primmg.copy()
    f.make_supercell(fsize)
    write_structure(primordialstructure,'POSCAR_prim3')
    final = finite_size_scale('POSCAR_Mdefect', ssize, 'POSCAR_prim3', fsize, psize)
    write_structure(f,'POSCAR_Mg3expected')
    write_structure(final,'POSCAR_Mg3')
