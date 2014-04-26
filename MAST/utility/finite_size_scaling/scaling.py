#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Wei Xie
# Last updated: 2014-04-25
##############################################################
import sys, getopt
from numpy import *
import pymatgen as mg
import shutil
from pymatgen.core.structure import Structure


def read(inp):
    my_struct = {}
    struct = mg.read_structure(inp)
    for i in range(len(struct)):
        try: my_struct[str(struct[i].specie)]
        except: my_struct[str(struct[i].specie)] = [] # New element specie added
        my_struct[str(struct[i].specie)].append([struct[i].a,struct[i].b,struct[i].c])
            
    return my_struct

def make_super(size,filename):
    shutil.copyfile('POSCAR_primitive', filename)
    struct = mg.read_structure(filename)
    Structure.make_supercell(struct,size)
    mg.write_structure(struct,filename)
        
def corrctReflection(relax,defect): 
    # Correct periodically reflected atoms according to the perfect cell.
    # relax is the CONTCAR (after relaxation) and defect is the POSCAR (before relaxation).
    shift = array([0,0,0])
    for ele in relax.keys():
        for i in range(len(relax[ele])):
            for j in range(3):
                shift[j] = shift[j] + relax[ele][i][j] - defect[ele][i][j]
    shift = shift/len(mg.read_structure('POSCAR_defect'))
            
    for ele in relax.keys(): 
        for i in range(len(relax[ele])):
            for j in range(3):
                relax[ele][i][j]=relax[ele][i][j]-shift[j]
                if abs(relax[ele][i][j]-defect[ele][i][j])>abs(relax[ele][i][j]+defect[ele][i][j]-1):
                    if abs(relax[ele][i][j]-defect[ele][i][j]-1)<abs(relax[ele][i][j]-defect[ele][i][j]+1):
                        relax[ele][i][j] = relax[ele][i][j] - 1
                    else:
                        relax[ele][i][j] = relax[ele][i][j] + 1
    return relax

def rescale(struct,size,Size): 
    # Necessary for fractional coordinates.
    for ele in struct.keys():
        for i in range(len(struct[ele])):
            for k in range(3):
                struct[ele][i][k] = struct[ele][i][k]*size[k]/Size[k]
    return struct

def delete(supercell,chunk):
    for ele in supercell.keys():
        total=len(supercell[ele])
        i=0
        while i<total:
            for n in range(len(chunk[ele])):
                if abs(supercell[ele][i][0] - chunk[ele][n][0])<1e-6 and abs(supercell[ele][i][1] - chunk[ele][n][1])<1e-6 and abs(supercell[ele][i][2] - chunk[ele][n][2])<1e-6:
                    supercell[ele].remove(supercell[ele][i]) # Removing the small cell from supercell.  
                    i = i - 1            
                    total=total-1
            i = i + 1
                                                 
    return supercell

def embed(supercell,relax):
    for ele in relax.keys():
        for i in range(len(relax[ele])):
            try: supercell[ele]
            except:  supercell[ele] = []
            supercell[ele].append(relax[ele][i])
    return supercell

def write(out,supercell):
    fp = open(out,'w')
    fp.write('Final Structure\n')
    i = 0
    for line in open('POSCAR_super','r'):
        if i>0 and i<=4: fp.write(line)
        i = i + 1
    for ele in supercell.keys():
        fp.write(ele + ' ')
    fp.write('\n')
    for ele in supercell.keys():
        fp.write(str(len(supercell[ele])) + ' ')
    fp.write('\n')
    fp.write('Direct\n')
    for ele in supercell.keys():
        for i in range(len(supercell[ele])):
            fp.write(str(supercell[ele][i][0])+' '+str(supercell[ele][i][1])+' '+str(supercell[ele][i][2])+'\n')
    fp.close()




#this section reads command line args
argv = sys.argv[1:]
opts, args = getopt.getopt(argv,"s")
target_size = [1,1,1]
primordial_size = [1,1,1]
for i in range(3):
    target_size[i] = int(args[i])
    primordial_size[i] = int(args[i+3])

make_super(target_size,'POSCAR_super')
supercell = read('POSCAR_super')

make_super(primordial_size,'POSCAR_chunk')
chunk = read('POSCAR_chunk')
chunk = rescale(chunk,primordial_size,target_size)

supercell = delete(supercell,chunk)

relax = read('POSCAR_relax')
defect = read('POSCAR_defect')
relax = corrctReflection(relax,defect)
relax = rescale(relax,primordial_size,target_size)

supercell = embed(supercell,relax)
write('POSCAR_final',supercell)


