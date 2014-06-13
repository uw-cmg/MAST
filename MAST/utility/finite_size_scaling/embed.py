#!/usr/bin/env python
#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Wei Xie
# Last updated: 2014-04-21
##############################################################
import sys
import getopt 
import os
import shutil
import errno
import numpy as np
from scipy import stats

import pymatgen as mg
from pymatgen.analysis import ewald
from pymatgen.io import vaspio

def CenterMassnReflectionCorr(relax,defect): 
    # Correct periodically reflected atoms according to the perfect cell.
    # relax is the CONTCAR (after relaxation) and defect is the POSCAR (before relaxation).
    shift = np.array([0,0,0])
    for ele in relax.keys():
        for i in range(len(relax[ele])):
            for j in range(3):
                shift[j] = shift[j] + relax[ele][i][j] - defect[ele][i][j]
    #shift = shift/len(mg.read_structure('POSCAR_defect'))
    shift = shift/relax.composition.num_atoms 
            
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
    
def rescale(struct,size,Size): 
    # Necessary for fractional coordinates.
    for ele in struct.keys():
        for i in range(len(struct[ele])):
            for k in range(3):
                struct[ele][i][k] = struct[ele][i][k]*size[k]/Size[k]
    return struct
