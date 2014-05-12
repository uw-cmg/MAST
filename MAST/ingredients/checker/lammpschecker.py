##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Amy Kacmarowski
# Last updated: 2014-04-25
##############################################################
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar
from MAST.utility import dirutil
from MAST.utility.mastfile import MASTFile
from MAST.utility import MASTError
from MAST.ingredients.checker import BaseChecker
from MAST.ingredients.checker import VaspChecker
import os
import logging
import pymatgen
import numpy as np
import time
class LammpsChecker(BaseChecker):
    """LAMMPS checker functions
        Mostly structure functions right now.
    """
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseChecker.__init__(self, allowed_keys, **kwargs)

    def forward_final_structure_file(self, mydir=""):
        """Get the final structure from LAMMPS dump file.
            The dump file must be named "dump.out"
            Args:
                mydir <str>: Directory. If no directory is given,
                    use current ingredient directory given
                    as keyword 'name'
        """
        if mydir == "":
            mydir = self.keywords['name']
        dumpname = "dump.out"
        curdir = self.keywords['name']
        #import data
        #import dump
        #import xyz
        import ase
        from ase.io import read, write
        dumppath = os.path.join(curdir, dumpname)
        if not os.path.isfile(dumppath):
            raise MASTError(self.__class__.__name__,"No dump file at %s" % dumppath)
        #dumpobject = dump.dump(dumppath)
        #dumpobject.map()
        #xyzobject = xyz.xyz(dumpobject)
        #xyzobject.many("temp_many_output")
        #myfiles = os.listdir(curdir)
        #myfiles.sort()
        #xyzfiles=list()
        #for myfile in myfiles:
        #    if "temp_many_output" in myfile:
        #        xyzfiles.append(myfile)
        #lastxyz = os.path.join(mydir, xyzfiles[-1])
        #aseatoms = ase.io.read(lastxyz)
        aseatomsobj = ase.io.read(dumppath)
        newpath = os.path.join(curdir, "CONTCAR")
        ase.io.write(newpath,aseatomsobj,"vasp", direct=True, sort=True, vasp5=True)
        self.copy_a_file(mydir, "CONTCAR", "POSCAR")
