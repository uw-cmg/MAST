import pymatgen
import os
import shutil
import sys
import unittest
import time
import filecmp
from filecmp import dircmp

from pymatgen.io.vaspio import *
from MAST.ingredients.pmgextend import vasp_extensions
from MAST.utility.dirutil import *

class TestPmgextend(unittest.TestCase):
    def setUp(self):
        scripts = get_mast_install_path()
        os.chdir(os.path.join(scripts, 'test','pmgextend_test'))
        self.myposcar = Poscar.from_file("test_poscar")
        self.mypotcar = Potcar(symbols=self.myposcar.site_symbols,
            functional="PBE", sym_potcar_map=None)

    def tearDown(self):
        pass
    
    def test_enmax(self):
        maxenmax = vasp_extensions.get_max_enmax_from_potcar(self.mypotcar)
        self.assertEquals(maxenmax,400)

    def test_nelectrons(self):
        elecs = vasp_extensions.get_total_electrons(self.myposcar, 
            self.mypotcar)
        self.assertEquals(elecs,24)
