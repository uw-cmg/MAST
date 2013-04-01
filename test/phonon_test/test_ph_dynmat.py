import pymatgen
import os
import shutil
import sys
import unittest
import time
import filecmp
from filecmp import dircmp

from MAST.ingredients.frozenphonons import FrozenPhonons
from MAST.utility import MASTError
from pymatgen.io.vaspio import Poscar
from MAST.ingredients.pmgextend import vasp_extensions

class TestDynmat(unittest.TestCase):
    def setUp(self):
        scripts = os.getenv("SCRIPTPATH")
        os.chdir(os.path.join(scripts, 'test','phonon_test'))

    def tearDown(self):
        pass

    def test_dynmat(self):
        myPos=Poscar.from_file("dynPOSCAR")
        mydir = os.getcwd()
        hess = vasp_extensions.combine_dynmats(myPos, mydir)
        self.assertRaises(MASTError, echo, "hello")
