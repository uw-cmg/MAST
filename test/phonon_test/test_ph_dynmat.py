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
from MAST.utility.dirutil import *

class TestDynmat(unittest.TestCase):
    def setUp(self):
        scripts = get_mast_install_path()
        os.chdir(os.path.join(scripts, 'test','phonon_test'))

    def tearDown(self):
        pass

    def test_dynmat(self):
        myPos=Poscar.from_file("dynPOSCAR")
        mydir = os.getcwd()
        freqs = vasp_extensions.combine_dynmats(myPos, mydir)
        myfreqs=[]
        for freq in freqs:
            myfreqs.append(str(freq))
        print myfreqs
        testfreqs=['5.41762305941e-08', '8.78201163933e-08', '4.0288407072', '4.29528709031', '5.51344211515', '7.23134055466', '7.74865476988', '7.89305465977', '9.36609281415', '10.6187089182', '10.6301691736', 'nan']
        self.assertEqual(myfreqs, testfreqs)
