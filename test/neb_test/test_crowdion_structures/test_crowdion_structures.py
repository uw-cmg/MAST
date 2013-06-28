import pymatgen
import MAST
import os
import shutil
import sys
import unittest
import time
from filecmp import dircmp

from MAST.ingredients.neb import NEB
from MAST.utility import MASTError
from MAST.utility.dirutil import *



class TestCrowdionStructures(unittest.TestCase):
    def setUp(self):
        scripts=get_mast_install_path()
        os.chdir(os.path.join(scripts,'test','neb_test',
                'test_crowdion_structures')) 
        if os.path.exists('output'):
            self.tearDown()
        self.progkeys = dict()
        self.progkeys['neblines'] = dict()


    def tearDown(self):
        if os.path.exists('output'):
            shutil.rmtree('output')

    def test_structure_sorting(self):
        def1str = pymatgen.io.smartio.read_structure("POSCAR_defectgroup1")
        def2str = pymatgen.io.smartio.read_structure("POSCAR_defectgroup2")
        self.progkeys['neblines']["1-2"] = list()
        self.progkeys['neblines']["1-2"].append(["Cr","0.3 0 0", "0 0 0"])
        self.progkeys['neblines']["1-2"].append(["Ni", "0.6 0 0", "0.3 0 0"])
        myneb=NEB(name="neb_1-2", program_keys = self.progkeys)
        myneb.keywords['structure'] = pymatgen.io.vaspio.Poscar.from_file('POSCAR_bulk').structure
        sorted1=myneb.sort_structure_and_neb_lines(def1str,0)
        sorted2=myneb.sort_structure_and_neb_lines(def2str,1)
        os.mkdir('output')
        pymatgen.io.smartio.write_structure(sorted1,"output/POSCAR_sorted1")
        pymatgen.io.smartio.write_structure(sorted2,"output/POSCAR_sorted2")
        mydiff = dircmp('output','expected').diff_files
        self.assertEquals(len(mydiff),0)
