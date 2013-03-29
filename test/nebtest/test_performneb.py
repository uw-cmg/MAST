import pymatgen
import os
import shutil
import sys
import unittest
import time

from MAST.ingredients.performneb import PerformNEB
class TestPerformNEB(unittest.TestCase):
    def setUp(self):
        os.chdir("//home/tam/bin/git/MAST4pymatgen/test/nebtest")
        os.mkdir('nebtest_neb10-11')
        shutil.copy('copyparent1','nebtest_neb10-11/parent_path_10')
        shutil.copy('copyparent2','nebtest_neb10-11/parent_path_11')

    def tearDown(self):
        shutil.rmtree('nebtest_neb10-11')

    def test_unsupported_program(self):
        NEBing=PerformNEB(name="nebtest_neb10-11", program="no program", program_keys={'images':3})
        image_list=NEBing.do_interpolation(NEBing.get_parent_paths())
        self.assertEqual(image_list, None)

    def test_supported_program(self):
        NEBing=PerformNEB(name="nebtest_neb10-11", program="vasp", program_keys={'images':3})
        image_list=NEBing.do_interpolation(NEBing.get_parent_paths())
        self.assertNotEqual(image_list, None)

    def test_get_parent_paths(self):
        NEBing=PerformNEB(name="nebtest_neb10-11", program="vasp", program_keys={'images':3})
        self.assertEqual(NEBing.get_parent_paths(),['ep1','ep2'])

    def test_is_not_complete(self):
        NEBing=PerformNEB(name="nebtest_neb10-11", program="vasp", program_keys={'images':3})
        NEBing.write_files()
        self.assertEqual(NEBing.is_complete(),False)

    def test_is_partially_complete(self):
        NEBing=PerformNEB(name="nebtest_neb10-11", program="vasp", program_keys={'images':3})
        NEBing.write_files()
        shutil.copy(os.getcwd() + "/OUTCAR", "nebtest_neb10-11/01")
        time.sleep(1)
        shutil.copy(os.getcwd() + "/OUTCAR", "nebtest_neb10-11/02")
        time.sleep(1)
        self.assertEqual(NEBing.is_complete(),False)

    def test_is_complete(self):
        NEBing=PerformNEB(name="nebtest_neb10-11", program="vasp", program_keys={'images':3})
        NEBing.write_files()
        shutil.copy(os.getcwd() + "/OUTCAR", "nebtest_neb10-11/01")
        time.sleep(1)
        shutil.copy(os.getcwd() + "/OUTCAR", "nebtest_neb10-11/02")
        time.sleep(1)
        shutil.copy(os.getcwd() + "/OUTCAR", "nebtest_neb10-11/03")
        time.sleep(1)
        self.assertEqual(NEBing.is_complete(),True)



