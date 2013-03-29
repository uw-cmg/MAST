import pymatgen
import os
import shutil
import sys
import unittest
import time
from filecmp import dircmp

from MAST.ingredients.performneb import PerformNEB
class TestPerformNEB(unittest.TestCase):
    def setUp(self):
        os.chdir("//home/tam/bin/git/MAST4pymatgen/test/nebtest")
        os.mkdir('nebtest_neb10-11')
        shutil.copy('parent_path_10','nebtest_neb10-11')
        shutil.copy('parent_path_11','nebtest_neb10-11')

    def tearDown(self):
        shutil.rmtree('nebtest_neb10-11')
        try:
            shutil.rmtree('nebtest_neb10-11_image01_stat')
            shutil.rmtree('nebtest_neb10-11_image02_stat')
            shutil.rmtree('nebtest_neb10-11_image03_stat')
        except OSError:
            pass

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
        self.assertEqual(NEBing.get_parent_paths(),['nebtest_defect10_stat','nebtest_defect11_stat'])

    def test_interpolated_images(self):
        NEBing=PerformNEB(name="nebtest_neb10-11", program="vasp", program_keys={'images':3})
        NEBing.write_files()
        difflist = dircmp("nebtest_neb10-11","expected_interpolation_write").diff_files
        print difflist
        self.assertEqual(len(difflist), 0)

    def test_is_not_complete(self):
        NEBing=PerformNEB(name="nebtest_neb10-11", program="vasp", program_keys={'images':3})
        NEBing.write_files()
        self.assertEqual(NEBing.is_complete(),False)

    def test_is_partially_complete(self):
        NEBing=PerformNEB(name="nebtest_neb10-11", program="vasp", program_keys={'images':3})
        NEBing.write_files()
        shutil.copy(os.getcwd() + "/copyoutcar1", "nebtest_neb10-11/01/OUTCAR")
        time.sleep(1)
        shutil.copy(os.getcwd() + "/copyoutcar2", "nebtest_neb10-11/02/OUTCAR")
        time.sleep(1)
        self.assertEqual(NEBing.is_complete(),False)

    def test_is_complete(self):
        NEBing=PerformNEB(name="nebtest_neb10-11", program="vasp", program_keys={'images':3})
        NEBing.write_files()
        self.force_complete()
        self.assertEqual(NEBing.is_complete(),True)

    def force_complete(self):
        shutil.copy(os.getcwd() + "/copyoutcar01", "nebtest_neb10-11/01/OUTCAR")
        shutil.copy(os.getcwd() + "/copyoutcar02", "nebtest_neb10-11/02/OUTCAR")
        shutil.copy(os.getcwd() + "/copyoutcar03", "nebtest_neb10-11/03/OUTCAR")
        time.sleep(1)
        shutil.copy(os.getcwd() + "/copycontcar01", "nebtest_neb10-11/01/CONTCAR")
        shutil.copy(os.getcwd() + "/copycontcar02", "nebtest_neb10-11/02/CONTCAR")
        shutil.copy(os.getcwd() + "/copycontcar03", "nebtest_neb10-11/03/CONTCAR")
        time.sleep(1)

    def test_update_children(self):
        childname="nebtest_neb10-11_image"
        NEBing=PerformNEB(name="nebtest_neb10-11", program="vasp", program_keys={'images':3}, child_dict={childname+"01_stat":"",childname+"02_stat":"",childname+"03_stat":""})
        for key in NEBing.keywords['child_dict'].iterkeys():
            os.makedirs(key)
        NEBing.write_files()
        self.force_complete()
        NEBing.update_children()
        difflist = dircmp("nebtest_neb10-11_image01_stat","expected_image01_stat").diff_files
        self.assertEqual(len(difflist), 0)
        difflist = dircmp("nebtest_neb10-11_image02_stat","expected_image02_stat").diff_files
        self.assertEqual(len(difflist), 0)
        difflist = dircmp("nebtest_neb10-11_image03_stat","expected_image03_stat").diff_files
        self.assertEqual(len(difflist), 0)
