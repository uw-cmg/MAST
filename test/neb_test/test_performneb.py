import pymatgen
import os
import shutil
import sys
import unittest
import time
from filecmp import dircmp

from MAST.ingredients.performneb import PerformNEB
from MAST.utility import MASTError
from MAST.utility.dirutil import *

class TestPerformNEB(unittest.TestCase):
    def setUp(self):
        scripts=get_mast_install_path()
        os.chdir(os.path.join(scripts,'test','neb_test'))
        if os.path.exists('nebtest_neb_10-11'):
            self.tearDown()
        os.mkdir('nebtest_neb_10-11')
        shutil.copy('nebtest_defect_10_stat/CONTCAR','nebtest_neb_10-11/parent_structure_10')
        shutil.copy('nebtest_defect_11_stat/CONTCAR','nebtest_neb_10-11/parent_structure_11')
        shutil.copy('nebtest_defect_10_stat/OSZICAR','nebtest_neb_10-11/parent_energy_10')
        shutil.copy('nebtest_defect_11_stat/OSZICAR','nebtest_neb_10-11/parent_energy_11')
        self.NEBing=PerformNEB(name="nebtest_neb_10-11", program="vasp", program_keys={'images':3,'mast_kpoints':[2,2,2,"M"],'mast_xc':"PBE",'mast_setmagmom':"5 1 5",'neblines':["Mg, 0 0 0.1, 0.75 0 0"]})

    def tearDown(self):
        if os.path.exists('nebtest_neb_10-11'):
            shutil.rmtree('nebtest_neb_10-11')
        if os.path.exists('nebtest_neb_10-11_image01_stat'):
            shutil.rmtree('nebtest_neb_10-11_image01_stat')
            shutil.rmtree('nebtest_neb_10-11_image02_stat')
            shutil.rmtree('nebtest_neb_10-11_image03_stat')
    
    def test_unsupported_program(self):
        self.NEBing.keywords['program']='no program'
        self.assertRaises(MASTError, self.NEBing.write_files)

    def test_get_my_labels(self):
        labels = self.NEBing.get_my_labels()
        print labels
        llist=['10','11']
        self.assertEqual(labels,llist)

    def test_get_parent_structures(self):
        structures=self.NEBing.get_parent_structures()
        from pymatgen.io.vaspio import Poscar
        manstr=[]
        manstr.append(Poscar.from_file('nebtest_defect_10_stat/CONTCAR').structure)
        manstr.append(Poscar.from_file('nebtest_defect_11_stat/CONTCAR').structure)
        self.assertEqual(self.NEBing.get_parent_structures(),manstr)

    def test_interpolated_images(self):
        self.NEBing.write_files()
        difflist = dircmp("nebtest_neb_10-11","expected_interpolation_write").diff_files
        #shutil.copytree("nebtest_neb_10-11","nebtest_neb_10-11_save")
        print difflist
        self.assertEqual(len(difflist), 0)

    def test_is_not_complete(self):
        self.NEBing.write_files()
        self.assertEqual(self.NEBing.is_complete(),False)

    def test_is_partially_complete(self):
        self.NEBing.write_files()
        shutil.copy(os.getcwd() + "/copyoutcar01", "nebtest_neb_10-11/01/OUTCAR")
        shutil.copy(os.getcwd() + "/copyoutcar02", "nebtest_neb_10-11/02/OUTCAR")
        time.sleep(1)
        self.assertEqual(self.NEBing.is_complete(),False)

    def test_is_complete(self):
        self.NEBing.write_files()
        self.force_complete()
        self.assertEqual(self.NEBing.is_complete(),True)

    def force_complete(self):
        shutil.copy(os.getcwd() + "/copyoutcar01", "nebtest_neb_10-11/01/OUTCAR")
        shutil.copy(os.getcwd() + "/copyoutcar02", "nebtest_neb_10-11/02/OUTCAR")
        shutil.copy(os.getcwd() + "/copyoutcar03", "nebtest_neb_10-11/03/OUTCAR")
        shutil.copy(os.getcwd() + "/copycontcar01", "nebtest_neb_10-11/01/CONTCAR")
        shutil.copy(os.getcwd() + "/copycontcar02", "nebtest_neb_10-11/02/CONTCAR")
        shutil.copy(os.getcwd() + "/copycontcar03", "nebtest_neb_10-11/03/CONTCAR")
        time.sleep(1)

    def test_update_children(self):
        childname="nebtest_neb_10-11_image"
        self.NEBing.keywords['child_dict']={childname + '01_stat':"",
                childname + '02_stat':"", childname + '03_stat':""}
        for key in self.NEBing.keywords['child_dict'].iterkeys():
            os.makedirs(key)
        self.NEBing.write_files()
        self.force_complete()
        self.NEBing.update_children()
        difflist = dircmp("nebtest_neb_10-11_image01_stat","expected_image01_stat").diff_files
        self.assertEqual(len(difflist), 0)
        difflist = dircmp("nebtest_neb_10-11_image02_stat","expected_image02_stat").diff_files
        self.assertEqual(len(difflist), 0)
        difflist = dircmp("nebtest_neb_10-11_image03_stat","expected_image03_stat").diff_files
        self.assertEqual(len(difflist), 0)
