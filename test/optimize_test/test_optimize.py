import pymatgen
import os
import shutil
import sys
import unittest
import time
import filecmp
from filecmp import dircmp

from MAST.ingredients.optimize import Optimize
from MAST.utility import MASTError
from pymatgen.io.vaspio import Poscar
from MAST.utility.dirutil import *

class TestOptimize(unittest.TestCase):
    def setUp(self):
        scripts = get_mast_install_path()
        os.chdir(os.path.join(scripts, 'test','optimize_test'))
        os.mkdir('test_opt1')

    def tearDown(self):
        shutil.rmtree('test_opt1')
        try:
            shutil.rmtree('test_opt2')
        except OSError:
            pass
    
    def test_unsupported_program(self):
        myOpt=Optimize(name="test_opt1", program="no program", program_keys={'isif':3})
        self.assertRaises(MASTError, myOpt.write_files)

    def test_write_files_opt1(self):
        myOpt=Optimize(name="test_opt1", program="vasp", program_keys={'ibrion':2})
        myOpt.keywords['structure'] = Poscar.from_file('expected_opt1/POSCAR').structure 
        myOpt.write_files()
        difflist = dircmp("test_opt1","expected_write").diff_files
        #shutil.copytree("test_opt1","test_opt1_save")
        print difflist
        self.assertEqual(len(difflist), 0)

    def test_is_not_complete(self):
        myOpt=Optimize(name="test_opt1", program="vasp", program_keys={'ibrion':2})
        self.assertEqual(myOpt.is_complete(),False)

    def test_is_complete(self):
        myOpt=Optimize(name="test_opt1", program="vasp", program_keys={'ibrion':2})
        shutil.copy("expected_opt1/OUTCAR","test_opt1")
        time.sleep(1)
        self.assertEqual(myOpt.is_complete(),True)

    def test_is_ready(self):
        myOpt=Optimize(name="test_opt1", program="vasp", program_keys={'ibrion':2})
        self.assertEqual(myOpt.is_ready_to_run(),False)
        myOpt.keywords['structure'] = Poscar.from_file('expected_opt1/POSCAR').structure 
        myOpt.write_files()
        self.assertEqual(myOpt.is_ready_to_run(),True)

    def test_update_children(self):
        os.mkdir("test_opt2")
        myOpt=Optimize(name="test_opt1", program="vasp", program_keys={'ibrion':2},child_dict={'test_opt2':""})
        shutil.copy("expected_opt1/OUTCAR","test_opt1")
        shutil.copy("expected_opt1/CONTCAR","test_opt1")
        myOpt.update_children()
        self.assertTrue(filecmp.cmp("test_opt1/CONTCAR","test_opt2/POSCAR"))
