import os
import time
import unittest
from unittest import SkipTest
import filecmp
from filecmp import dircmp
import MAST
from MAST.ingredients.checker.vaspchecker import VaspChecker
import shutil
import pymatgen
import numpy as np
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import MASTFile
testname ="workflow_test"
testdir = dirutil.get_test_dir(testname)
import subprocess
from MAST.test.workflow_test import workflow_setup

class TestWorkflows(unittest.TestCase):
    """Test Workflows
    """
    def setUp(self):
        return
    def tearDown(self):
        testlist=list()
        testlist.append("simple_optimization.inp")
        for testname in testlist:
            shortname = testname.split(".")[0]
            try:
                pass
                #os.remove("output_%s" % shortname)
                #os.remove("submit_%s.sh" % shortname)
            except OSError:
                pass
        myfiles = os.listdir(testdir)
        for myfile in myfiles:
            if "output_workflow_testing" in myfile:
                pass
                #os.remove(myfile)
        return

    def test_simple_optimization(self):
        [mystatus, my_test_dir]=workflow_setup.generic_submit("simple_optimization.inp")
        if mystatus == "Unfinished":
            self.assertTrue(False)
            return
        elif mystatus == "Completed":
            recipedir = workflow_setup.get_finished_recipe_dir(my_test_dir)
            myfile=MASTFile(os.path.join(recipedir,"SUMMARY.txt"))
            okays=0
            for myline in myfile.data:
                if "defect_int1_stat" in myline:
                    myenergy=myline.split()[-1]
                    self.assertEquals(myenergy, "-13.953")
                    okays=okays+1
                if "defect_vac1_stat" in myline:
                    myenergy=myline.split()[-1]
                    self.assertEquals(myenergy, "-10.623")
                    okays=okays+1
                if "defect_sub1_stat" in myline:
                    myenergy=myline.split()[-1]
                    self.assertEquals(myenergy, "-20.145")
                    okays=okays+1
            self.assertEquals(okays, 3)
            return
        else:
            self.assertTrue(False)
            return
    def test_neb_with_phonons(self):
        [mystatus, my_test_dir]=workflow_setup.generic_submit("neb_with_phonons.inp")
        if mystatus == "Unfinished":
            self.assertTrue(False)
            return
        elif mystatus == "Completed":
            recipedir = workflow_setup.get_finished_recipe_dir(my_test_dir)
            myfile=MASTFile(os.path.join(recipedir,"diffcoeff_utility","Diffusivity.txt"))
            okays=0
            for myline in myfile.data:
                if myline[0:3] == "1.0":
                    myval = myline.split()[1].strip()
                    compareval = "%3.1E" % 8.09778E-7
                    self.assertEquals("%3.1E" % float(myval),compareval)
                    #Not sure about 6.198170E-07.
                    #use results from 2014 pre structure-indexing
                    #and 2017 post-structure-indexing fixes, in savedoff
                    okays=okays+1
            self.assertEquals(okays, 1)
            return
        else:
            self.assertTrue(False)
            return
