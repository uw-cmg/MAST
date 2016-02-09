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
#oldcontrol = os.getenv("MAST_CONTROL")
#oldrecipe = os.getenv("MAST_RECIPE_PATH")
#oldscratch = os.getenv("MAST_SCRATCH")
#print "Old directories:"
#print oldcontrol
#print oldrecipe
#print oldscratch
testdir = dirutil.get_test_dir(testname)


class TestWorkflows(unittest.TestCase):
    """Test Workflows
    """
    def setUp(self):
        pass
        return
    def tearDown(self):
        testlist=list()
        testlist.append("simple_optimization.inp")
        for testname in testlist:
            shortname = testname.split(".")[0]
            try:
                pass
                #os.remove("output_%s" % shortname)
            except OSError:
                pass
        myfiles = os.listdir(testdir)
        for myfile in myfiles:
            if "output_workflow_testing" in myfile:
                os.remove(myfile)
        return
    def test_none(self):
        self.assertTrue(True)
