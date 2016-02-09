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
            except OSError:
                pass
        myfiles = os.listdir(testdir)
        for myfile in myfiles:
            if "output_workflow_testing" in myfile:
                os.remove(myfile)
        return

    def test_simple_optimization(self):
        mystatus=workflow_setup.generic_submit("simple_optimization.inp")
        if mystatus == "Unfinished":
            self.assertTrue(False)
            return
        elif mystatus == "Completed":
            #do more checks here
            self.assertTrue(True)
            return
        else:
            self.assertTrue(False)
            return
