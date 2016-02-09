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
        os.chdir(testdir)
        if not os.path.exists("childdir"):
            os.mkdir("childdir")
        shutil.copy("files/metadata.txt","childdir")

    def tearDown(self):
        for fname in ["POSCAR","XDATCAR","DYNMAT","OSZICAR","DYNMAT_combined","KPOINTS","POTCAR","INCAR","WAVECAR","CHGCAR","POSCAR_no_sd","XDATCAR_combined","CONTCAR","metadata.txt"]:
            try:
                os.remove("childdir/%s" % fname)
            except OSError:
                pass
        os.rmdir("childdir")

    def test_none(self):
        self.assertTrue(True)
