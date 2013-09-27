import os
import time
import unittest
from unittest import SkipTest
import filecmp
from filecmp import dircmp
from MAST.mast import MAST
import MAST.utility.dirutil

testname ="mast_test"
oldcontrol = os.getenv("MAST_CONTROL")
oldrecipe = os.getenv("MAST_RECIPE_PATH")
oldscratch = os.getenv("MAST_SCRATCH")
print "Old directories:"
print oldcontrol
print oldrecipe
print oldscratch
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)


class TestMAST(unittest.TestCase):
    """Test MAST (mast -i input.inp)
    """
    def setUp(self):
        os.environ['MAST_CONTROL'] = testdir
        os.environ['MAST_RECIPE_PATH'] = testdir
        os.environ['MAST_SCRATCH'] = testdir
        os.chdir(testdir)

    def tearDown(self):
        os.environ['MAST_CONTROL'] = oldcontrol
        os.environ['MAST_RECIPE_PATH'] = oldrecipe
        os.environ['MAST_SCRATCH'] = oldscratch

    def test_check_indep_loop_no_loops(self):
        raise SkipTest
        print "New paths:"
        print os.getenv("MAST_CONTROL")
        print os.getenv("MAST_RECIPE_PATH")
        print os.getenv("MAST_SCRATCH")
        mymast = MAST.mast.MAST(inputfile="basic_test.inp",outputfile="output.inp")
        mymast.check_independent_loops()
    def test_set_up_recipe(self):
        mymast = MAST.mast.MAST(inputfile="basic_test.inp",
                                outputfile="output.inp")
        mymast.set_input_options()
        print "INPUT OPTIONS: ", mymast.input_options
        mymast.set_class_attributes()
        mymast.make_working_directory()
        print "WORKING DIRECTORY:", mymast.working_directory
        mymast.create_recipe_metadata()
        mymast.copy_input_file()
        mymast.parse_recipe_template()
        mymast.create_recipe_plan()
        mymast.create_archive_files()

