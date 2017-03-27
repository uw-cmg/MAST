import os
import time
import unittest
from unittest import SkipTest
import filecmp
from filecmp import dircmp
import MAST
import shutil
from MAST.utility import dirutil
from MAST.controllers.mastinput import MASTInput as controllerMAST
testname ="mast_test"
testdir = dirutil.get_test_dir(testname)
oldcontrol = os.getenv("MAST_CONTROL")
oldscratch = os.getenv("MAST_SCRATCH")
print "Old directories:"
print oldcontrol
print oldscratch


class TestMAST(unittest.TestCase):
    """Test MAST (mast -i input.inp)
    """
    def setUp(self):
        os.environ['MAST_CONTROL'] = testdir
        os.environ['MAST_SCRATCH'] = testdir
        os.chdir(testdir)

    def tearDown(self):
        basicdirs = MAST.utility.dirutil.walkdirs(testdir, 1,1,"*Basic*")
        for basicdir in basicdirs:
            shutil.rmtree(basicdir)
        os.environ['MAST_CONTROL'] = oldcontrol
        os.environ['MAST_SCRATCH'] = oldscratch

    def test_check_indep_loop_no_loops(self):
        raise SkipTest
        print "New paths:"
        print os.getenv("MAST_CONTROL")
        print os.getenv("MAST_SCRATCH")
        mymast = MAST.mast.MAST(inputfile="basic_test.inp",outputfile="output.inp")
        mymast.check_independent_loops()
    def test_set_up_recipe(self):
        mymast = controllerMAST(inputfile="basic_test.inp",
                                outputfile="output.inp")
        mymast.set_input_options()
        #print "INPUT OPTIONS: ", mymast.input_options
        mymast.set_class_attributes()
        mymast.make_working_directory()
        #print "WORKING DIRECTORY:", mymast.working_directory
        mymast.create_recipe_metadata()
        mymast.copy_input_file()
        mymast.parse_recipe_template()
        mymast.create_recipe_plan()
        mymast.create_archive_files()
        myfolders = MAST.utility.dirutil.walkdirs(testdir,1,1,"*Basic*")
        myfiles = MAST.utility.dirutil.walkfiles(myfolders[-1],1,1)
        subfolders = MAST.utility.dirutil.walkdirs(myfolders[-1],1,1)
        lacksfolders=0
        for folder in ['opt1','opt2']:
            if not os.path.join(myfolders[-1],folder) in subfolders:
                lacksfolders = lacksfolders + 1
        lacksfiles=0
        for file in ['archive_input_options.txt','metadata.txt',
                        'status.txt',
                        'archive_recipe_plan.txt', 'input.inp']:
            if not os.path.join(myfolders[-1],file) in myfiles:
                lacksfiles = lacksfiles + 1
        self.assertEqual(lacksfiles,0)
        self.assertEqual(lacksfolders,0)

