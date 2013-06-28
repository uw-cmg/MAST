import pymatgen
import os
import shutil
import sys
import unittest
import time
import filecmp
from filecmp import dircmp
from nose import SkipTest

from MAST.parsers import InputParser
from MAST.parsers import InputPythonCreator
from MAST.utility import MASTError
from MAST.utility import MASTFile
from MAST.utility.dirutil import *

class TestInput(unittest.TestCase):
    def setUp(self):
        scripts=get_mast_install_path()
        os.chdir(os.path.join(scripts,'test','input_test'))
    
    def tearDown(self):
        if os.path.isfile("input_python_created"):
            os.remove("input_python_created")

    def test_inputparser(self):
        raise SkipTest
        myip = InputParser(inputfile='mast.inp')
        myoptions=myip.parse()
        print myoptions
        self.assertEqual(True,True)

    def test_inputpythoncreator(self):
        raise SkipTest
        myip = InputParser(inputfile='mast.inp')
        myoptions=myip.parse()
        myipc = InputPythonCreator(input_options=myoptions) 
        mylines=myipc.print_input_options()
        myfile = MASTFile()
        myfile.data = mylines
        myfile.to_file("./input_python_created")
        #print mylines
        self.assertTrue(filecmp.cmp("input_python_created","test_input_python_created"))

