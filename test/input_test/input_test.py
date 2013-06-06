import pymatgen
import os
import shutil
import sys
import unittest
import time
from filecmp import dircmp
from nose import SkipTest

from MAST.parsers import InputParser
from MAST.parsers import InputPythonCreator
from MAST.utility import MASTError
from MAST.utility.dirutil import *

class TestInput(unittest.TestCase):
    def setUp(self):
        scripts=get_mast_install_path()
        os.chdir(os.path.join(scripts,'test','input_test'))
    
    def tearDown(self):
        pass

    def test_inputparser(self):
        raise SkipTest
        myip = InputParser(inputfile='mast.inp')
        myoptions=myip.parse()
        print myoptions
        self.assertEqual(True,True)

    def test_inputpythoncreator(self):
        myip = InputParser(inputfile='mast.inp')
        myoptions=myip.parse()
        myipc = InputPythonCreator(optiondict=myoptions) 
        myipc.print_input_options()
        self.assertEqual(0,0)

