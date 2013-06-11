import pymatgen
import os
import shutil
import sys
import unittest
import time
from filecmp import dircmp
from nose import SkipTest

from MAST.parsers.indeploopinputparser import IndepLoopInputParser
from MAST.utility import MASTError
from MAST.utility.dirutil import *

class TestIndepLoopInput(unittest.TestCase):
    def setUp(self):
        scripts=get_mast_install_path()
        os.chdir(os.path.join(scripts,'test','loopingtest'))
    
    def tearDown(self):
        pass

    def test_scan(self):
        myilp = IndepLoopInputParser(inputfile='mast2.inp')
        loopdict=myilp.scan_for_indep_loop()
        #print loopdict
        looplines=myilp.prepare_looped_lines(loopdict)
        print looplines
        loopdata=myilp.prepare_looped_datasets(looplines)
        #print loopdata
        myilp.create_input_files(loopdata)
    
    def test_noloop(self):
        myilp = IndepLoopInputParser(inputfile='mast.inp')
        loopdict=myilp.scan_for_indep_loop()
        #print loopdict
        looplines=myilp.prepare_looped_lines(loopdict)
        print looplines
        loopdata=myilp.prepare_looped_datasets(looplines)
        #print loopdata
        myilp.create_input_files(loopdata)
