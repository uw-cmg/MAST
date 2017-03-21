"""Tests for Diffusion Coefficient Tool"""

from MAST.utility import MASTError
from MAST.controllers.mastmon import MASTMon
import unittest
import numpy as np
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
import subprocess
from MAST.utility import dirutil
from MAST.utility import MASTFile
from MAST.utility import Metadata
from MAST.utility.DiffusionCoefficient import DiffCoeff
import shutil
testname="diffusion_coefficient_test"
testdir = dirutil.get_test_dir(testname)
workdir = os.path.join(testdir,'workdir')
if not os.path.isdir(workdir):
    os.mkdir(workdir)
old_control = os.getenv("MAST_CONTROL")
old_archive = os.getenv("MAST_ARCHIVE")
old_scratch = os.getenv("MAST_SCRATCH")

class TestDiffusionCoefficientTool(unittest.TestCase):
    def __init__(self, methodName='runTest'):
        unittest.TestCase.__init__(self, methodName)
        self.tname = self.id().split(".")[-1]
        self.wdir = os.path.join(workdir, self.tname)
        return

    def setUp(self):
        os.chdir(testdir)
        if not os.path.isdir(self.wdir):
            os.mkdir(self.wdir)
        return

    def tearDown(self):
        #pass
        #return
        for removedir in [self.wdir]:
            if os.path.isdir(removedir):
                shutil.rmtree(removedir, ignore_errors=True)
        os.environ["MAST_CONTROL"] = old_control
        os.environ["MAST_SCRATCH"] = old_scratch
        os.environ["MAST_ARCHIVE"] = old_archive

    def test__init__(self):
        raise SkipTest
        return
    
    def test_tool_without_wrapper(self):
        """Test diffusion coefficient tool without the wrapper
        """
        #set up paths
        rwdir = os.path.join(self.wdir, 'diffcoeff_utility')
        shutil.copytree(os.path.join(testdir,'files','diffcoeff_utility'), 
                        rwdir)
        outputname = os.path.join(rwdir, 'Diffusivity.txt')
        plotname = os.path.join(rwdir, 'Diffusivity.png')
        os.remove(outputname) #remove test output
        os.remove(plotname) #remove test output
        print(os.listdir(rwdir))
        #run the utility
        dctool_input = os.path.join(rwdir, 'diffcoeff_input.txt')
        curdir=os.getcwd()
        os.chdir(rwdir)
        DiffCoeff(dctool_input).calculatingD()
        os.chdir(curdir)
        #check output.
        print "output in: %s" % (outputname)
        with open(outputname, 'r') as opfile:
            olines = opfile.readlines()
        comparename = os.path.join(testdir,'files','diffcoeff_utility',
                                    'Diffusivity.txt')
        with open(comparename, 'r') as cfile:
            clines = cfile.readlines()
        print clines
        print olines
        self.assertEqual(len(clines),len(olines))
        for cidx in range(0, len(clines))[-5:]: #compare last few numeric lines
            [ctemp, cval] = clines[cidx].strip().split()
            [otemp, oval] = olines[cidx].strip().split()
            self.assertEqual(ctemp, otemp)
            self.assertEqual(cval, oval)
        return
    
    def test_tool_with_wrapper(self):
        """Test diffusion coefficient tool with the wrapper
        """
        #set up paths
        rwdir = os.path.join(self.wdir, 'diffcoeff_utility')
        shutil.copytree(os.path.join(testdir,'files','diffcoeff_utility'), 
                        rwdir)
        outputname = os.path.join(rwdir, 'Diffusivity.txt')
        plotname = os.path.join(rwdir, 'Diffusivity.png')
        os.remove(outputname) #remove test output
        os.remove(plotname) #remove test output
        print(os.listdir(rwdir))
        #run the utility
        curdir=os.getcwd()
        os.chdir(rwdir)
        mycommand="mast_diffusion_coefficient -i diffcoeff_input.txt"
        dproc = subprocess.Popen(mycommand, shell=True,
                        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        dproc.wait()
        print(dproc.communicate())
        os.chdir(curdir)
        #check output
        print "output in: %s" % (outputname)
        with open(outputname, 'r') as opfile:
            olines = opfile.readlines()
        comparename = os.path.join(testdir,'files','diffcoeff_utility',
                                    'Diffusivity.txt')
        with open(comparename, 'r') as cfile:
            clines = cfile.readlines()
        print clines
        print olines
        self.assertEqual(len(clines),len(olines))
        for cidx in range(0, len(clines))[-5:]: #compare last few numeric lines
            [ctemp, cval] = clines[cidx].strip().split()
            [otemp, oval] = olines[cidx].strip().split()
            self.assertEqual(ctemp, otemp)
            self.assertEqual(cval, oval)
        return
    
    def test_bcc_9freq(self):
        """Test BCC 9 frequency
        """
        #set up paths
        rwdir = os.path.join(self.wdir, 'diffcoeff_utility')
        shutil.copytree(os.path.join(testdir,'files','diffcoeff_utility_bcc_9freq_WAg'), 
                        rwdir)
        outputname = os.path.join(rwdir, 'Diffusivity.txt')
        plotname = os.path.join(rwdir, 'Diffusivity.png')
        os.remove(outputname) #remove test output
        #os.remove(plotname) #remove test output
        print(os.listdir(rwdir))
        #run the utility
        dctool_input = os.path.join(rwdir, 'diffcoeff_input.txt')
        curdir=os.getcwd()
        os.chdir(rwdir)
        DiffCoeff(dctool_input).calculatingD()
        os.chdir(curdir)
        #check output.
        print "output in: %s" % (outputname)
        with open(outputname, 'r') as opfile:
            olines = opfile.readlines()
        comparename = os.path.join(testdir,'files',
                            'diffcoeff_utility_bcc_9freq_WAg',
                            'Diffusivity.txt')
        with open(comparename, 'r') as cfile:
            clines = cfile.readlines()
        print clines
        print olines
        self.assertEqual(len(clines),len(olines))
        for cidx in range(0, len(clines))[-5:]: #compare last few numeric lines
            [ctemp, cval] = clines[cidx].strip().split()
            [otemp, oval] = olines[cidx].strip().split()
            self.assertEqual(ctemp, otemp)
            self.assertEqual(cval, oval)
        return
