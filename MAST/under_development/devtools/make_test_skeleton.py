#!/usr/bin/env python
#TTM 20130128
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
"""Create a unit test framework using function names in a file

Writes to a new file.
Args:
    [$1] <str>: script name for which to create unit tests
"""
from MAST.utility import MASTFile
import os
import sys

if len(sys.argv) < 2:
    print "Call with a script name."
    sys.exit()
filename=""
filename=sys.argv[1]
basename=""
myfile = MASTFile(filename)
basename = os.path.basename(filename)[:-3]
newfile = MASTFile()
defline=""
baseline=""
newdefline=""
newfile.data.append("\"\"\"Tests for " + basename.title() + "\"\"\"" + '\n')
newfile.data.append("\n")
newfile.data.append("from MAST.." +basename+ " import "+ basename.title()+ "\n")
newfile.data.append("\n")
newfile.data.append("import unittest\n")
newfile.data.append("from unittest import SkipTest\n")
newfile.data.append("import os\n")
newfile.data.append("import time\n")
newfile.data.append("import MAST\n")
newfile.data.append("import pymatgen\n")
newfile.data.append("from MAST.utility import dirutil\n")
newfile.data.append("\n")
newfile.data.append("testname=\"<test_folder_here>\"\n")
newfile.data.append("testdir = os.path.join(os.getenv(\"MAST_INSTALL_PATH\"),\'test\',testname)\n")
newfile.data.append("\n")
newfile.data.append("class Test"+ basename.title() + "(unittest.TestCase):" + "\n")
newfile.data.append("\n")
newfile.data.append("    def setUp(self):" + '\n')
newfile.data.append("        os.chdir(testdir)\n")
newfile.data.append("\n")
newfile.data.append("    def tearDown(self):" + '\n')
newfile.data.append("        pass" + '\n')
newfile.data.append("\n")
for line in myfile.data:
    if "def " in line.strip()[0:4]:
        defline = line
        print defline
        if not "self" in defline:
            partone = myfile.get_segment(defline,"def ","(").strip()
            parttwo = myfile.get_segment(defline,"(",":")
            parttwo.strip()

            newfile.data.append("    def test_" + partone + "(self):" + '\n')
            newfile.data.append("        raise SkipTest\n")
            newfile.data.append("        #self.testclass." + partone + "(" + parttwo+'\n')
            newfile.data.append("\n")

        else:
            partone = myfile.get_segment(defline,"def ","(self").strip()
            parttwo = myfile.get_segment(defline,"self,",":")
            if parttwo == None:
                parttwo = myfile.get_segment(defline,"self",":")
            parttwo.strip()

            newfile.data.append("    def test_" + partone + "(self):" + '\n')
            newfile.data.append("        raise SkipTest\n")
            newfile.data.append("        #self.testclass." + partone + "(" + parttwo+'\n')
            newfile.data.append("\n")
newfile.to_unique_file(os.getcwd(), "test_" + basename, ".py")

