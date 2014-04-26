##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import sys
import fnmatch
import time
import subprocess
from MAST.utility import MASTError
from MAST.utility.metadata import Metadata
"""General file utilities that do not involve opening a MASTFile
"""

def grepme(filename="", grepstr="", lastlines=""):
    """Use the 'grep' function to grep through a file.
        Args:
            filename <str>: Filename (full path preferred)
            grepstr <str>: String to grep.
            lastlines <str or int>: Search only last X lines.
                Leave blank to grep through entire file.
        Returns:
            Grep results as a list of strings
    """
    if lastlines == "":
        grepcmd = 'grep "%s" %s' % (grepstr, filename)
    else:
        lastlines = str(lastlines)
        grepcmd = 'tail -n %s %s | grep "%s"' % (lastlines, filename, grepstr)
    #print "GREP CMD: %s" % grepcmd
    grepproc = subprocess.Popen(grepcmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    grepcomm = grepproc.communicate()[0]
    grepresults = grepcomm.split("\n")
    grepresults.remove("") #Remove trailing carriage return
    return grepresults

