#!/usr/bin/env python
# Purpose: Modify and get information out of submission scripts.
# Authors: Tam Mayeshiba
# 11/9/12 created
import os
import sys
import math #TA to calc num_nodes
from MAST.utility import fileutil
from MAST.utility import MASTError

"""
Script commands
"""

def get_script(templatepath):
    if not os.path.isfile(templatepath):
        raise MASTError("submit script_commands",
            "Could not find file at " + templatepath)
    if not templatepath[-3:] == '.sh':
        raise MASTError("submit script_commands", 
        "Script name does not end with .sh to signify a submission script.")
    script=fileutil.MASTFile(templatepath)
    if script.data == []:
        raise MASTError("submit script_commands",
        "Script empty at " + templatepath)
    return script

def modify_jobname(templatepath, writepath, jobname):
    """
        Change the jobname in the submission script.
        INPUTS:
            writepath <str> = full path to submission script file
            jobname <str> = new job name to use
        OUTPUTS:
            Modifies the file at writepath to use the new jobname.
    """
    myscript = get_script(templatepath)
    
    #begin platform specific:
    namenum=myscript.get_line_match_number("PBS -N")
    myscript.modify_file_by_line_number(namenum, "R", "#PBS -N " + jobname + "\n")
    #end platform specific
    
    myscript.to_file(writepath)
    return True
