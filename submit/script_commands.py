#!/usr/bin/env python
# Purpose: Modify and get information out of submission scripts.
# Authors: Tam Mayeshiba
# 11/9/12 created
import os
import sys
import math #TA to calc num_nodes
from MAST.utility.mastfile import MASTFile
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
    script=MASTFile(templatepath)
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

def write_submit_script(keywords):
    """This script example is built on the following ingredient keywords,
        and may require significant customization.
        mast_processors
        mast_ppn
        mast_queue
        mast_nodes
        mast_exec
        mast_walltime
        mast_memory
    """
    import os
    #set defaults
    name = os.path.basename(keywords['name'])
    try:
        mast_processors = str(keywords['program_keys']['mast_processors'])
    except KeyError:
        mast_processors = "8"
    try:
        mast_ppn = str(keywords['program_keys']['mast_ppn'])
    except KeyError:
        mast_ppn = "8"
    try:
        mast_queue = str(keywords['program_keys']['mast_queue'])
    except KeyError:
        mast_queue = 'default'
    try:
        mast_nodes = str(keywords['program_keys']['mast_nodes'])
    except KeyError:
        mast_nodes = "1"
    try:
        mast_exec = str(keywords['program_keys']['mast_exec'])
    except KeyError:
        mast_exec = "mpiexec vasp"
    try:
        mast_walltime = str(keywords['program_keys']['mast_walltime'])
    except KeyError:
        mast_walltime = "24"
    try:
        mast_memory = str(keywords['program_keys']['mast_memory'])
    except KeyError:
        mast_memory = "1000"
    #create submission script
    myscript = MASTFile()
    myscript.data.append("#!/bin/bash" + "\n")
    myscript.data.append("#PBS -N " + name + "\n")
    myscript.data.append("#PBS -l nodes=" + mast_nodes + ":ppn=" + mast_ppn + ",pvmem=" + mast_memory + "mb" + "\n")
    myscript.data.append("#PBS -l walltime=" + mast_walltime + "\n")
    myscript.data.append("NN=`cat $PBS_NODEFILE | wc -l`" + "\n")
    myscript.data.append('echo "Processors received = "$NN' + "\n")
    myscript.data.append('echo "script running on host `hostname`"' + "\n")
    myscript.data.append('cd $PBS_O_WORKDIR' + "\n")
    myscript.data.append('echo "PBS_NODEFILE"' + "\n")
    myscript.data.append('cat $PBS_NODEFILE' + "\n")
    myscript.data.append(mast_exec + "\n")
    myscript.to_file(os.path.join(keywords['name'] + '/submit.sh'))
 

def write_submit_script_from_template(keywords):
    if not ('script' in keywords['program_keys'].keys()):
        templatename = 'submitscript.sh'
    else:
        templatename = keywords['program_keys']['script']
    templatepath = os.path.join(dirutil.get_mast_install_path(),
                    'submit',templatename)
    bname = os.path.basename(keywords['name'])
    wpath = keywords['name'] + '/submit.sh'
    #print wpath
    #print bname
    modify_jobname(templatepath, wpath, bname)
    return

