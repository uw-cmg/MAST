#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
# Purpose: Modify and get information out of submission scripts.
# 11/9/12 created

import os
import sys
import math #TA to calc num_nodes
import importlib

from MAST.utility.mastfile import MASTFile
from MAST.utility import MASTError

mast_control = os.getenv("MAST_CONTROL")
platform_file = MASTFile("%s/set_platform" % mast_control)
mast_platform = platform_file.data[0].strip()

#mysc = "submit.platforms.script_commands_%s" % mast_platform
#my_script_commands = importlib.import_module(mysc)



"""
Script commands
"""

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
        mast_walltime = str(keywords['program_keys']['mast_walltime']) + ":00:00"
    except KeyError:
        mast_walltime = "24:00:00"
    try:
        mast_memory = str(keywords['program_keys']['mast_memory'])
    except KeyError:
        mast_memory = "1000"
    #create submission script
    newkey = dict()
    newkey['mast_memory'] = mast_memory
    newkey['mast_exec'] = mast_exec
    newkey['mast_walltime'] = mast_walltime
    newkey['mast_queue'] = mast_queue
    newkey['mast_ppn'] = mast_ppn
    newkey['mast_nodes'] = mast_nodes
    newkey['mast_processors'] = mast_processors
    newkey['mast_name'] = name
    
    my_template = MASTFile("%s/platforms/%s/submit_template.sh" % (mast_control, mast_platform))
    newdata = list()
    for myline in my_template.data:
        for mykey in newkey.keys():
            querykey = "?" + mykey + "?"
            if querykey in myline:
                myline = myline.replace(querykey, newkey[mykey])
        newdata.append(myline)
    filled_template = MASTFile()
    filled_template.data = list(newdata)
    filled_template.to_file("%s/submit.sh" % keywords['name']) #use full path here
    return
