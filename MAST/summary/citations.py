##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import sys
import time
from MAST.utility import fileutil
from MAST.utility import MASTFile
from MAST.utility import dirutil
def get_citations(recipedir):
    """Use mast_exec lines from the recipe directory in order
        to assign a CITATIONS section to the SUMMARY.txt file.
        Args:
            recipedir <str>: Recipe directory
        Returns:
            linelist <list of str>: Lines to add to SUMMARY.txt
    """
    mylist=list()
    mylist=['mast','pymatgen','spglib'] #need to use; removed python as standard
    myinput = os.path.join(recipedir, "input.inp")
    mastexecs = fileutil.grepme(myinput, "mast_exec")
    for mastexec in mastexecs:
        if 'vasp' in mastexec.lower():
            mylist.append('vasp')
            mastxcs = fileutil.grepme(myinput, "mast_xc")
            if len(mastxcs) >= 1:
                mylist.append('vasp_pps')
                for mastxc in mastxcs:
                    if ("pw" in mastxc.lower()) or ("paw" in mastxc.lower()):
                        mylist.append('vasp_paw')
        #if 'phon' in mastexec.lower():
        #    mylist.append('phon')
        if 'structopt' in mastexec.lower():
            mylist.append('structopt')
    citationpath = os.path.join(dirutil.get_mast_install_path(),"MAST","summary","citations")
    citationfiles = os.listdir(citationpath)
    linelist = list()
    for listitem in mylist:
        for citationfile in citationfiles:
            if citationfile[:-3] == listitem: #all except number
                cfpath = os.path.join(citationpath, citationfile)
                cffile = MASTFile(cfpath)
                linelist.extend(cffile.data)
    return linelist
