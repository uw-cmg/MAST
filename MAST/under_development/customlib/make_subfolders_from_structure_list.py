##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba - file should be removed when structopt integrated; look for edits by Amy Kaczmarowski
# Last updated: 2014-04-25
##############################################################
import os
import sys
import logging
import pymatgen
from MAST.utility import loggerutils
from MAST.utility import MASTFile
from pymatgen.io.vasp import Poscar

def main(mydir="", fname="", startat=0, childdir=""):
    """This method makes POSCAR-containing subfolders from a 
        list of structures,
        in childdir if a child directory is given, or in 
        mydir otherwise.
        File fname must reside in directory mydir.
        Args:
            mydir <str>: Ingredient directory
            fname <str>: File containing list of structures
            startat <str, will be converted to int>: 
                0 (default) - start subfolders at 00
                1 - start subfolders at 01
            childdir <str>: Child directory (optional)
        If xyz files, fname is required to be of the form:
            "xyz"
            <float a> <float b> <float c for box size>
            fileone.xyz
            filetwo.xyz
            ...
        If POSCAR or other structure files, fname may simply 
        start with file names:
            POSCAR_01
            POSCAR_02
    """
    logger = logging.getLogger(mydir)
    logger = loggerutils.add_handler_for_recipe(mydir, logger)

    if childdir == "":
        childdir = mydir
    fpath = os.path.join(mydir, fname)
    if not os.path.isfile(fpath):
        logger.error("No file found at %s" % fpath)
        return "File not found. No effect."
    str_list_file = MASTFile(fpath)
    if str_list_file.data[0].strip().lower() == "xyz":
        namelist = list(str_list_file.data[2:])
        structure_list = get_xyz_structures(mydir, str_list_file.data[1], namelist)
    else:
        namelist = list(str_list_file.data)
        structure_list = get_structures(mydir, namelist)
    if not (len(structure_list) == len(namelist)):
        logger.error("Not all structures in file %s found." % fname)
        return "Not all structures found. No effect."
    strct = int(startat)
    for one_structure in structure_list:
        subname = str(strct).zfill(2)
        subpath = os.path.join(childdir, subname)
        if not os.path.isdir(subpath):
            os.mkdir(subpath)
        one_poscar = Poscar(one_structure)
        pospath = os.path.join(subpath, "POSCAR")
        if os.path.isfile(pospath):
            logger.error("POSCAR file already exists at %s" % subpath)
        else:
            one_poscar.write_file(os.path.join(subpath, "POSCAR"))
        strct = strct + 1
    return "Wrote all POSCAR files to subdirectories in %s" % childdir

def get_xyz_structures(mydir, boxline, namelist):
    """Get xyz structures from list
        Args:
            mydir <str>: Ingredient directory
            boxline <str>: Box size line, example:
                2.0 3.0 4.0
                for an orthorhombic box with lattice vector
                a = 2.0, b = 3.0, and c = 4.0
            namelist <list: List of file names of structures
        Returns:
            return_structures <list of Structure objects>
    """
    logger = logging.getLogger(mydir)
    logger = loggerutils.add_handler_for_recipe(mydir, logger)
    return_structures = list()
    boxline = boxline.strip()
    boxsplit = boxline.split()
    boxa = float(boxsplit[0])
    boxb = float(boxsplit[1])
    boxc = float(boxsplit[2])
    
    for strline in namelist:
        tryline = strline.strip()
        if os.path.isfile(tryline):
            str_path = tryline
        else:
            trypath = os.path.join(mydir, tryline)
            if os.path.isfile(trypath):
                str_path = trypath
            else:
                logger.error("No file found at %s" % tryline)
                continue
        one_mol = pymatgen.io.smartio.read_mol(str_path)
        one_structure = one_mol.get_boxed_structure(boxa, boxb, boxc)
        return_structures.append(one_structure)
    return return_structures

def get_structures(mydir, namelist):
    """Get structures from list
        Args:
            mydir <str>: Ingredient directory
            namelist <list: List of file names of structures
        Returns:
            return_structures <list of Structure objects>
    """
    logger = logging.getLogger(mydir)
    logger = loggerutils.add_handler_for_recipe(mydir, logger)
    return_structures = list()
    
    for strline in namelist:
        tryline = strline.strip()
        if os.path.isfile(tryline):
            str_path = tryline
        else:
            trypath = os.path.join(mydir, tryline)
            if os.path.isfile(trypath):
                str_path = trypath
            else:
                logger.error("No file found at %s" % tryline)
                continue
        one_structure = pymatgen.io.smartio.read_structure(str_path)
        return_structures.append(one_structure)
    return return_structures

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print "Not enough arguments given. Returning."
        sys.exit()
    myreturn = main(sys.argv[1],sys.argv[2],sys.argv[3])
    print myreturn
    sys.exit()
