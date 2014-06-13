##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import sys
from MAST.utility import Metadata
from MAST.utility import fileutil
from MAST.utility import MASTFile
from MAST.ingredients.checker import VaspChecker
from MAST.ingredients.checker import VaspNEBChecker
def main(ingname=""):
    """Get the last lattice from an ingredient.
        Args:
            ingname <str>: Ingredient name, full path
        Returns:
            <string>: "Last lattice (Angstroms);scale, a, b, c", ex:
                "Last lattice (Angstroms); scale 1, 4.01 0.2 0, 0 4.01 0, 
                        0.1 0.1 4.01"
    """
    trymeta = "%s/metadata.txt" % ingname
    frontstr = "Last lattice (Angstroms);"
    if os.path.isfile(trymeta):
        mymeta = Metadata(metafile=trymeta)
        myprogram = mymeta.read_data("program")
    else:
        myprogram = "None"
    if myprogram in ['vasp','vasp_neb']:
        if os.path.isdir("%s/01" % ingname):
            tryfile = "%s/01/CONTCAR" % ingname
        else:
            tryfile = "%s/CONTCAR" % ingname
        cfile = MASTFile(tryfile)
        if len(cfile.data) == 0:
            return "%s N/A" % frontstr
        scale = cfile.data[1].strip()
        avec = cfile.data[2].strip()
        bvec = cfile.data[3].strip()
        cvec = cfile.data[4].strip()
        return "%s scale %s, %s, %s, %s" % (frontstr, scale, avec, bvec, cvec)
    else:
        return "%s N/A" % frontstr
