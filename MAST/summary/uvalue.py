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
from MAST.ingredients.checker import VaspChecker
from MAST.ingredients.checker import VaspNEBChecker
def main(ingname=""):
    """Get the U value from an ingredient.
        Args:
            ingname <str>: Ingredient name, full path
        Returns:
            <string>: "LDAUU, LDAUJ (eV);x x x, x x x", ex:
                "LDAUU, LDAUJ (eV); 0 5 0, 0 1 0"
    """
    trymeta = "%s/metadata.txt" % ingname
    frontstr = "LDAUU, LDAUJ (eV);"
    if os.path.isfile(trymeta):
        mymeta = Metadata(metafile=trymeta)
        myprogram = mymeta.read_data("program")
    else:
        myprogram = "None"
    if myprogram in ['vasp','vasp_neb']:
        if os.path.isdir("%s/01" % ingname):
            tryfile = "%s/01/OUTCAR" % ingname
        else:
            tryfile = "%s/OUTCAR" % ingname
        grepldauu = fileutil.grepme(tryfile, "LDAUU")
        if grepldauu == []:
            return "%s N/A" % frontstr
        ldauuparse = grepldauu[0].split("=")[1].strip()
        grepldauj = fileutil.grepme(tryfile, "LDAUJ")
        ldaujparse = grepldauj[0].split("=")[1].strip()
        return "%s %s, %s" % (frontstr, ldauuparse, ldaujparse)
    else:
        return "%s N/A" % frontstr
