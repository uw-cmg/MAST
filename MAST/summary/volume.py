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
    """Get the last volume value from an ingredient.
        Args:
            ingname <str>: Ingredient name, full path
        Returns:
            <string>: "Last volume (Angstroms^3);0.000", ex:
                "Last volume (Angstroms^3); 324.456"
    """
    trymeta = "%s/metadata.txt" % ingname
    frontstr = "Last volume (Angstroms^3);"
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
        grepvol = fileutil.grepme(tryfile, "volume of cell")
        if grepvol == []:
            return "%s N/A" % frontstr
        myvol = grepvol[-1].split(":")[1].strip()
        return "%s %s" % (frontstr, myvol)
    else:
        return "%s N/A" % frontstr
