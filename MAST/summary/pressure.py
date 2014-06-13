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
    """Get the last pressure value from an ingredient.
        Args:
            ingname <str>: Ingredient name, full path
        Returns:
            <string>: "Last pressure (kbar);0.00", ex:
                "Last pressure (kbar);-23.55"
    """
    trymeta = "%s/metadata.txt" % ingname
    frontstr = "Last pressure (kbar);"
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
        greppress = fileutil.grepme(tryfile, "external pressure")
        #"external pressure = -1.97 kB Pullay stress = 0.00 kB"
        if greppress == []:
            return "%s N/A" % frontstr
        mypress = greppress[-1].strip().split()[3]
        return "%s %s" % (frontstr, mypress)
    else:
        return "%s N/A" % frontstr
