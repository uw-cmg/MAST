##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import sys
from MAST.utility import Metadata
from MAST.utility import dirutil
from MAST.ingredients.checker import VaspChecker
def main(ingname=""):
    """Get the energy from an ingredient.
        Args:
            ingname <str>: Ingredient name, full path
        Returns:
            <string>: "energy (eV);<energy as a string>", ex:
                "energy (eV); 3.0"
                Returns last E0 energy for a VASP run.
                Returns last E0 energy for all images for a VASP neb.
                Returns "N/A" otherwise.
    """
    trymeta = "%s/metadata.txt" % ingname
    if os.path.isfile(trymeta):
        mymeta = Metadata(metafile=trymeta)
        myprogram = mymeta.read_data("program")
    else:
        myprogram = "None"
    if 'induce' in ingname: #skip inducedefect ingredients
        myprogram = "None"
    if 'vasp' in myprogram:
        if os.path.isdir("%s/01" % ingname):
            estr = "energies (eV)"
            for subdir in dirutil.immediate_subdirs(ingname):
                mychecker = VaspChecker(name="%s/%s" % (ingname, subdir))
                estr = estr + ";%3.3f" % mychecker.get_energy_from_energy_file()
            return estr
        else:
            mychecker = VaspChecker(name=ingname)
            return "energy (eV);%3.3f" % mychecker.get_energy_from_energy_file()
    else:
        return "energy (eV);N/A"
