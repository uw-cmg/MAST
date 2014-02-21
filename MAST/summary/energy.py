import os
import sys
from MAST.utility import Metadata
from MAST.ingredients.checker import VaspChecker
from MAST.ingredients.checker import VaspNEBChecker
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
    if myprogram == 'vasp':
        if os.path.isdir("%s/01" % ingname):
            myprogram = 'vasp_neb' #probably a NEB static
        else:
            mychecker = VaspChecker(name=ingname)
            return "energy (eV);%3.3f" % mychecker.get_energy_from_energy_file()
    if myprogram == 'vasp_neb':
        mychecker = VaspNEBChecker(name=ingname)
        return "energy (eV);%s" % mychecker.get_energy_from_energy_file()
    else:
        return "energy (eV);N/A"
