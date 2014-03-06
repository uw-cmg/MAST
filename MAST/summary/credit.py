import os
import sys
import time
from MAST.utility import fileutil
def credit(recipedir):
    """Use mast_exec lines from the recipe directory in order
        to assign a CITATIONS section to the SUMMARY.txt file.
        Args:
            recipedir <str>: Recipe directory
        Returns:
            linelist <list of str>: Lines to add to SUMMARY.txt
    """
    myinput = os.path.join(recipedir, "input.inp")
    mastexecs = fileutil.grepme(myinput, "mast_exec")
    linelist=list()
    add_vasp=False
    add_phon=False
    add_pymatgen=True #MAST uses so much pymatgen it should always be cited
    add_mast=True
    for mastexec in mastexecs:
        if 'vasp' in mastexec.lower():
            add_vasp = True
        if 'phon' in mastexec.lower():
            add_phon = True
    linelist.append("==================================")
    linelist.append("----------- CITATIONS ------------")
    if add_mast:
        linelist.append("Cite MAST here.")
    if add_pymatgen:
        linelist.append("Cite pymatgen here.")
    if add_phon:
        linelist.append("Cite PHON here.")
    if add_vasp:
        linelist.append("Cite VASP here.")

