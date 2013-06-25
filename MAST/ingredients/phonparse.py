import numpy as np
import pymatgen
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar
from pymatgen.io.vaspio import Kpoints
from pymatgen.io.vaspio import Potcar
from pymatgen.io.vaspio import Incar

from MAST.utility.mastobj import MASTObj
from MAST.ingredients.baseingredient import BaseIngredient
from MAST.ingredients.pmgextend import vasp_extensions
from MAST.ingredients.optimize import Optimize

import os
import shutil
#import pdb
#TTM

class PhonParse(Optimize):
    """Phonon parser using PHON
    """
    
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)
        self.keywords['program'] = 'phon'

