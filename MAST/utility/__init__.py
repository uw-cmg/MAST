# utility/__init__.py

import numpy as np
import pymatgen as pmg

from MAST.utility.mastobj import MASTObj
from MAST.utility.masterror import MASTError
from MAST.utility.inputoptions import InputOptions
from MAST.utility.mastfile import MASTFile
from MAST.utility.metadata import Metadata

def MAST2Structure(lattice=None, coordinates=None, atom_list=None, coord_type='fractional'):
    """Helper function for converting input options into a pymatgen Structure object"""
    if (coord_type.lower() == 'fractional'):
        return pmg.Structure(lattice, atom_list, coordinates)
    elif (coord_type.lower() == 'cartesian'):
        return pmg.Structure(lattice, atom_list, coordinates,
                             coords_are_cartesian=True)
