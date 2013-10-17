import os

from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import Metadata
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifParser

class BaseError(MASTObj):
    """Base error handling class. This class switches between
        program-specific error-handling functions.
        For example, phon-related functions or vasp-related
        functions get their own error class.
    """
    
    def __init__(self, allowed_keys, **kwargs):
        allowed_keys_base = dict()
        allowed_keys_base.update(allowed_keys) 
        MASTObj.__init__(self, allowed_keys_base, **kwargs)
   
    def loop_through_errors(self):
        raise NotImplementedError
