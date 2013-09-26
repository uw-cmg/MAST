import os

from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import Metadata
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifParser

class BaseChecker(MASTObj):
    """Base checker class. This class switches between
        program-specific functions.
        For example, phon-related functions or vasp-related
        functions get their own checker class.
    """
    
    def __init__(self, allowed_keys, **kwargs):
        allowed_keys_base = dict()
        allowed_keys_base.update(allowed_keys) 
        MASTObj.__init__(self, allowed_keys_base, **kwargs)
    
    def is_complete(self):
        raise NotImplementedError
    
    def get_coordinates_only_structure_from_input(self):
        """Get coordinates-only structures from mast_coordinates
            ingredient keyword
            Args:
                keywords <dict>: ingredient keywords
            Returns:
                coordstrucs <list>: list of Structure objects
        """
        coordposlist=self.keywords['program_keys']['mast_coordinates']
        coordstrucs=list()
        coordstruc=None
        for coordpositem in coordposlist:
            if ('poscar' in os.path.basename(coordpositem).lower()):
                coordstruc = Poscar.from_file(coordpositem).structure
            elif ('cif' in os.path.basename(coordpositem).lower()):
                coordstruc = CifParser(coordpositem).get_structures()[0]
            else:
                error = 'Cannot build structure from file %s' % coordpositem
                raise MASTError(self.__class__.__name__, error)
            coordstrucs.append(coordstruc)
        return coordstrucs
