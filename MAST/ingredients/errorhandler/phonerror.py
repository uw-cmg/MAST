<<<<<<< HEAD
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# PHON error handler is obsolete, as PHON is no longer being supported.
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
=======
>>>>>>> 1e7ff1933f44df0bd90ea0325109658e94e0222c
from custodian.vasp import handlers
from pymatgen.core.structure import Structure
import inspect
import os
import logging
from MAST.ingredients.errorhandler import BaseError
class PhonError(BaseError):
    """PHON error-handling functions 
    """
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseError.__init__(self, allowed_keys, **kwargs)

    
