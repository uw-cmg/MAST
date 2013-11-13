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

        logging.basicConfig(filename="%s/mast.log" % os.getenv("MAST_CONTROL"), level=logging.DEBUG)
        self.logger = logging.getLogger(__name__)

    def loop_through_errors(self):
        """Loop through all errors in the error handlers.
        """
        errct = 0
        return errct
    
