import os
import time
import shutil
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

    def check_output_frozen(self, filetocheck):
        """Check if the output file is frozen.
            Args:
                filetocheck <str>: name of output file
            Return:
                True if the ingredient's output file has not been updated in
                    ingredient keyword mast_frozenseconds seconds
                    (default 7200)
                False otherwise
                False if no output file exists
        """
        raise MASTError("Use custodian instead")
        fullfile=os.path.join(self.keywords['name'],filetocheck)
        if not os.path.isfile(fullfile):
            self.logger.warning("No file at %s for frozen check." % fullfile)
            return False
        fstat=os.stat(fullfile)
        last_modified=fstat.st_mtime
        current_time=time.now()
        lag_time = current_time - last_modified
        frozen_time=0
        if not 'mast_frozenseconds' in self.program_keys.keys():
            frozen_time=7200
        else:
            frozen_time=float(self.program_keys['mast_frozenseconds'])
        if lag_time > frozen_time:
            return True
        else:
            return False
