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
from custodian.custodian import ErrorHandler, backup
class MASTFrozenJobErrorHandler(ErrorHandler):
    """Check if a job is frozen. Prepare the job for resubmission.
    """
    def __init__(self, output_filename="", timeout=3600, archivelist=list()):
        """
        Detects an error when the output file has not been updated
        in timeout seconds. Perturbs structure and restarts
        Copied from custodian.vasp.handlers.FrozenJobErrorHandler,
            except no perturbation of the structure.
            Args: (archivelist is different from the custodian version)
                output_filename <str>: output filename
                timeout <int>: timeout seconds
                archivelist <list of str>: list of file names to archive
            Returns:
                Archives files to error.#.tar.gz
        """
        self.output_filename = output_filename
        self.timeout = timeout
        self.archivelist = list(archivelist)

    def check(self): 
        st = os.stat(self.output_filename) 
        if time.time() - st.st_mtime > self.timeout: 
            return True

    def correct(self): 
        backup(self.archivelist)
        actions=list()
        actions.append("Archived files %s" % self.archivelist)
        return {"errors": ["MAST Frozen job"], "actions": actions}

    @property
    def is_monitor(self): return True

    @property
    def to_dict(self): return {"@module": self.__class__.__module__, "@class": self.__class__.__name__, "output_filename": self.output_filename, "timeout": self.timeout}
