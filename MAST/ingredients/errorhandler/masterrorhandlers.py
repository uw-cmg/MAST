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
    def __init__(self, output_filename="", timeout=3600, archivelist=list(), copyfromlist=list(), copytolist=list()):
        """
        Detects an error when the output file has not been updated
        in timeout seconds. Perturbs structure and restarts
        Copied from custodian.vasp.handlers.FrozenJobErrorHandler,
            except no perturbation of the structure.
            Args: (archivelist is different from the custodian version)
                output_filename <str>: output filename
                timeout <int>: timeout seconds
                archivelist <list of str>: list of file names to archive
                copyfromlist <list of str>: list of file names to copy from
                copytolist <list of str>: list of file names to copy to; one-to
                        -one correspondence with copyfromlist (e.g. copy from 
                        CONTCAR to POSCAR for VASP)
            Returns:
                Archives files to error.#.tar.gz
        """
        self.output_filename = output_filename
        self.timeout = timeout
        self.archivelist = list(archivelist)
        self.copyfromlist = list(copyfromlist)
        self.copytolist = list(copytolist)

    def check(self): 
        st = os.stat(self.output_filename) 
        if time.time() - st.st_mtime > self.timeout: 
            return True

    def correct(self): 
        actions=list()
        backup(self.archivelist)
        actions.append("Archived files %s" % self.archivelist)
        for file in self.archivelist:
            if (not file in self.copyfromlist) and (not file in self.copytolist):
                os.remove(file)
            else:
                cindex = self.copyfromlist.index(file)
                os.rename(file, self.copytolist[cindex])
                actions.append("Copied file %s to %s" % (self.copyfromlist[cindex], self.copytolist[cindex]))
                self.archivelist.remove(self.copytolist[cindex])
        return {"errors": ["MAST Frozen job"], "actions": actions}

    @property
    def is_monitor(self): return True

    @property
    def to_dict(self): return {"@module": self.__class__.__module__, "@class": self.__class__.__name__, "output_filename": self.output_filename, "timeout": self.timeout}
