import os
import time
import shutil
import pymatgen
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import Metadata
from MAST.utility import MASTFile
from MAST.utility import loggerutils
from submit import queue_commands
from submit import script_commands
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifParser
from custodian.custodian import ErrorHandler, backup
import logging

class VaspReachedNSWErrorHandler(ErrorHandler):
    """Check if a VASP job has reached NSW.
        Typically it should converge long before NSW.
        If NSW = 0 , ignore this error.
    """
    def __init__(self, ingpath, archivelist=list(), copyfromlist=list(), copytolist=list()):
        """
            Args: (archivelist is different from the custodian version)
                ingpath <str>: ingredient path
                archivelist <list of str>: list of file names to archive
                copyfromlist <list of str>: list of file names to copy from
                copytolist <list of str>: list of file names to copy to; one-to
                        -one correspondence with copyfromlist (e.g. copy from 
                        CONTCAR to POSCAR for VASP)
            Returns:
                Archives files to error.#.tar.gz
        """
        self.ingpath = ingpath
        self.archivelist = list(archivelist)
        self.copyfromlist = list(copyfromlist)
        self.copytolist = list(copytolist)
        self.logger = logging.getLogger(self.ingpath)
        self.logger = loggerutils.add_handler_for_recipe(self.ingpath, self.logger)
    
    def check(self):
        ifile = "%s/INCAR" % self.ingpath
        if os.path.isfile(ifile):
            myincar = pymatgen.io.vaspio.Incar.from_file(ifile)
        else:
            self.logger.warning("No INCAR file found at %s" % ifile)
            return False
        nsw = myincar.get('NSW',0)
        if nsw == 0:
            self.logger.info("NSW is 0 or not set.")
            return False
        ofile = "%s/vasprun.xml" % self.ingpath
        if os.path.isfile(ofile):
            myvasprun = pymatgen.io.vaspio.Vasprun(ofile)
            if myvasprun.nionic_steps >= nsw:
                self.logger.info("Ionic steps %s is at NSW %s" % (myvasprun.nionic_steps, nsw))
                return True
        else:
            self.logger.error("No vasprun.xml file found at %s/vasprun.xml" % self.ingpath)
            return False

    def correct(self): 
        actions=list()
        backup(self.archivelist)
        actions.append("Archived files %s" % self.archivelist)
        for file in self.archivelist:
            if (not file in self.copyfromlist) and (not file in self.copytolist): #if it is a copy-from, don't want to delete it, and if it is a copy-to, it will be either overwritted by a copy-from file, or preserved
                if os.path.isfile(file):
                    os.remove(file)
        for file in self.copyfromlist:
            cindex = self.copyfromlist.index(file)
            if os.path.isfile(file):
                if os.stat(file).st_size > 0:
                    os.rename(file, self.copytolist[cindex])
                    actions.append("Copied file %s to %s" % (self.copyfromlist[cindex], self.copytolist[cindex]))
                else:
                    actions.append("Skipped file copy of %s to %s because %s was empty." % (self.copyfromlist[cindex],self.copytolist[cindex],self.copyfromlist[cindex]))
            else:
                actions.append("Skipped file copy of %s to %s because %s did not exist." % (self.copyfromlist[cindex],self.copytolist[cindex],self.copyfromlist[cindex]))
        return {"errors": ["MAST VASP at NSW error"], "actions": actions}

    @property
    def is_monitor(self): return True

    @property
    def to_dict(self): return {"@module": self.__class__.__module__, "@class": self.__class__.__name__, "output_filename": self.output_filename, "timeout": self.timeout}

