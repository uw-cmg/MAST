import os
import time
import shutil
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import Metadata
from MAST.utility import MASTFile
from submit import queue_commands
from submit import script_commands
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifParser
from custodian.custodian import ErrorHandler, backup

class MASTWalltimeErrorHandler(ErrorHandler):
    """Check if a job has exceeded its walltime."""
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
    
    def check(self):
        errfilepath = queue_commands.get_job_error_file(self.ingpath)
        if errfilepath == None:
            return False
        errfile = MASTFile(errfilepath)
        for errline in errfile.data:
            if 'walltime' in errline.lower():
                return True
    
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
        return {"errors": ["MAST exceeded walltime error"], "actions": actions}

    @property
    def is_monitor(self): return True

    @property
    def to_dict(self): return {"@module": self.__class__.__module__, "@class": self.__class__.__name__, "output_filename": self.output_filename, "timeout": self.timeout}

class MASTMemoryErrorHandler(ErrorHandler):
    """Check if a job has insufficient virtual memory.
        If found, increases number of nodes and processors; typically
        the virtual memory requested is already at a max (at least for
        CMG queues).
    """
    def __init__(self, ingpath):
        """
            Args:
                ingpath <str>: ingredient path
            Returns:
                modifies submission script to add more nodes
        """
        self.ingpath = ingpath
    
    def check(self):
        errfilepath = queue_commands.get_job_error_file(self.ingpath)
        if errfilepath == None:
            return False
        errfile = MASTFile(errfilepath)
        for errline in errfile.data:
            if 'insufficient virtual memory' in errline.lower():
                return True
    
    def correct(self): 
        actions=list()
        multiplier = 4
        if ['mast_nodes'] in self.keywords['program_keys'].keys():
            currnodes = self.keywords['program_keys']['mast_nodes']
            newnodes = int(currnodes) * multiplier
            self.keywords['program_keys']['mast_nodes'] = newnodes
            actions.append("Multiplied mast_nodes by %i" % (multiplier, newnodes))
        if ['mast_processors'] in self.keywords['program_keys'].keys():
            currprocs = self.keywords['program_keys']['mast_processors']
            newprocs = int(currprocs) * multiplier
            self.keywords['program_keys']['mast_processors'] = newprocs
            actions.append("Multiplied mast_processors by %i" % (multiplier, newprocs))
        script_commands.write_submit_script(self.keywords)
        actions.append("Wrote new submission script.")
        return {"errors": ["MAST insufficient virtual memory error"], "actions": actions}

    @property
    def is_monitor(self): return True

    @property
    def to_dict(self): return {"@module": self.__class__.__module__, "@class": self.__class__.__name__, "output_filename": self.output_filename, "timeout": self.timeout}


