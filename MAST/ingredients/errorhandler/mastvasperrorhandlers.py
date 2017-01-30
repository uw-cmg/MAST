##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import time
import shutil
import pymatgen
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import fileutil
from MAST.utility import Metadata
from MAST.utility import MASTFile
from MAST.utility import loggerutils
from MAST.submit import queue_commands
from MAST.submit import script_commands
from pymatgen.io.vasp import Incar
from custodian.custodian import ErrorHandler
try:
   from custodian.custodian import backup
except:
   from custodian.utils import backup
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
        self.logger = loggerutils.get_mast_logger(self.ingpath)
    
    def check(self):
        ifile = "%s/INCAR" % self.ingpath
        if os.path.isfile(ifile):
            myincar = Incar.from_file(ifile)
        else:
            self.logger.warning("No INCAR file found at %s" % ifile)
            return False
        nsw = myincar.get('NSW',0)
        if nsw == 0:
            self.logger.info("NSW is 0 or not set.")
            return False
        ofile = "%s/OUTCAR" % self.ingpath
        if not (os.path.isfile(ofile)):
            ofile = "%s/01/OUTCAR" % self.ingpath
        if not (os.path.isfile(ofile)):
            self.logger.error("No OUTCAR file at %s or %s/01." % (self.ingpath, self.ingpath))
            return False
        myiters = fileutil.grepme(ofile, "Iter")
        if len(myiters) == 0:
            self.logger.warning("No Iterations in OUTCAR at %s" % ofile)
            return False
        lastentry = myiters[-1]
        ionichalf = lastentry.split('(')[0]
        lastiter = ionichalf.split()[-1]
        lastiter = int(lastiter)
        if lastiter >= nsw:
            self.logger.info("Last ionic step %s is at NSW %s" % (lastiter, nsw))
            return True
        return False
        #ofile = "%s/vasprun.xml" % self.ingpath
        #if os.path.isfile(ofile):
        #    myvasprun = pymatgen.io.vasp.Vasprun(ofile)
        #    if myvasprun.nionic_steps >= nsw:
        #        self.logger.info("Ionic steps %s is at NSW %s" % (myvasprun.nionic_steps, nsw))
        #        return True
        #else:
        #    self.logger.error("No vasprun.xml file found at %s/vasprun.xml" % self.ingpath)
        #    return False

    def correct(self): 
        actions=list()
        backup(self.archivelist)
        actions.append("Archived files %s" % self.archivelist)
        for myfile in self.archivelist:
            if (not myfile in self.copyfromlist) and (not myfile in self.copytolist): #if it is a copy-from, don't want to delete it, and if it is a copy-to, it will be either overwritted by a copy-from file, or preserved
                if os.path.isfile(myfile):
                    os.remove(myfile)
        for myfile in self.copyfromlist:
            cindex = self.copyfromlist.index(myfile)
            if os.path.isfile(myfile):
                if os.stat(myfile).st_size > 0:
                    os.rename(myfile, self.copytolist[cindex])
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

