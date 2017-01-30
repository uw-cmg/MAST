##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import time
import shutil
import logging
import inspect
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import Metadata
from MAST.utility import loggerutils
from MAST.ingredients.errorhandler import masterrorhandlers 
from MAST.ingredients.errorhandler import mastvasperrorhandlers 
from pymatgen.core.structure import Structure
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
        self.logger = loggerutils.get_mast_logger(self.keywords['name'])
  

    def get_error_handlers(self):
        """Get error handler list from custodian.vasp.handlers.
            This is a dynamic function because the Materials Project is
            continuously adding new handlers.
            Implement additional errors within the specific program type
            e.g. vasperror.py
        """
        handlerlist = list()
        #handlerlist = inspect.getmembers(handlers)
        handlerlist.extend(inspect.getmembers(masterrorhandlers))
        handlerlist.extend(inspect.getmembers(mastvasperrorhandlers))
        handlerdict=dict()
        for htry in handlerlist:
            hname = htry[0]
            if 'ErrorHandler' in hname and not (hname == 'ErrorHandler'):
                handlerdict[hname]=htry[1]
        return handlerdict
    
    def set_handler_inputs(self):
        """We may want to override the files that are being checked for
            each specific program type (vasp, phon, generic, etc.)
            Implement additional errors within the specific program type
            e.g. vasperror.py
        """
        handler_input_d=dict()
        handler_input_d['MASTMemoryErrorHandler']=[self.keywords['name'],self.keywords]
        return handler_input_d

    def loop_through_errors(self):
        """Loop through all errors in the error handlers.
        """
        mydir = self.keywords['name']
        errct = 0
        handlerdict = self.get_error_handlers() 
        handler_input_d = self.set_handler_inputs()
        myerror = ""
        os.chdir(mydir)
        self.logger.info("Checking errors for: %s" % mydir)
        for hname in handlerdict.keys():
            hinputs=""
            if hname in handler_input_d.keys():
                hinputs = handler_input_d[hname]
            else: #ignore errors whose inputs are not listed
                self.logger.warning("Skipping error check %s because no input was specified." % hname)
                continue 
            #if hinputs == "mast_skip":
            #    self.logger.info("Skipping %s" % hname)
            #    pass
            self.logger.info("Checking for %s with inputs %s." % (hname, hinputs))
            if len(hinputs) == 0:
                myerror = handlerdict[hname]()
            elif len(hinputs) == 1:
                myerror = handlerdict[hname](hinputs[0])
            elif len(hinputs) == 2:
                myerror = handlerdict[hname](hinputs[0],hinputs[1])
            elif len(hinputs) == 3:
                myerror = handlerdict[hname](hinputs[0],hinputs[1],hinputs[2])
            elif len(hinputs) == 4:
                myerror = handlerdict[hname](hinputs[0],hinputs[1],hinputs[2],hinputs[3])
            elif len(hinputs) == 5:
                myerror = handlerdict[hname](hinputs[0],hinputs[1],hinputs[2],hinputs[3],hinputs[4])
            else:
                raise MASTError(self.__class__.__name__, "Error %s has too many inputs (> 5)" % hname)
            if myerror.check():
                self.logger.error("%s Error found in directory %s!" % (hname, self.keywords['name']))
                errct = errct + 1
                if "mast_auto_correct" in self.keywords['program_keys'].keys():
                    if str(self.keywords['program_keys']['mast_auto_correct']).strip()[0].lower() == 'f':
                        self.logger.error("mast_auto_correct is false for this ingredient. Skipping correction.")
                        continue
                c_dict = myerror.correct()
                self.logger.error("Errors: %s" % c_dict["errors"])
                self.logger.error("Actions taken for %s: %s" % (self.keywords['name'],c_dict["actions"]))
            else:
                self.logger.info("%s No error found." % hname)
                pass
        return errct

    def clean_up_directory(self, alist=list()):
        """Clean up old output files. (Custodian handlers to not necessarily
            remove the old OUTCAR and OSZICAR files.)
            IMPORTANT: Any output to be archived must have been archived by
                the error handlers.
            Args:
                alist <list of str>: List of file names to remove
        """
        for aitem in alist:
            apath = os.path.join(self.keywords['name'], aitem)
            if os.path.isfile(apath):
                self.logger.info("Removing old file %s" % apath)
                os.remove(apath)
        return
