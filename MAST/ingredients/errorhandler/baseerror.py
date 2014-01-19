import os
import time
import shutil
import logging
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import Metadata
from MAST.utility import loggerutils
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
        self.logger = logging.getLogger(self.keywords['name'])
        self.logger = loggerutils.add_handler_for_recipe(self.keywords['name'], self.logger)
  

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
        handler_input_d['MASTMemoryErrorHandler']=[self.keywords['name']]
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
            self.logger.info("Checking for %s" % hname)
            if len(hinputs) == 0:
                myerror = handlerdict[hname]()
            else:
                myerror = handlerdict[hname](hinputs)
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

