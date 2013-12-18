from custodian.vasp import handlers
from pymatgen.core.structure import Structure
import logging
import inspect
import os
from MAST.ingredients.errorhandler import BaseError
from MAST.ingredients.errorhandler import VaspError
from MAST.ingredients.errorhandler import masterrorhandlers
class VaspNEBError(BaseError):
    """VASP NEB error-handling functions (wraps custodian)
    """
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseError.__init__(self, allowed_keys, **kwargs)


    def get_error_handlers(self):
        """Get error handler list from custodian.vasp.handlers.
            This is a dynamic function because the Materials Project is
            continuously adding new handlers.
            This function may change when NEB-specific error handlers
            become available (see MIT)
        """
        raise NotImplementedError()
        handlerlist = inspect.getmembers(handlers)
        handlerlist.extend(inspect.getmembers(masterrorhandlers))
        handlerdict=dict()
        for htry in handlerlist:
            hname = htry[0]
            if 'ErrorHandler' in hname and not (hname == 'ErrorHandler'):
                handlerdict[hname]=htry[1]
        return handlerdict


    def set_handler_inputs(self):
        """We may want to override the files that are being checked, which are
            usually vasp.out or vasprun.xml.
            Use key as error handler class name (like 'VaspErrorHandler').
            Use mast_skip to skip the error.
            This function may change when NEB-specific error handlers
            become available (see MIT)
        """
        raise NotImplementedError()
        handler_input_d=dict()
        handler_input_d['MASTFrozenJobErrorHandler']=["OUTCAR",1000,["OUTCAR","OSZICAR","CONTCAR","POSCAR"],["CONTCAR"],["POSCAR"]]

        return handler_input_d

    def loop_through_errors(self):
        """Loop through all errors in the error handlers.
            Currently uses only VASP single-run error handlers,
            per image.
        """
        return 0 #VaspErrors are not meant to be in image folders.
        mydir = self.keywords['name']
        errct = 0
        #handlerdict = self.get_error_handlers() 
        #handler_input_d = self.set_handler_inputs()
        self.logger.info("Checking errors for: %s" % mydir)
        imct = 1
        while imct <= self.keywords['program_keys']['mast_neb_settings']['images']:
            imstr = str(imct).zfill(2)
            self.logger.info("Checking errors for: %s" % imstr)
            impath = os.path.join(mydir, imstr)
            imhandler = VaspError(name=impath, program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])
            imerrct = imhandler.loop_through_errors()
            errct = errct + imerrct
            imct = imct + 1
        return errct
