##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
from custodian.vasp import handlers
from pymatgen.core.structure import Structure
import logging
import inspect
import os
from MAST.ingredients.errorhandler import BaseError
from MAST.ingredients.errorhandler import masterrorhandlers 
from MAST.utility import MASTFile
#from MAST.utility import loggerutils
class VaspError(BaseError):
    """VASP error-handling functions (wraps custodian)
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
        """
        handlerdict = BaseError.get_error_handlers(self)
        handlerlist = inspect.getmembers(handlers)
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
        """
        handler_input_d=dict()
        handler_input_d['VaspErrorHandler']=["OUTCAR"]
        handler_input_d['DentetErrorHandler']=["OUTCAR"]
        handler_input_d['UnconvergedErrorHandler']=["mast_skip"]
        #handler_input_d['UnconvergedErrorHandler']=["vasprun.xml"]
        handler_input_d['PoscarErrorHandler']=["OUTCAR"]
        handler_input_d['TripleProductErrorHandler']=["OUTCAR"]
        #handler_input_d['FrozenJobErrorHandler']=["mast_skip"]
        #handler_input_d['NonConvergingErrorHandler']=["mast_skip"]
        #handler_input_d['MeshSymmetryErrorHandler']=["mast_skip"]
        handler_input_d['MASTWalltimeErrorHandler']=[self.keywords['name'],["OUTCAR","OSZICAR","CONTCAR","POSCAR","XDATCAR"],["CONTCAR"],["POSCAR"]]
        handler_input_d['MASTMemoryErrorHandler']=[self.keywords['name'],self.keywords,["OUTCAR","OSZICAR","CONTCAR","XDATCAR"]]
        handler_input_d['VaspReachedNSWErrorHandler']=[self.keywords['name'],["OUTCAR","OSZICAR","CONTCAR","POSCAR","XDATCAR"],["CONTCAR"],["POSCAR"]]

        return handler_input_d

    def clean_up_directory(self):
        """Clean up directory."""
        return BaseError.clean_up_directory(self, ["OUTCAR","OSZICAR"])
