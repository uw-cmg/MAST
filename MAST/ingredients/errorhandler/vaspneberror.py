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


    def set_handler_inputs(self):
        """We may want to override the files that are being checked, which are
            usually vasp.out or vasprun.xml.
            Use key as error handler class name (like 'VaspErrorHandler').
            Use mast_skip to skip the error.
        """
        handler_input_d=dict()
        numim = self.keywords['program_keys']['mast_neb_settings']['images']
        numim = int(numim)
        archlist = list()
        clist = list()
        plist = list()
        for imct in range(1, numim+1):
            imstr = str(imct).zfill(2)
            archlist.append("%s/OUTCAR" % imstr)
            archlist.append("%s/OSZICAR" % imstr)
            archlist.append("%s/CONTCAR" % imstr)
            archlist.append("%s/POSCAR" % imstr)
            archlist.append("%s/XDATCAR" % imstr)
            clist.append("%s/CONTCAR" % imstr)
            plist.append("%s/POSCAR" % imstr)
        handler_input_d['MASTWalltimeErrorHandler']=[self.keywords['name'],archlist, clist, plist]
        handler_input_d['MASTMemoryErrorHandler']=[self.keywords['name'],archlist]

        return handler_input_d

