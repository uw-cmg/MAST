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
        """
        handler_input_d=dict()
        handler_input_d['VaspErrorHandler']=["OUTCAR"]
        handler_input_d['DentetErrorHandler']=["OUTCAR"]
        handler_input_d['UnconvergedErrorHandler']=["vasprun.xml"]
        handler_input_d['PoscarErrorHandler']=["OUTCAR"]
        handler_input_d['TripleProductErrorHandler']=["OUTCAR"]
        handler_input_d['FrozenJobErrorHandler']="mast_skip"
        handler_input_d['NonConvergingErrorHandler']="mast_skip"
        handler_input_d['MeshSymmetryErrorHandler']="mast_skip"
        frozensec=21000
        if "mast_frozen_seconds" in self.keywords['program_keys'].keys():
            frozensec = int(self.keywords['program_keys']['mast_frozen_seconds'])
        handler_input_d['MASTFrozenJobErrorHandler']=["OUTCAR",frozensec,["OUTCAR","OSZICAR","CONTCAR","POSCAR"],["CONTCAR"],["POSCAR"]]

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
            if hinputs == "mast_skip":
                self.logger.info("Skipping %s" % hname)
                pass
            else:
                self.logger.info("Checking for %s" % hname)
                if len(hinputs) == 0:
                    myerror = handlerdict[hname]()
                elif len(hinputs) < 2:
                    myerror = handlerdict[hname](hinputs[0])
                elif len(hinputs) < 3:
                    myerror = handlerdict[hname](hinputs[0],hinputs[1])
                elif len(hinputs) < 4:
                    myerror = handlerdict[hname](hinputs[0],hinputs[1],hinputs[2])
                elif len(hinputs) < 5:
                    myerror = handlerdict[hname](hinputs[0],hinputs[1],hinputs[2],hinputs[3])
                elif len(hinputs) < 6:
                    myerror = handlerdict[hname](hinputs[0],hinputs[1],hinputs[2],hinputs[3],hinputs[4])
                else:
                    raise MASTError(self.__class__.__name__,"Error %s has too many inputs (more than 5)" % hname)
                if myerror.check():
                    self.logger.error("%s Error found in directory %s! Attempting to correct." % (hname, self.keywords['name']))
                    errct = errct + 1
                    if "mast_auto_correct" in self.keywords['program_keys'].keys():
                        if str(self.keywords['program_keys']['mast_auto_correct']).strip()[0].lower() == 'f':
                            errpath = os.path.join(os.path.dirname(mydir), "MAST_ERROR")
                            if not os.path.isfile(errpath):
                                errfile = MASTFile()
                            else:
                                errfile = MASTFile(errpath)
                            errfile.data.append("Error %s found in error handling for ingredient %s! Writing MAST_ERROR file for recipe.\n" % (hname, mydir))
                            errfile.to_file(errpath)
                        else:
                            c_dict = myerror.correct()
                            self.logger.error("Errors: %s" % c_dict["errors"])
                            self.logger.error("Actions taken for %s: %s" % (self.keywords['name'],c_dict["actions"]))
                    else: 
                        c_dict = myerror.correct()
                        self.logger.error("Errors: %s" % c_dict["errors"])
                        self.logger.error("Actions taken for %s: %s" % (self.keywords['name'],c_dict["actions"]))
                else:
                    self.logger.info("%s No error found." % hname)
                    pass
        return errct
