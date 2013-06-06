from custodian.vasp import handlers
import inspect
import os

def get_error_handlers():
    """Get error handler list from custodian.vasp.handlers.
        This is a dynamic function because the Materials Project is
        continuously adding new handlers.
    """
    handlerlist = inspect.getmembers(handlers)
    handlerdict=dict()
    for htry in handlerlist:
        hname = htry[0]
        if 'ErrorHandler' in hname and not (hname == 'ErrorHandler'):
            handlerdict[hname]=htry[1]
    return handlerdict


def set_handler_inputs():
    """We may want to override the files that are being checked, which are
        usually vasp.out or vasprun.xml.
        Use key as error handler class name (like 'VaspErrorHandler').
        Use mast_skip to skip the error.
    """
    handler_input_d=dict()
    handler_input_d['VaspErrorHandler']="OUTCAR"
    handler_input_d['DentetErrorHandler']="OUTCAR"
    handler_input_d['UnconvergedErrorHandler']="vasprun.xml,7200"
    handler_input_d['PoscarErrorHandler']="OUTCAR"
    handler_input_d['TripleProductErrorHandler']="OUTCAR"
    handler_input_d['FrozenJobErrorHandler']="OUTCAR"
    handler_input_d['NonConvergingErrorHandler']="mast_skip"
    handler_input_d['MeshSymmetryErrorHandler']="mast_skip"

    return handler_input_d

def loop_through_errors(mydir):
    """Loop through all errors in the error handlers.
    """
    errct = 0
    handlerdict = get_error_handlers() 
    handler_input_d = set_handler_inputs()
    myerror = ""
    os.chdir(mydir)
    print "Checking errors for: ",mydir
    for hname in handlerdict.keys():
        hinputs=""
        if hname in handler_input_d.keys():
            hinputs = handler_input_d[hname]
        if hinputs == "mast_skip":
            print "Skipping ",hname
        else:
            print "Checking for ",hname
            myerror = handlerdict[hname](hinputs)
            if myerror.check():
                print hname, "Error found! Attempting to correct."
                errct = errct + 1
                myerror.correct()
            else:
                print hname,": No error found."
    return errct
