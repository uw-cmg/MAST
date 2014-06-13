##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Glen Jenness
# Last updated: 2013-07-03
##############################################################
import sys

"""Helper function for handling errors within MAST
   Modified from:
       http://www.doughellmann.com/articles/how-tos/python-exception-handling/index.html

   Example of use:
        error = 'File %s does not exist' % filename
        MASTerror(self.__class__.__name__, error)
"""

#def MASTError(classname, errormsg):
#    try:
#        raise RuntimeError('%s' % errormsg)
#        return 0
#    except Exception, err:
#        message = '\n\nError detected in MAST class %s!\n' % classname
#        message += 'Please fix the error and rerun mast.\n'
#        message += 'Error message:\n'
#        message += '    %s\n' % errormsg
#        message += '\n\n'
#
#        sys.stderr.write(message)
#        return 1

class MASTError(Exception):
    """Base class for exceptions in MAST."""
    def __init__(self, classname, errormsg):
        message = '\n\nError detected in MAST class %s!\n' % classname
        message += 'Please fix the error and rerun MAST.\n'
        message += 'Error message:\n'
        message += '    %s\n' % errormsg
        message += '\n\n'
        self.value = message
        sys.stderr.write(message)
    def __str__(self):
        return self.value

class FileNotFoundError(MASTError):
    """Exception raised for file not found in MAST."""

    def __init__(self, classname, filepath):
        errormsg = "File not found at " + filepath
        MASTError.__init__(self, classname, errormsg)

