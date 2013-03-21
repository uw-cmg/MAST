import sys

"""Helper function for handling errors within MAST
   Modified from:
       http://www.doughellmann.com/articles/how-tos/python-exception-handling/index.html

   Example of use:
        error = 'File %s does not exist' % filename
        MASTerror(self.__class__.__name__, error)
"""

def MASTError(classname, errormsg):
    try:
        raise RuntimeError('%s' % errormsg)
        return 0
    except Exception, err:
        message = '\n\nError detected in MAST class %s!\n' % classname
        message += 'Please fix the error and rerun mast.\n'
        message += 'Error message:\n'
        message += '    %s\n' % errormsg
        message += '\n\n'

        sys.stderr.write(message)
        return 1

