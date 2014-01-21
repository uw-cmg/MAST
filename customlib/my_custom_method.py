#!/usr/bin/env python
import os
import sys

def main(arg1="No argument given"):
    """This is a custom method."""
    returnme = "%s example" % arg1
    return "Return value: %s" % returnme

if __name__ == '__main__':
    if len(sys.argv) > 1:
        arg1 = sys.argv[1]
        myreturn = main(arg1)
    else:
        myreturn = main()
    print myreturn
    sys.exit()
