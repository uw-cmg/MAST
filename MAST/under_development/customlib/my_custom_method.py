##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import sys

def main(arg1="No argument given"):
    """This is a custom method."""
    returnme = "%s example" % arg1
    return "Return value: %s" % returnme

if __name__ == '__main__':
    #This function is meant to be called by the recipe plan,
    #    as part of the mast_xxx_method keywords.
    #In addition to inputs specified in the input file,
    #    the method will also get the ingredient directory, 
    #    as the first input.
    #If the method is called from the mast_update_children_method
    #    keyword, it will get the child directory as the 
    #    last input.
    #For example:
    #    mast_write_method write_some_function.py bulk.xyz
    #    mast_update_children_method myfunction.py eval.txt
    #The first example will have sys.argv as follows:
    #    sys.argv[1] = <string, of ingredient directory>
    #    sys.argv[2] = "bulk.xyz" <string>
    #The second example will have sys.argv as follows:
    #    sys.argv[1] = <string, of ingredient directory>
    #    sys.argv[2] = "eval.txt" <string>
    #    sys.argv[3] = <string, of child directory>
    if len(sys.argv) > 1:
        arg1 = sys.argv[1]
        myreturn = main(arg1)
    else:
        myreturn = main()
    print myreturn
    sys.exit()
