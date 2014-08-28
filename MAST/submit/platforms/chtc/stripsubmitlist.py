#!/bin/env python
############
# TTM 2014-08-28 nshd_mast
############
import os
import sys
#import optparse
#from MAST.utility.mastfile import MASTFile
#from MAST.utility import dirutil
#from MAST.recipe.recipeplan import RecipePlan

def main():
    """Strip $TEMP directory out of the submission list."""
    homedir=os.getenv("HOME")
    subname="%s/submitlist" % os.getenv("MAST_CONTROL")
    sublist=open(subname,"rb")
    mylines = sublist.readlines()
    sublist.close()
    newlines=list()
    newline=""
    for myline in mylines:
        newline=""
        if "/var/lib/" in myline:
            mysplit = myline.split("/")
            ingname = mysplit[-1]
            recipename = mysplit[-2]
            newline = "%s/MAST/SCRATCH/%s/%s\n" % (homedir, recipename, ingname)
            newlines.append(newline)
        else:
            newlines.append(myline)
    sublist=open(subname,"wb")
    sublist.writelines(newlines)
    sublist.close()
    return

if __name__ == "__main__":
    main()
