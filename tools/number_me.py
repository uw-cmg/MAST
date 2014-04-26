#!/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
#TTM 13Mar2014
def doc():
    """This routine numbers all the lines in a text file
        and produces a new text file.
        e.g. 
            apple
            banana
            carrot
        becomes
            1.apple
            2.banana
            3.carrot
        Call as: python number_me.py original_file file_to_produce
        
        python number_me.py unnumbered numbered
    """
import os
import sys

if len(sys.argv) < 3:
    print "Not enough inputs."
    print doc.__doc__
    sys.exit(-1)

myname = sys.argv[1]
if not os.path.isfile(myname):
    print "No such file %s" % myname
    sys.exit(-1)
newname = sys.argv[2]
if os.path.isfile(newname):
    print "Cannot overwrite existing file %s" % newname
    sys.exit(-1)

openfile = open(myname, 'rb')
mylines = openfile.readlines()
openfile.close()
lct=0
numberedlines=list()
for myline in mylines:
    lct = lct + 1
    numline = "%i.%s" % (lct, myline)
    numberedlines.append(numline)
newfile = open(newname, 'wb')
newfile.writelines(numberedlines)
newfile.close()
print "Done! Check output in %s" % newname
