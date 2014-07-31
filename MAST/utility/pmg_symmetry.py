#!/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-07-31
##############################################################
import optparse
import pymatgen
import dirutil

def main():
    """Takes a POSCAR file and finds its symmetry based on
        Pymagen functions.
    """
    parser = optparse.OptionParser()
    parser.add_option('-t', '--tolerance', dest='tolerance', default='0.0001',
                      help="Tolerance for Pymatgen's symmetry-finding functions. Default 0.0001. Set to larger number for more variable structures.")
    parser.add_option('-i', '--input', dest='input', default=None,
                      help="POSCAR file name")
    
    parser.add_option('-d', '--debug', dest='debug', default="off",
                        help="Show debug information. Default 'off'. Set to 'on' to debug")
    
    (myopt, myarg) = parser.parse_args()
    if myopt.debug == "on":
        print myopt
        print myarg

    if myopt.input == None:
        print "An input file of VASP POSCAR-type is needed."
        print "Please use -i <POSCAR-type file> and try again."
        return
     
    print "---------------------------------------------------"
    print "Welcome to the MAterials Simulation Toolkit (MAST)"
    print "    Pymatgen symmetry-finder wrapper tool"
    print "Version: " + dirutil.get_version()
    print "Installed in: " + dirutil.get_mast_install_path()
    print "*** Pymatgen is produced by the Materials Project at"
    print "materialsproject.org." 
    print "---------------------------------------------------"
    
    mystructure=""

    try:
        mystructure = pymatgen.io.vaspio.Poscar.from_file(myopt.input).structure
    except ValueError:
        print "Structure could not be obtained from %s." % myopt.input
        print "Please check the file contents and try again."
        return
    if (mystructure == None) or (mystructure == ""):
        print "Structure could not be obtained from %s." % myopt.input
        print "Please check the file contents and try again."
        return

    mysf = pymatgen.symmetry.finder.SymmetryFinder(mystructure, float(myopt.tolerance))

    print "Tolerance used: ", myopt.tolerance
    print "Spacegroup:     ", mysf.get_spacegroup_symbol()
    print "Spacegroup No.: ", mysf.get_spacegroup_number()
    print "Hall No.:       ", mysf.get_hall()
    print "Point group:    ", mysf.get_point_group()

    return

if __name__ == '__main__':
    main()
