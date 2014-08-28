#!/usr/bin/env python
#TTM 2014-08-27
from pymatgen.matproj.rest import MPRester
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter
from pymatgen.io.vaspio import Poscar
import sys
import os

def main(api="", queryid=""):
    """Get VASP inputs for Materials Project structure
        Args:
            api <str>: Materials Project API key
            queryid <str>: Materials Project ID of the structure
        Returns:
            creates a folder named with that mpid, and including
            some VASP input files.
    """
    if api == "":
        print "Must have an API key from materialsproject.org"
        return None
    if queryid == "":
        print "No MP structure ID given. Exiting."
        return None
    rest_adapter = MPRester(api)
    entries=list()
    proplist=list()
    proplist.append('pretty_formula')
    proplist.append('structure')
    proplist.append('potcar')
    proplist.append('material_id')
    
    myentry = rest_adapter.mpquery(criteria={'material_id':queryid}, properties=proplist)
    if len(myentry) == 0:
        print "Could not find entry for %s as material_id. Trying entry_id." % queryid
        myentry = rest_adapter.mpquery(criteria={'entry_id':queryid}, properties=proplist)
    if len(myentry) == 0:
        print "Could not find entry for %s" % queryid
        return None
    entries.extend(myentry)

    workdir = os.getcwd()
    from pymatgen.io.vaspio_set import MITVaspInputSet, MPVaspInputSet
    for entry in entries: 
        mpvis = MPVaspInputSet()
        myname = str(entry['pretty_formula'])
        #print entry['structure'].composition
        #print entry['structure'].formula
        #myname = entry['pretty_formula']
        myname = myname.replace("(","_").replace(")","_")
        myname = myname + "_" + entry['material_id']
        os.mkdir(myname)
        os.chdir(myname)
        mystructure = entry['structure']
        if mystructure.num_sites <= 10:
            mystructure.make_supercell([2,2,2])
        #mystructure.perturb(0.01)
        incar = mpvis.get_incar(mystructure)
        incar.write_file("INCAR")
        potcar = mpvis.get_potcar(mystructure)
        potcar.write_file("POTCAR")
        #potcar_symbols = mpvis.get_potcar_symbols(mystructure) 
        myposcar=Poscar(mystructure)
        mykpoints=mpvis.get_kpoints(mystructure)
        mykpoints.write_file("KPOINTS")


        myposcar.write_file("POSCAR")
        os.chdir(workdir)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Not enough inputs."
        print main.__doc__
        sys.exit()
    main(sys.argv[1],sys.argv[2])
