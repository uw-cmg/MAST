#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba, Ryan Jacobs
# Last updated: 2015-01-14
#
##############################################################
#TTM 2014-03-27
#use pymatgen.phasediagram.pdanalyzer get_decomposition(pymatgen.core.composition Composition, and using an open element O2 quaternary phase diagram) to get compositions and percentages for some composition (e.g. LaMn0.4Fe0.6O3)
from pymatgen.matproj.rest import MPRester
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.pdmaker import GrandPotentialPhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter
from pymatgen.phasediagram.pdanalyzer import PDAnalyzer
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.entries import GrandPotPDEntry
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
import numpy as np

def main(comp="La0.5Sr0.5MnO3", energy=-43.3610, ostart="", oend="", ostep=""):
    """Get energy above hull for a composition
        Args:
            comp <str>: Composition in string form
            energy <float>: Energy PER FORMULA UNIT of composition given
            (Leave the following arguments blank for a non-grand potential
                phase diagram.)
            ostart <float>: Starting oxygen chemical potential. 
            oend <float>: Ending oxygen chemical potential. 
            ostep <float>: Step for oxygen chemical potential
        Returns:
            Prints to screen
    """        
    #a = MPRester("<YOUR_MPREST_API_KEY_HERE>")
    a = MPRester("wfmUu5VSsDCvIrhz")
    
    mycomp=Composition(comp)
    print "Composition: ", mycomp
    myenergy=energy
    print "Energy: ", myenergy
    myPDEntry = PDEntry(mycomp, myenergy)

    elements = mycomp.elements
    ellist = map(str, elements)
    
    chemsys_entries = a.get_entries_in_chemsys(ellist)
    #For reference: other ways of getting entries
    #entries = a.mpquery(criteria={'elements':{'$in':['La','Mn'],'$all':['O']},'nelements':3})
    #entries = a.mpquery(criteria={'elements':{'$in':['La','Mn','O'],'$all':['O']}},properties=['pretty_formula'])
    #entries = a.get_entries_in_chemsys(['La', 'Mn', 'O', 'Sr'])
   
    if ostart=="": #Regular phase diagram
        entries = list(chemsys_entries)
        entries.append(myPDEntry)
        pd = PhaseDiagram(entries)
        #plotter = PDPlotter(gppd)
        #plotter.show()
        ppda = PDAnalyzer(pd)
        eabove=ppda.get_decomp_and_e_above_hull(myPDEntry)
        print "Energy above hull: ", eabove[1]
        print "Decomposition: ", eabove[0]
        return eabove
    else: #Grand potential phase diagram
        orange = np.arange(ostart, oend+ostep, ostep) #add ostep because otherwise the range ends before oend
        for o_chem_pot in orange:
            entries = list(chemsys_entries)
            myGrandPDEntry = GrandPotPDEntry(myPDEntry,{Element('O'): float(o_chem_pot)}) #need grand pot pd entry for GPPD
            entries.append(myGrandPDEntry)
            gppd = GrandPotentialPhaseDiagram(entries,{Element('O'): float(o_chem_pot)})
            gppda = PDAnalyzer(gppd)
            geabove=gppda.get_decomp_and_e_above_hull(myGrandPDEntry, True)
            print "******** Decomposition for mu_O = %s eV ********" % o_chem_pot
            print "%30s%1.4f" % ("mu_O: ",o_chem_pot)
            print "%30s%1.4f" % ("Energy above hull (eV): ",geabove[1])
            decomp=geabove[0]
            #print "Decomp: ", decomp
            print "%30s" % "Decomposition: "
            for dkey in decomp.keys():
                print "%30s:%1.4f" % (dkey.composition,decomp[dkey])
    return

if __name__=="__main__":
    import sys
    lensys = len(sys.argv)
    if lensys < 2:
        comp="La0.5Sr0.5MnO3"
    else:
        comp=sys.argv[1]
    if lensys < 3:
        energy=-43.3610
    else:
        energy=float(sys.argv[2])
    if lensys < 4:
        ostart=""
    else:
        ostart=float(sys.argv[3])
    if lensys < 5:
        oend=""
    else:
        oend=float(sys.argv[4])
    if lensys < 6:
        ostep=""
    else:
        ostep=float(sys.argv[5])
    main(comp, energy, ostart, oend, ostep)
