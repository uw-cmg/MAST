#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
#TTM 2014-03-27
#use pymatgen.phasediagram.pdanalyzer get_decomposition(pymatgen.core.composition Composition, and using an open element O2 quaternary phase diagram) to get compositions and percentages for some composition (e.g. LaMn0.4Fe0.6O3)
from pymatgen.matproj.rest import MPRester
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.pdmaker import GrandPotentialPhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter
from pymatgen.phasediagram.pdanalyzer import PDAnalyzer
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element

def get_decomp(o_chem_pot, mycomp, verbose=1):
    """Get decomposition from open phase diagram
        Args:
            o_chem_pot <float>: Oxygen chemical potential
            mycomp <pymatgen Composition>: Composition
            verbose <int>: 1 - verbose (default)
                           0 - silent
        Returns:
            decomposition string
    """        
    a = MPRester("<YOUR_MPREST_API_KEY_HERE>")
    elements = mycomp.elements
    ellist = map(str, elements)
    entries = a.get_entries_in_chemsys(ellist)
    #entries = a.get_entries_in_chemsys(['La', 'Mn', 'O', 'Fe'])
    pd = PhaseDiagram(entries)
    gppd = GrandPotentialPhaseDiagram(entries,{Element('O'): float(o_chem_pot)})
    print gppd
    #plotter = PDPlotter(gppd)
    #plotter.show()

    gppda = PDAnalyzer(gppd)
    #mychempots = gppda.get_composition_chempots(mycomp)
    #print "My chem pots:"
    #print mychempots
    mydecompgppd = gppda.get_decomposition(mycomp)
    #pdentry = PDEntry(mycomp, 0)
    #print "Decomp and energy:"
    #decompandenergy = gppda.get_decomp_and_e_above_hull(pdentry)
    #print decompandenergy
    #mydecomppd = pda.get_decomposition(mycomp)
    #print "Mn profile:"
    #mnprof= gppda.get_element_profile(Element('Mn'),mycomp)
    #print mnprof

    if verbose:
        for (entry,amount) in mydecompgppd.iteritems():
            print "%s: %3.3f" % (entry.name, amount)
            #mymurangegppd = gppda.getmu_range_stability_phase(Composition(entry.name),Element('O'))
            #print mymurangegppd
        #for (entry,amount) in mydecomppd.iteritems():
        #    print "%s: %3.3f" % (entry.name, amount)
        print ""
    return mydecompgppd

mycomp = Composition(La=1,Mn=0.4,Fe=0.6,O=3)
import numpy as np
import sys
myrange = np.arange(-12,1,0.25)

for ocp in myrange:
    get_decomp(ocp, mycomp)
