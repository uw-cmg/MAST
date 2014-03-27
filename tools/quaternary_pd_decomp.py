#!/usr/bin/env python
#TTM 2014-03-27
#use pymatgen.phasediagram.pdanalyzer get_decomposition(pymatgen.core.composition Composition, and using an open element O2 quaternary phase diagram) to get compositions and percentages for some composition (e.g. LaMn0.4Fe0.6O3)
from pymatgen.matproj.rest import MPRester
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.pdmaker import GrandPotentialPhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter
from pymatgen.phasediagram.pdanalyzer import PDAnalyzer
from pymatgen.core.composition import Composition

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
    a = MPRester("przN3JsK1grF3Q0w")
    elements = mycomp.elements
    ellist = map(str, elements)
    entries = a.get_entries_in_chemsys(ellist)
    #entries = a.get_entries_in_chemsys(['La', 'Mn', 'O', 'Fe'])
    pd = PhaseDiagram(entries)
    gppd = GrandPotentialPhaseDiagram(entries,{'O': float(o_chem_pot)})
    print gppd
    #plotter = PDPlotter(gppd)
    #plotter.show()

    gppda = PDAnalyzer(gppd)
    mychempots = gppda.get_composition_chempots(mycomp)
    print mychempots
    mydecompgppd = gppda.get_decomposition(mycomp)
    #mymurangegppd = gppda.getmu_range_stability_phase(mycomp, 'O')
    #print mymurangegppd
    pda = PDAnalyzer(pd)
    #mymurangepd = pda.getmu_range_stability_phase(mycomp, 'O')
    #print mymurangepd
    mydecomppd = pda.get_decomposition(mycomp)
    if verbose:
        for (entry,amount) in mydecompgppd.iteritems():
            print "%s: %3.3f" % (entry.name, amount)
        for (entry,amount) in mydecomppd.iteritems():
            print "%s: %3.3f" % (entry.name, amount)
        print ""
    return mydecompgppd

mycomp = Composition(La=1,Mn=0.4,Fe=0.6,O=3)
for ocp in [-10, -5, 0, 5, 10]:
    get_decomp(ocp, mycomp)
