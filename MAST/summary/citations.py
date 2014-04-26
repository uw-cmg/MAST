##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import sys
import time
from MAST.utility import fileutil
citation_dict=dict()
citation_dict['mast']=list()
citation_dict['mast'].append("MAterials Simulation Toolkit (MAST). Developed at the University of Wisconsin-Madison under NSF award number 1148011. First use cited in T. Angsten, T. Mayeshiba, H. Wu, and D. Morgan, New J. Phys. 16, 015018 (2014). doi:10.1088/1367-2630/16/1/015018")
citation_dict['vasp']=list()
citation_dict['vasp'].append("Vienna Ab initio Simulation Package (VASP). https://www.vasp.at")
citation_dict['vasp'].append("G. Kresse and J. Hafner, Phys. Rev. B 47 , 558 (1993).")
citation_dict['vasp'].append("G. Kresse and J. Hafner, Phys. Rev. B 49 , 14 251 (1994).")
citation_dict['vasp'].append("G. Kresse and J. Furthmueller, Comput. Mat. Sci. 6 , 15 (1996).")
citation_dict['vasp'].append("G. Kresse and J. Furthmueller, Phys. Rev. B 54 , 11 169 (1996).")
citation_dict['vasp_pps']=list()
citation_dict['vasp_pps'].append("G. Kresse and J. Hafner, J. Phys.: Condens. Matt. 6, 8245 (1994)")
citation_dict['vasp_paw']=list()
citation_dict['vasp_paw'].append("G. Kresse and D. Joubert, Phys. Rev. 59 , 1758 (1999).")
citation_dict['phon']=list()
citation_dict['phon'].append("PHON. D. Alfe`, Comput. Phys. Commun. 180, 2622.2633 (2009). http://www.homepages.ucl.ac.uk/~ucfbdxa/phon/")
citation_dict['pymatgen']=list()
citation_dict['pymatgen'].append("Pymatgen. http://pymatgen.org")
citation_dict['pymatgen'].append("Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier, Michael Kocher, Shreyas Cholia, Dan Gunter, Vincent Chevrier, Kristin A. Persson, Gerbrand Ceder. Python Materials Genomics (pymatgen) : A Robust, Open-Source Python Library for Materials Analysis. Computational Materials Science, 2013, 68, 314.319. doi:10.1016/j.commatsci.2012.10.028")
citation_dict['spglib']=list()
citation_dict['spglib'].append("Spglib. Atsushi Togo, http://spglib.sourceforge.net")
citation_dict['spglib'].append("Spglib project acknowledges Yusuke Seto for the Crystallographic database and Dimitar Pashov for the fortran interface.")
###More pymatgen citations may be added depending on which modules are additionally used by MAST.
citation_dict['python']=list()
citation_dict['python'].append("Python. http://www.python.org")
citation_dict['structopt']=list()
citation_dict['structopt'].append("The 'structopt' package for implementing genetic algorithms on clusters was created 2012-2014 at the University of Wisconsin-Madison by Amy Kaczmarowski and tested by Hyunseok Ko and Min Yu, under the research direction of Professor Dane Morgan.")
def get_citations(recipedir):
    """Use mast_exec lines from the recipe directory in order
        to assign a CITATIONS section to the SUMMARY.txt file.
        Args:
            recipedir <str>: Recipe directory
        Returns:
            linelist <list of str>: Lines to add to SUMMARY.txt
    """
    mylist=list()
    mylist=['mast','pymatgen','spglib','python'] #need to use
    myinput = os.path.join(recipedir, "input.inp")
    mastexecs = fileutil.grepme(myinput, "mast_exec")
    linelist=list()
    for mastexec in mastexecs:
        if 'vasp' in mastexec.lower():
            mylist.append('vasp')
            mastxcs = fileutil.grepme(myinput, "mast_xc")
            if len(mastxcs) > 1:
                mylist.append('vasp_pps')
                for mastxc in mastxcs:
                    if ("pw" in mastxc.lower()) or ("paw" in mastxc.lower()):
                        mylist.append('vasp_paw')
        if 'phon' in mastexec.lower():
            mylist.append('phon')
        if 'structopt' in mastexec.lower():
            mylist.append('structopt')
    linelist.append("==================================")
    linelist.append("----------- CITATIONS-------------")
    linelist.append("==================================")
    cited=list()
    for listitem in mylist:
        if not listitem in cited:
            if not listitem in citation_dict.keys():
                linelist.append("NEED TO CITE %s BUT NO ENTRY IS IN CITATION DICT.")
            else:
                linelist.append('======================')
                linelist.append("    %s" % listitem)
                linelist.append('======================')
                for citeline in citation_dict[listitem]:
                    linelist.append(citeline)
            cited.append(listitem)
    return linelist
