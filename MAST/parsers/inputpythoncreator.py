##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# This file is obsolete. 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os

import numpy as np
import pymatgen as pmg

from MAST.utility import InputOptions
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import MAST2Structure
from MAST.utility.mastfile import MASTFile

ALLOWED_KEYS = {\
                 'input_options'    : (InputOptions, None, 'Input file name'),\
                 'filename'         : (str, "input.py", 'Script file name'),\
               }

class InputPythonCreator(MASTObj):
    """Creates a runnable python file from a *.inp file
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
    
    def write_script(self, dirname="", fname=""):
        """Write the python input script, created from the *.inp input file.
            Args:
                dirname <str>: directory in which to write
                fname <str>: filename to write
            Returns:
                filename <str>: Full file name of the input script created.
        """
        mylines = self.print_input_options()
        if dirname == "":
            dirname = os.getcwd()
        if fname == "":
            fname = 'input.py'
        filename = os.path.join(dirname, fname)
        inputpy = MASTFile()
        inputpy.data = mylines
        inputpy.to_file(filename)
        return filename


    def addnewlines(self, linestoadd):
        """Appends newline to each string in a list
            Args:
                linestoadd <list of str>: lines to add newline to.
        """
        newline = "\n"
        newlist=list()
        for myline in linestoadd:
            try:
                myline = myline + newline
            except:
                raise MASTError(self.__class__.__name__,"Error with " + str(myline))
            newlist.append(myline)
        return newlist

    def print_input_options(self):
        inputoptions = self.keywords['input_options']
        pln=list()
        pln.append("from MAST.utility import InputOptions")
        pln.append("import numpy as np")
        pln.append("#MAST INPUT OPTIONS")
        pln.append("inputoptions = InputOptions()")

        for sectionname,secvalue in inputoptions.options.iteritems():
            pln.append(" ")
            pln.append("###################################")
            pln.append("#Section: %3s" % sectionname)
            pln.append("####################################")
            pln.append(" ")
            pln.extend(inputoptions.print_python_section(sectionname))
            pln.append("inputoptions.options['"+sectionname+"'] = " + sectionname)
        pln.extend(self.print_create_structure_section())
        pln.extend(self.print_mast_command_section())
        return self.addnewlines(pln)

    def print_create_structure_section(self):
        """Print the structure creation section.
            This is a little complex and may not belong here.

            Returns:
                pln <list of str>: list of lines to print out
        """
        pln=list()
        pln.append(" ")
        pln.append("###############################################")
        pln.append("#Structure Creation Section")
        pln.append("################################################")
        pln.append(" ")
        iops = self.keywords['input_options']
        strposfile = iops.get_item('structure','posfile')
        if strposfile is None:
            pln.append("iopscoords=inputoptions.get_item('structure','coordinates')")
            pln.append("iopslatt=inputoptions.get_item('structure','lattice')")
            pln.append("iopsatoms=inputoptions.get_item('structure','atom_list')")
            pln.append("iopsctype=inputoptions.get_item('structure','coord_type')")
            pln.append("from MAST.utility import MAST2Structure")
            pln.append("structure = MAST2Structure(lattice=iopslatt,")
            pln.append("    coordinates=iopscoords, atom_list=iopsatoms,")
            pln.append("    coord_type=iopsctype)")
        elif ('poscar' in strposfile.lower()):
            pln.append("from pymatgen.io.vaspio import Poscar")
            pln.append("structure = Poscar.from_file('" +strposfile+ "').structure")
        elif ('cif' in strposfile.lower()):
            #from pymatgen.io.cifio import CifParser
            #structure = CifParser(strposfile).get_structures()[0]
            pln.append("from pymatgen.io.cifio import CifParser")
            pln.append("structure = CifParser('"+strposfile+"').get_structures()[0]")
        else:
            error = 'Cannot build structure from file %s' % strposfile
            raise MASTError(self.__class__.__name__, error)
        pln.append("inputoptions.update_item('structure','structure',structure)")
        return pln

    def print_mast_command_section(self):
        """Print the mast commands section.
            Returns:
                pln <list of str>: list of lines to append
        """
        pln=list()
        pln.append(" ")
        pln.append("###############################################")
        pln.append("#MAST Command Section")
        pln.append("################################################")
        pln.append(" ")
        inputoptions = self.keywords['input_options']
        recipe_base_name = os.path.basename(inputoptions.options['recipe']['recipe_file'])
        recipe_base_name = recipe_base_name.split('.')[0]
        import time
        timestamp = time.strftime('%Y%m%d%H%M%S')
        self.name = recipe_base_name + '_' + timestamp + '.py'
        
        pln.append("from MAST.mast import MAST")
        pln.append("mast_obj = MAST()")
        pln.append("mast_obj.start_from_input_options(inputoptions)")
        return pln

