############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
import os

import numpy as np
import pymatgen as pmg

from MAST.utility import InputOptions
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import MAST2Structure


ALLOWED_KEYS = {\
                 'optiondict'    : (InputOptions, None, 'Input file name'),\
               }

class InputPythonCreator(MASTObj):
    """Creates a runnable python file from a *.inp file
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.optdict=dict(self.keywords['optiondict'].options) 
        self.typeerror=0

    def print_input_options(self, lotsofspace=1):
        print "#MAST INPUT OPTIONS"

        for sectionname,secvalue in self.optdict.iteritems():
            if lotsofspace:
                print ""
                print "###################################"
            print "#Section: %3s" % sectionname
            if lotsofspace:
                print "####################################"
                print ""
            #print "TTM DEBUG: straight print---------------------"
            #print secvalue
            #print "TTM DEBUG: New print:------------------------------"
            self.print_multilevel_dict(secvalue)
        if self.typeerror > 0:
            raise MASTError(self.__class__.__name__,
                "Unsupported types detected in input file. Check .py file")
            
    def print_multilevel_dict(self, mydict):
        """Print a multi-level dictionary.
        """
        for key1, val1 in mydict.iteritems():
            if key1 == 'structure':
                print "#ignoring created structure"
                continue
            if not type(key1) == str:
                errorstr = "One of the keys in section names is not a string."
                raise MASTError(self.__class__.__name__, errorstr)
            if type(val1) == str:
                print key1 + " = " + "'" + val1 + "'"
            elif type(val1) == int:
                print key1 + " = ", val1
            elif type(val1) == float:
                print key1 + " = ", val1
            elif type(val1) == bool:
                print key1 + " = ", val1
            elif val1 == None:
                print key1 + " = None"
            elif type(val1) == np.ndarray:
                print key1 + " = " + "np.array(", val1.tolist(), ")"
            elif type(val1) == list:
                print key1 + " = list()"
                for myitem in val1:
                    if type(myitem) == str:
                        print key1 + ".append(" + "'" + myitem + "'" + ")"
                    else:
                        print key1 + ".append(" + myitem + ")"
            elif type(val1) == dict:
                print key1 + " = dict()"
                self.print_secondlevel_dict(key1, val1)
            else:
                print key1 + " UNSUPPORTED TYPE: ", type(val1)
                self.typeerror += 1
                print key1 + " = ", val1

    def print_secondlevel_dict(self, key1, seconddict):
        """Print a second-level dictionary. Key1 will be a string."""
        for key2, val2 in seconddict.iteritems():
            header=""
            if type(key2) == str:
                header = key1 + "[" + "'" + key2 + "'" + "]"
            else:
                header = key1 + "[" + key2 + "]"
            if type(val2) == str:
                print header + " = " + "'" + val2 + "'"
            elif type(val2) == int:
                print header + " = ", val2
            elif type(val2) == float:
                print header + " = ", val2
            elif type(val2) == bool:
                print header + " = ", val2
            elif val2 == None:
                print header + " = None"
            elif type(val2) == np.ndarray:
                print header + " = " + "np.array(", val2.tolist(), ")"
            elif type(val2) == list:
                print header + " = list()"
                for myitem in val2:
                    if type(myitem) == str:
                        print header + ".append(" + "'" + myitem + "'" + ")"
                    else:
                        print header + ".append(" + myitem + ")"
            elif type(val2) == dict:
                print header + " = dict()"
                self.print_thirdlevel_dict(header, key2, val2)
            else:
                print header + "UNSUPPORTED TYPE: ", type(val2)
                self.typeerror += 1
                print header + " = ", val2

    def print_thirdlevel_dict(self, oldheader, key2, thirddict):
        """Print a second-level dictionary. Key1 will be a string."""
        for key3, val3 in thirddict.iteritems():
            header=oldheader
            if type(key3) == str:
                header = header + "[" + "'" + key3 + "'" + "]"
            else:
                header = header + "[" + key3 + "]"
            if type(val3) == str:
                print header + " = " + "'" + val3 + "'"
            elif type(val3) == int:
                print header + " = ", val3
            elif type(val3) == float:
                print header + " = ", val3
            elif type(val3) == bool:
                print header + " = ", val3
            elif val3 == None:
                print header + " = None"
            elif type(val3) == np.ndarray:
                print header + " = " + "np.array(", val3.tolist(), ")"
            elif type(val3) == list:
                print header + " = list()"
                for myitem in val3:
                    if type(myitem) == str:
                        print header + ".append(" + "'" + myitem + "'" + ")"
                    else:
                        print header + ".append(" + myitem + ")"
            elif type(val3) == dict:
                print header + " = dict()"
                print header + " = " + val3
            else:
                print header + "UNSUPPORTED TYPE: ", type(val3)
                self.typeerror += 1
                print header + " = ", val3


