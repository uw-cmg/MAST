############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
from __future__ import print_function
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
        if not 'mydir' in self.keywords.keys():
            mydir = os.getcwd()
        else:
            mydir == self.keywords['mydir']
        self.myfile = open(os.path.join(mydir, "input.py"), "wb")

    def print_input_options(self, lotsofspace=1):
        print("#MAST INPUT OPTIONS", file=self.myfile)

        for sectionname,secvalue in self.optdict.iteritems():
            if lotsofspace:
                print("", file=self.myfile)
                print("###################################", file=self.myfile)
            print("#Section: %3s" % sectionname, file=self.myfile)
            if lotsofspace:
                print("####################################", file=self.myfile)
                print("", file=self.myfile)
            #print "TTM DEBUG: straight print---------------------"
            #print secvalue
            #print "TTM DEBUG: New print:------------------------------"
            self.print_multilevel_dict(secvalue)
        self.print_buffet_section()
        self.myfile.close()
        if self.typeerror > 0:
            raise MASTError(self.__class__.__name__,
                "Unsupported types detected in input file. Check .py file")
            
    def print_multilevel_dict(self, mydict):
        """Print a multi-level dictionary.
        """
        for key1, val1 in mydict.iteritems():
            if key1 == 'structure':
                print("#ignoring created structure", file=self.myfile)
                continue
            if not type(key1) == str:
                errorstr = "One of the keys in section names is not a string."
                raise MASTError(self.__class__.__name__, errorstr)
            if type(val1) == str:
                mystr = key1 + " = " + "'" + val1 + "'"
                print(mystr, file=self.myfile)
            elif type(val1) in [int, bool, float]:
                mystr = key1 + " = ", val1
                print(key1,"=",val1, file=self.myfile)
            elif val1 == None:
                mystr = key1 + " = None"
                print(mystr, file=self.myfile)
            elif type(val1) == np.ndarray:
                print(key1,"=","np.array(",val1.tolist(),")",file=self.myfile)
            elif type(val1) == list:
                mystr = key1 + " = list()"
                print(mystr, file=self.myfile)
                for myitem in val1:
                    if type(myitem) == str:
                        mystr = key1 + ".append(" + "'" + myitem + "'" + ")"
                    else:
                        mystr = key1 + ".append(" + myitem + ")"
                    print(mystr, file=self.myfile)
            elif type(val1) == dict:
                mystr = key1 + " = dict()"
                print(mystr, file=self.myfile)
                self.print_secondlevel_dict(key1, val1)
            else:
                print(key1,"UNSUPPORTED TYPE",type(val1), file=self.myfile)
                self.typeerror += 1
                print(key1,"=",val1, file=self.myfile)

    def print_secondlevel_dict(self, key1, seconddict):
        """Print a second-level dictionary. Key1 will be a string."""
        for key2, val2 in seconddict.iteritems():
            header=""
            if type(key2) == str:
                header = key1 + "[" + "'" + key2 + "'" + "]"
            else:
                header = key1 + "[" + key2 + "]"
            if type(val2) == str:
                mystr = header + " = " + "'" + val2 + "'"
                print(mystr, file=self.myfile)
            elif type(val2) in [int,float,bool]:
                print(header,"=",val2, file=self.myfile)
            elif val2 == None:
                mystr = header + " = None"
                print(mystr, file=self.myfile)
            elif type(val2) == np.ndarray:
                print(header,"=","np.array(",val2.tolist(),")",file=self.myfile)
            elif type(val2) == list:
                mystr = header + " = list()"
                print(mystr, file=self.myfile)
                for myitem in val2:
                    if type(myitem) == str:
                        mystr = header + ".append(" + "'" + myitem + "'" + ")"
                    else:
                        mystr = header + ".append(" + myitem + ")"
                    print(mystr, file=self.myfile)
            elif type(val2) == dict:
                mystr = header + " = dict()"
                print(mystr, file=self.myfile)
                self.print_thirdlevel_dict(header, key2, val2)
            else:
                print(header,"UNSUPPORTED TYPE",type(val2), file=self.myfile)
                self.typeerror += 1
                print(header,"=",val2, file=self.myfile)

    def print_thirdlevel_dict(self, oldheader, key2, thirddict):
        """Print a second-level dictionary. Key1 will be a string."""
        for key3, val3 in thirddict.iteritems():
            header=oldheader
            if type(key3) == str:
                header = header + "[" + "'" + key3 + "'" + "]"
            else:
                header = header + "[" + key3 + "]"
            if type(val3) == str:
                mystr = header + " = " + "'" + val3 + "'"
                print(mystr, file=self.myfile)
            elif type(val3) in [int, float, bool]:
                print(header,"=",val3, file=self.myfile)
            elif val3 == None:
                mystr = header + " = None"
                print(mystr, file=self.myfile)
            elif type(val3) == np.ndarray:
                print(header,"=","np.array(",val3.tolist(),")",file=self.myfile)
            elif type(val3) == list:
                mystr = header + " = list()"
                print(mystr, file=self.myfile)
                for myitem in val3:
                    if type(myitem) == str:
                        mystr = header + ".append(" + "'" + myitem + "'" + ")"
                    else:
                        mystr = header + ".append(" + myitem + ")"
                    print(mystr, file=self.myfile)
            elif type(val3) == dict:
                mystr = header + " = dict()"
                print(mystr, file=self.myfile)
                print(header,"=",val3, file=self.myfile)
            else:
                print(header,"UNSUPPORTED TYPE",type(val3), file=self.myfile)
                self.typeerror += 1
                print(header,"=",val3, file=self.myfile)

    def print_buffet_section(self):
        """Print the buffet section lines. Looping is accomplished
            through an "independent_loop_key" option key.
        """
        print("", file=self.myfile)
        print("###############################################", 
                    file=self.myfile)
        print("#Buffet Command Section", file=self.myfile)
        print("################################################",
                    file=self.myfile)
        print("", file=self.myfile)
        print("from MAST.recipe.recipeinput import RecipeInput",
                    file=self.myfile)
        print("from MAST.buffet.buffetmanager import Buffet",
                    file=self.myfile)
        recipe_base_name = os.path.basename(self.optdict['recipe']['recipe_file'])
        recipe_base_name = recipe_base_name.split('.')[0]
        import time
        timestamp = time.strftime('%Y%m%dH%M%S')
        
        if "independent_loop_key" in self.optdict.keys():
            #for run_idx in xrange(len(SYSTEM_NAME)):
            #    recipe_input = RecipeInput(recipe_name="independent_recipe%s" % run_idx)

            #   #add input parameters
            #    recipe_input.addMastInfo(PROGRAM, SYSTEM_NAME[run_idx])
            #    recipe_input.addStructureInfoFromCoords(COORD_TYPE, COORDINATES[run_idx], LATTICE)
            #    recipe_input.addGlobalIngredientsInfo(ING_GLOBAL_PARAM)
            #    recipe_input.addIngredientsInfo(OPTIMIZE_ING, OPTIMIZE_INFO)
            #    recipe_input.addRecipeInfo(RECIPE_FILE)

            #    #add recipe input to buffet and cook it
            #    buffet_obj = Buffet(name="independent_recipe_%s" % SYSTEM_NAME[run_idx])
            #    buffet_obj.addRecipe(recipe_input)
            #    buffet_obj.cook()
            pass
        
        else:
            sysname = self.optdict['mast']['system_name']
            myname = recipe_base_name + "_" + sysname + "_" + timestamp
            print("recname = '" + myname + "'", file=self.myfile)
            print("recipe_input = RecipeInput(recipe_name=recname)",
                    file=self.myfile)
            print("recipe_input.addMastInfo(program, systemname)",
                    file=self.myfile)
            print("recipe_input.addStructureInfoFromCoords(coord_type, coordinates, lattice)", file=self.myfile) 
            #need to append element info
            print("recipe_input.addGlobalIngredientsInfo(global)",
                    file = self.myfile)
            print("recipe_input.addIngredientsInfo(optimize)",
                    file = self.myfile) 
            #need to be able to add all ingredients
            print("recipe_input.addRecipeInfo(recipe_file)",
                    file = self.myfile)
            print("buffet_obj = Buffet(name='buffet_'+recname)",
                    file=self.myfile)
            print("buffet_obj.addRecipe(recipe_input)", file=self.myfile)
            print("buffet_obj.cook()", file=self.myfile)

