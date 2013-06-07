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

class InputOptions:
    """Stores options information in dictionary format for retrieval
        Dictionary is organised as sections and key, value pairs 
            within a section
        Attributes:
            options <dict of sections which is a dict of key/values>: 
                stores options which consist of sections and key/values 
                under the sections
        e.g.: section mast, keyword 'program'
                inputoptions.options['mast']['program'] = 'vasp'
    """

    def __init__(self):
        """Constructor for InputOptions
        Initialize the dictionary
        """
        self.options = dict()

    def set_item(self, section, key, value):
        """Adds the key/value dict item into the option under the 
        respective section
        """
        section_dict = self.options.setdefault(section, dict())
        if 'path' in key:
	        section_dict.setdefault(key, os.path.expanduser(value))
        else:
	        section_dict.setdefault(key, value)

    def get_item(self, section, key, default_value=None):
        """Returns the value of the key under this section
        """
        return self.options.get(section, dict()).get(key, None)

    def reset(self):
        """Option to reset the dict values
        """
        self.options = dict()

    def __repr__(self):
        return str(self.options)
    


    def print_python_section(self, prepend, secname, myfile):
        """Print python commands that will recreate the InputOptions section
            when the commands are run in python.

            prepend <str>: prepend string, like input_option_varname.options
            secname <str>: section name
            myfile <open file>: file to which to write
        """
        if not type(secname) == str:
            errorstr = "Section name " + str(secname) + "must be a string."
            raise MASTError(self.__class__.__name__, errorstr)
        if not secname in self.options.keys():
            errorstr = "No section " + str(secname) + "in self.options."
            raise MASTError(self.__class__.__name__, errorstr)
        initheader=prepend + "['" + secname + "']"
        print(initheader + " = dict()", file=myfile)
        for key1, val1 in self.options[secname].iteritems():
            if key1 == 'structure':
                print("#ignoring created structure", file=myfile)
                continue
            if type(key1) == str:
                header = initheader + "[" + "'" + key1 + "'" + "]"
            else:
                header = initheader + "[" + key1 + "]"
            if not type(key1) == str:
                errorstr = "One of the keys in section names is not a string."
                raise MASTError(self.__class__.__name__, errorstr)
            if type(val1) == str:
                mystr = header + " = " + "'" + val1 + "'"
                print(mystr, file=myfile)
            elif type(val1) in [int, bool, float]:
                mystr = header + " = ", val1
                print(header,"=",val1, file=myfile)
            elif val1 == None:
                mystr = header + " = None"
                print(mystr, file=myfile)
            elif type(val1) == np.ndarray:
                print(header,"=","np.array(",val1.tolist(),")",file=myfile)
            elif type(val1) == list:
                mystr = header + " = list()"
                print(mystr, file=myfile)
                for myitem in val1:
                    if type(myitem) == str:
                        mystr = header + ".append(" + "'" + myitem + "'" + ")"
                    else:
                        mystr = header + ".append(" + myitem + ")"
                    print(mystr, file=myfile)
            elif type(val1) == dict:
                mystr = header + " = dict()"
                print(mystr, file=myfile)
                self.print_secondlevel_dict(header, val1, myfile)
            else:
                errstr="UNSUPPORTED TYPE" + str(type(val1))
                raise MASTError(self.__class__.__name__, errstr)

    def print_secondlevel_dict(self, initheader1, seconddict, myfile):
        """Print a second-level dictionary. 
            
            initheader1 <str>: prepending header
            seconddict <dict>: dictionary that was key1's value
            myfile <open file>: file for writing
        """
        for key2, val2 in seconddict.iteritems():
            if type(key2) == str:
                header = initheader1 + "[" + "'" + key2 + "'" + "]"
            else:
                header = initheader1 + "[" + key2 + "]"
            if type(val2) == str:
                mystr = header + " = " + "'" + val2 + "'"
                print(mystr, file=myfile)
            elif type(val2) in [int,float,bool]:
                print(header,"=",val2, file=myfile)
            elif val2 == None:
                mystr = header + " = None"
                print(mystr, file=myfile)
            elif type(val2) == np.ndarray:
                print(header,"=","np.array(",val2.tolist(),")",file=myfile)
            elif type(val2) == list:
                mystr = header + " = list()"
                print(mystr, file=myfile)
                for myitem in val2:
                    if type(myitem) == str:
                        mystr = header + ".append(" + "'" + myitem + "'" + ")"
                    else:
                        mystr = header + ".append(" + myitem + ")"
                    print(mystr, file=myfile)
            elif type(val2) == dict:
                mystr = header + " = dict()"
                print(mystr, file=myfile)
                self.print_thirdlevel_dict(header, val2, myfile)
            else:
                errstr="UNSUPPORTED TYPE" + str(type(val2))
                raise MASTError(self.__class__.__name__, errstr)

    def print_thirdlevel_dict(self, initheader2, thirddict, myfile):
        """Print a third-level dictionary. 
            
            initheader2 <str>: prepending header
            thirddict <dict>: dictionary that was key2's value
            myfile <open file>: file for writing
        """
        for key3, val3 in thirddict.iteritems():
            if type(key3) == str:
                header = initheader2 + "[" + "'" + key3 + "'" + "]"
            else:
                header = initheader2 + "[" + key3 + "]"
            if type(val3) == str:
                mystr = header + " = " + "'" + val3 + "'"
                print(mystr, file=myfile)
            elif type(val3) in [int, float, bool]:
                print(header,"=",val3, file=myfile)
            elif val3 == None:
                mystr = header + " = None"
                print(mystr, file=myfile)
            elif type(val3) == np.ndarray:
                print(header,"=","np.array(",val3.tolist(),")",file=myfile)
            elif type(val3) == list:
                mystr = header + " = list()"
                print(mystr, file=myfile)
                for myitem in val3:
                    if type(myitem) == str:
                        mystr = header + ".append(" + "'" + myitem + "'" + ")"
                    else:
                        mystr = header + ".append(" + myitem + ")"
                    print(mystr, file=myfile)
            elif type(val3) == dict:
                mystr = header + " = dict()"
                print(mystr, file=myfile)
                print(header,"=",val3, file=myfile)
            else:
                errstr="UNSUPPORTED TYPE" + str(type(val2))
                raise MASTError(self.__class__.__name__, errstr)
