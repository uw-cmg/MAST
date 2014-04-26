##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Glen Jenness
# Last updated: 2014-07-01
##############################################################
import os
import numpy as np
from MAST.utility import MASTError

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
        respective section.
            NO CHANGE if the key was already in the section.
        """
        section_dict = self.options.setdefault(section, dict())
        if 'path' in key:
	        section_dict.setdefault(key, os.path.expanduser(value))
        else:
	        section_dict.setdefault(key, value)

    def get_item(self, section, key=None, default_value=None):
        """Returns the value of the key under this section
        """
        if key is None:
            return self.options.get(section, dict())
        else:
            return self.options.get(section, dict()).get(key, None)

    def update_item(self, section, key, value):
        """Updates the key/value dict item in the option."""
        section_dict = self.options.setdefault(section, dict())
        if key in section_dict.keys():
            print "Changing value of %s, %s." % (section, key)
        if 'path' in key:
            section_dict.__setitem__(key, os.path.expanduser(value))
        else:
            section_dict.__setitem__(key, value)


    def get_sections(self):
        """Return the available sections in the input options.
        """
        return self.options.keys()

    def get_section_keys(self, section):
        """Return the available keys under a section."""
        if not section in self.get_sections():
            raise MASTError(self.__class__.__name__,
                section + " was not parsed from the input file.")
        return self.options[section].keys()

    def reset(self):
        """Option to reset the dict values
        """
        self.options = dict()

    def __repr__(self):
        rlines=""
        rlines=rlines + "Input options: \n"
        sections = self.get_sections()
        sections.sort()
        for section in sections:
            rlines = rlines + "*********************\n"
            rlines = rlines + "*   %s section\n" % section
            rlines = rlines + "*********************\n"
            skeys = self.get_section_keys(section)
            skeys.sort()
            for skey in skeys:
                rlines = rlines + "------------------\n" 
                rlines = rlines + "%s:\n" % skey
                rlines = rlines + "------------------\n" 
                rlines = rlines + "    %s\n" % str(self.get_item(section,skey))
        return rlines
    


    def print_python_section(self, secname):
        """Print python commands that will recreate the InputOptions section
            when the commands are run in python.

            secname <str>: section name
        """
        pln=list()
        if not type(secname) == str:
            errorstr = "Section name " + str(secname) + "must be a string."
            raise MASTError(self.__class__.__name__, errorstr)
        if not secname in self.options.keys():
            errorstr = "No section " + str(secname) + "in self.options."
            raise MASTError(self.__class__.__name__, errorstr)
        initheader=secname
        pln.append(initheader + " = dict()")
        for key1, val1 in self.options[secname].iteritems():
            if key1 == 'structure':
                pln.append("#ignoring created structure")
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
                pln.append(mystr)
            elif type(val1) in [int, bool, float]:
                mystr = header + " = " + str(val1)
                pln.append(mystr)
            elif val1 == None:
                mystr = header + " = None"
                pln.append(mystr)
            elif type(val1) == np.ndarray:
                pln.append(header + "=" + "np.array("+str(val1.tolist())+")")
            elif type(val1) == list:
                mystr = header + " = list()"
                pln.append(mystr)
                for myitem in val1:
                    if type(myitem) == str:
                        mystr = header + ".append(" + "'" + myitem + "'" + ")"
                        pln.append(mystr)
                    else:
                        pln.append(header + ".append(" + str(myitem) + ")")
            elif type(val1) == dict:
                mystr = header + " = dict()"
                pln.append(mystr)
                pln.extend(self.print_secondlevel_dict(header, val1))
            else:
                errstr="UNSUPPORTED TYPE" + str(type(val1))
                raise MASTError(self.__class__.__name__, errstr)
        return pln

    def print_secondlevel_dict(self, initheader1, seconddict):
        """Print a second-level dictionary. 
            
            initheader1 <str>: prepending header
            seconddict <dict>: dictionary that was key1's value
        """
        pln=list()
        for key2, val2 in seconddict.iteritems():
            if type(key2) == str:
                header = initheader1 + "[" + "'" + key2 + "'" + "]"
            else:
                header = initheader1 + "[" + key2 + "]"
            if type(val2) == str:
                mystr = header + " = " + "'" + val2 + "'"
                pln.append(mystr)
            elif type(val2) in [int,float,bool]:
                pln.append(header + "=" + str(val2))
            elif val2 == None:
                mystr = header + " = None"
                pln.append(mystr)
            elif type(val2) == np.ndarray:
                pln.append(header + "=" + "np.array(" + str(val2.tolist()) +")")
            elif type(val2) == list:
                mystr = header + " = list()"
                pln.append(mystr)
                for myitem in val2:
                    if type(myitem) == str:
                        mystr = header + ".append(" + "'" + myitem + "'" + ")"
                        pln.append(mystr)
                    else:
                        pln.append(header + ".append(" + str(myitem) +")")
            elif type(val2) == dict:
                mystr = header + " = dict()"
                pln.append(mystr)
                pln.extend(self.print_thirdlevel_dict(header, val2))
            else:
                errstr="UNSUPPORTED TYPE" + str(type(val2))
                raise MASTError(self.__class__.__name__, errstr)
        return pln

    def print_thirdlevel_dict(self, initheader2, thirddict):
        """Print a third-level dictionary. 
            
            initheader2 <str>: prepending header
            thirddict <dict>: dictionary that was key2's value
        """
        pln=list()
        for key3, val3 in thirddict.iteritems():
            if type(key3) == str:
                header = initheader2 + "[" + "'" + key3 + "'" + "]"
            else:
                header = initheader2 + "[" + key3 + "]"
            if type(val3) == str:
                mystr = header + " = " + "'" + val3 + "'"
                pln.append(mystr)
            elif type(val3) in [int, float, bool]:
                pln.append(header + " = " + str(val3))
            elif val3 == None:
                mystr = header + " = None"
                pln.append(mystr)
            elif type(val3) == np.ndarray:
                pln.append(header + "=" + "np.array(" + str(val3.tolist()) +")")
            elif type(val3) == list:
                mystr = header + " = list()"
                pln.append(mystr)
                for myitem in val3:
                    if type(myitem) == str:
                        mystr = header + ".append(" + "'" + myitem + "'" + ")"
                        pln.append(mystr)
                    else:
                        pln.append(header + ".append(" + str(myitem) +")")
            elif type(val3) == dict:
                mystr = header + " = dict()"
                pln.append(mystr)
                pln.extend(self.print_fourthlevel_dict(header, val3))
            else:
                errstr="UNSUPPORTED TYPE" + str(type(val3))
                raise MASTError(self.__class__.__name__, errstr)
        return pln
    
    def print_fourthlevel_dict(self, initheader2, thirddict):
        """Print a third-level dictionary. 
            
            initheader2 <str>: prepending header
            thirddict <dict>: dictionary that was key2's value
        """
        pln=list()
        for key3, val3 in thirddict.iteritems():
            if type(key3) == str:
                header = initheader2 + "[" + "'" + key3 + "'" + "]"
            else:
                header = initheader2 + "[" + key3 + "]"
            if type(val3) == str:
                mystr = header + " = " + "'" + val3 + "'"
                pln.append(mystr)
            elif type(val3) in [int, float, bool]:
                pln.append(header + " = " + str(val3))
            elif val3 == None:
                mystr = header + " = None"
                pln.append(mystr)
            elif type(val3) == np.ndarray:
                pln.append(header + "=" + "np.array(" + str(val3.tolist()) +")")
            elif type(val3) == list:
                mystr = header + " = list()"
                pln.append(mystr)
                for myitem in val3:
                    if type(myitem) == str:
                        mystr = header + ".append(" + "'" + myitem + "'" + ")"
                        pln.append(mystr)
                    else:
                        pln.append(header + ".append(" + str(myitem) +")")
            elif type(val3) == dict:
                mystr = header + " = dict()"
                pln.append(mystr)
                pln.append(header + "=" + str(val3))
            else:
                errstr="UNSUPPORTED TYPE" + str(type(val3))
                raise MASTError(self.__class__.__name__, errstr)
        return pln
