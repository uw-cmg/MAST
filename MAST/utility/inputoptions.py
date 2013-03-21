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

class InputOptions:
    """Stores options information in dictionary format for retrieval

    Dictionary is organised as sections and key, value pairs within a section

    Attributes:
        options <dict of sections 
                    which is a dict of key/values>: stores options which comprises
                                                    of sections and key/values under it
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

    def get_item(self, section, key):
        """Returns the value of the key under this section
        """
        return self.options.get(section, dict()).get(key)

    def reset(self):
        """Option to reset the dict values
        """
        self.options = dict()

    def __repr__(self):
        return str(self.options)
