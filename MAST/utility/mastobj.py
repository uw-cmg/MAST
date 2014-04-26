##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Glen Jenness
# Last updated: 2013-07-01
##############################################################

#import numpy as np

#import os
#import warnings
from MAST.utility.masterror import MASTError


class MASTObj(object):
    """Base MAST class.  Used to handle keyword input into MAST classes in a 
    uniform fashion.

    For examples see:
        MAST/interface/interface
        MAST/defects/defecttype

    Attributes:
        allowed_keys <dict>: Each key in this dictionary is a keyword for that
                             class, with each keyvalue being a tuple containing
                             the keyword type, default value, and doc string
        keywords <dict>: The **kwargs passed into the class
    """
    def __init__(self, allowed_keys, **kwargs):
        """
        Initializer for MASTobj.
        Calls:
           self.set_keywords

        Sets:
           self.keywords
           self.allowed_keys

        Returns:
           Nothing 
        """
        self.allowed_keys = allowed_keys
        self.keywords = self.set_keywords(**kwargs)

    def set_keywords(self, **kwargs):
        """
        Takes the self.allowed_keys dictionary, and sets appropiate defaults
        for keywords given in **kwargs and checks to make sure all keywords in
        **kwargs are valid

        Calls:
            Nothing

        Sets:
            Nothing

        Returns:
            Keywords
        """
        keywords = dict()

        for key, value in self.allowed_keys.items():
            keywords[key] = value[1]

        for key, value in kwargs.items():
            if key not in self.allowed_keys:
                error = 'Keyword %s for %s object not found' % \
                                   (key, self.__class__.__name__)
                MASTError(self.__class__.__name__, error)

#                raise RuntimeError('Keyword %s for %s object not found' % \
#                                   (key, self.__class__.__name__))

            if isinstance(value, self.allowed_keys[key][0]):
                keywords[key] = value
            else:
                error = 'Keyword %s value %s invalid; expected type %s, got type %s' % (key, str(value), self.allowed_keys[key][0], type(value))
                MASTError(self.__class__.__name__, error)
#                raise RuntimeError('Keyword %s value invalid' % key)

        return keywords

    def help(self, keyword):
        """Prints a help message listing all the allowed keywords and any relevant information"""
        if (keyword == 'all'):
            string = ('%-20s%-20s%-20s%s\n' % ('Keyword', 'Type', 'Default', 'Comment'))
            for key, value in self.allowed_keys.items():
                string += ('%-20s%-20s%-20s%s\n' % (key, str(value[0]), str(value[1]), value[2]))
            print string

