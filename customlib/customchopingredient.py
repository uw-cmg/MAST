import os
import shutil
import sys
import time
from MAST.ingredients.baseingredient import BaseIngredient
from pymatgen.core.structure import Structure
class CustomChopIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        """Please not modify this init method."""
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program': (str, str(), 'Program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)


    def my_custom_chop_method(self, optional_string=""):
        self.logger.info("Here is a custom function within a custom class.")
        self.logger.info("This function can access the ingredient's program keys.")
        self.logger.info("%s" % self.keywords['program_keys'])
        return "Hello! %s" % optional_string

