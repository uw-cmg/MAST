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
import glob
import inspect

from MAST import ingredients
from MAST.ingredients import baseingredient
from MAST.ingredients import optimize
from MAST.ingredients import performneb

class IngredientsLoader:
    def __init__(self):
        '''Loads the ingredients library dict'''
        self.ingredients_dict = {}

    def load_ingredients(self):
        '''Loads the ingredients dynamically from ingredients folder
           and returns the ingredients dict
        '''
        ingredients_dir = ingredients.__path__[0]
        cand_ingredients = map(lambda x: os.path.basename(x).replace(".py", ""),
                               glob.glob("%s/*.py" % ingredients_dir))

        for ingredient in cand_ingredients:
            if ingredient.startswith('_'):
                continue

            mod = __import__("MAST.ingredients.%s" % ingredient, globals(), locals(), [ingredient])

            for attr_name, attr in inspect.getmembers(mod):
                if attr_name.startswith('_'):
                    continue
                if inspect.isclass(attr):
                    for base_class in attr.__bases__:
                        if (base_class == baseingredient.BaseIngredient) or (base_class == optimize.Optimize) or (base_class == neb.NEB):
                        #if (base_class == ingredients.baseingredient.BaseIngredient) or (base_class == ingredients.optimize.Optimize) or (base_class == ingredients.performneb.PerformNEB):
                            self.ingredients_dict[attr_name.lower()] = attr

