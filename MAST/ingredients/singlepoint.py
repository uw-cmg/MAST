from MAST.ingredients.optimize import Optimize
from MAST.ingredients.baseingredient import BaseIngredient
from pymatgen.core.structure import Structure

class SinglePoint(Optimize):
    def __init__(self, **kwargs):
        Optimize.__init__(self, **kwargs)
