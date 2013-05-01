from MAST.ingredients.performneb import PerformNEB
from MAST.ingredients.baseingredient import BaseIngredient
from pymatgen.core.structure import Structure

class PerformNEBLowMesh(PerformNEB):
    def __init__(self, **kwargs):
        PerformNEB.__init__(self, **kwargs)
