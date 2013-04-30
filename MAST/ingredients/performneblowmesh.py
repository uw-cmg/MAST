from MAST.ingredients.performneb import PerformNEB
from MAST.ingredients.baseingredient import BaseIngredient
from pymatgen.core.structure import Structure

class PerformNEBLowMesh(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of optimization directory'),
            'program': (str, str(), 'DFT program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'child_dict': (dict, dict(), 'Dictionary of children'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)


    def is_complete(self):
        return PerformNEB.is_complete(self)

    def update_children(self):
        return PerformNEB.update_children(self)

    def write_files(self):
        return PerformNEB.write_files(self)

    def get_my_numbers(self):
        return PerformNEB.get_my_numbers(self)

    def get_parent_structures(self):
        return PerformNEB.get_parent_structures(self)

    def do_interpolation(self, parentstructures):
        return PerformNEB.do_interpolation(self, parentstructures)

    def place_parent_energy_files(self):
        return PerformNEB.place_parent_energy_files(self)

