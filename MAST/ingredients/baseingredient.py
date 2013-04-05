from MAST.utility import MASTObj
from MAST.utility import MASTError

class BaseIngredient(MASTObj):
    def __init__(self, allowed_keys, **kwargs):
        allowed_keys_base = dict()
        allowed_keys_base.update(allowed_keys) 
        MASTObj.__init__(self, allowed_keys_base, **kwargs)
        #self.logger    = logger #keep this space
        #self.structure = dict() #TTM 2013-03-27 structure is in allowed_keys
    
    def write_directory(self):
        try:
            os.makedirs(self.keywords['name'])
        except OSError:
            print "Directory exists."
            return
        return
   
    def get_structure_from_file(self, filepath):
        if self.keywords['program'] == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.get_structure_from_file(filepath)
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in get_structure_from_file)")

    def forward_parent_structure(self, parentpath, childpath, newname="POSCAR"):
        if self.keywords['program'] == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            vasp_checker.forward_parent_structure(parentpath, childpath, newname)
            return None
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in forward_parent_structure)")


    def write_files(self):
        '''writes the files needed as input for the jobs'''

    def is_complete(self):
        '''Function to check if Ingredient is ready'''
        if self.keywords['program'] == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.is_complete(self.keywords['name'])
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in is_complete)")

    def images_complete(self):
        '''Checks if all images in an NEB calculation are complete.'''
        if self.keywords['program'] == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.images_complete(self.keywords['name'],self.keywords['program_keys']['images'])
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in images_complete)")
