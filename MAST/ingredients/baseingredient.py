from MAST.utility import MASTObj

class BaseIngredient(MASTObj):
    def __init__(self, allowed_keys, **kwargs):
        allowed_keys_base = dict()
        allowed_keys_base.update(allowed_keys) 
        MASTObj.__init__(self, allowed_keys_base, **kwargs)
        #self.logger    = logger #keep this space
        #self.structure = dict() #TTM 2013-03-27 structure is in allowed_keys
    
    def write_directory(self):
        if os.path.exists(self.keywords['name']):
            print "Directory exists."
            return
        os.makedirs(self.keywords['name'])
        return
    
    def get_structure_from_parent(self, parentpath):
        if self.keywords['program'] == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.get_structure_from_parent(parentpath)
        else:
            print "Program not recognized (in get_structure_from_parent)"
            return None

    def forward_parent_structure(self, parentpath, childpath):
        if self.keywords['program'] == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            vasp_checker.forward_parent_structure(parentpath, childpath)
            return None
        else:
            print "Program not recognized (in forward_parent_structure)"
            return None

    def write_files(self):
        '''writes the files needed as input for the jobs'''

    def is_complete(self):
        '''Function to check if Ingredient is ready'''

    def images_complete(self):
        '''Checks if all images in an NEB calculation are complete.'''
        if self.keywords['program'] == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.images_complete(self.keywords['name'],self.keywords['program_keys']['images'])
        else:
            print "Program not recognized (in images_complete)"
            return None
