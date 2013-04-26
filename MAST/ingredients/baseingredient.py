from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
import os

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
    

    def directory_is_locked(self):
        return dirutil.directory_is_locked(self.keywords['name'])

    def lock_directory(self):
        return dirutil.lock_directory(self.keywords['name'])

    def unlock_directory(self):
        return dirutil.unlock_directory(self.keywords['name'])

    def wait_to_write(self):
        return dirutil.wait_to_write(self.keywords['name'])
   
    def is_ready_to_run(self):
        if self.directory_is_locked():
            return False
        if self.keywords['program'] == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.is_ready_to_run(self.keywords['name'])
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in is_complete)")
        
    def getpath(self):
        '''getpath returns the directory of the ingredient'''
        return self.keywords['name']
    
    def write_submit_script(self):
        if not ('script' in self.keywords['program_keys'].keys()):
            templatename = 'submitscript.sh'
        else:
            templatename = self.keywords['program_keys']['script']
        templatepath = os.path.join(dirutil.get_mast_install_path(),
                        'submit',templatename)
        bname = os.path.basename(self.keywords['name'])
        wpath = self.keywords['name'] + '/submit.sh'
        #print wpath
        #print bname
        from submit import script_commands
        script_commands.modify_jobname(templatepath, wpath, bname)
        return
