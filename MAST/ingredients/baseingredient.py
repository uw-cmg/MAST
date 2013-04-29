from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
import os
import subprocess

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
            return vasp_checker.images_complete(self.keywords['name'],\
                   self.keywords['program_keys']['images'])
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
        from submit import script_commands
        script_commands.write_submit_script(self.keywords)
        return
    
    def run(self, mode='noqsub', curdir=os.getcwd()):
        from submit import queue_commands 
        
        curdir = os.getcwd()
        os.chdir(self.keywords['name'])

        if mode == 'noqsub':
            programpath = queue_commands.direct_shell_command()
            p = subprocess.call([programpath])
            p.wait()
            
        elif mode == 'serial':
            queuesub = queue_commands.queue_submission_command()
            runme = subprocess.Popen(queuesub, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            runme.wait()
            # for scheduling other jobs
            #runme.wait()
        os.chdir(curdir)
        return

    def set_up_program_input(self):
        if self.keywords['program'] == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.set_up_program_input(self.keywords)
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in set_up_program_input)")

    def get_path_to_write_neb_parent_energy(self, parent):
        """Get path to write the NEB's parent energy file.
            parent = 1 for initial, 2 for final
        """
        if self.keywords['program'] == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.get_path_to_write_neb_parent_energy(self.keywords['name'], self.keywords['program_keys']['images'],parent)
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in get_path_to_write_neb_parent_energy)")

    def set_up_program_input_neb(self, image_structures):
        if self.keywords['program'] == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.set_up_program_input_neb(self.keywords, image_structures)
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in set_up_neb)")


# The following functions need to be defined by the child class:
    def write_files(self):
        '''writes the files needed as input for the jobs'''
        raise NotImplementedError

