import os
import subprocess
import time

from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import Metadata

class BaseIngredient(MASTObj):
    def __init__(self, allowed_keys, **kwargs):
        allowed_keys_base = dict()
        allowed_keys_base.update(allowed_keys) 
        MASTObj.__init__(self, allowed_keys_base, **kwargs)

        self.metafile = Metadata(metafile='%s/metadata.txt' % self.keywords['name'])
        #self.logger    = logger #keep this space
        #self.structure = dict() #TTM 2013-03-27 structure is in allowed_keys

    def write_directory(self):
        try:
            os.makedirs(self.keywords['name'])
            self.metafile.write_data('Directory created', time.asctime())
            self.metafile.write_data('Name', self.keywords['name'].split('/')[-1])
            self.metafile.write_data('Program', self.keywords['program'])
            self.metafile.write_data('Ingredient type', self.__class__.__name__)
            if 'mast_charge' in self.keywords['program_keys']:
                self.metafile.write_data('Charge', self.keywords['program_keys']['mast_charge'])
        except OSError:
            print "Directory exists."
            return
        return
  
    def get_structure_from_directory(self, dirname):
        if self.keywords['program'].lower() == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.get_structure_from_directory(dirname)
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in get_structure_from_directory)")

    def get_structure_from_file(self, filepath):
        if self.keywords['program'].lower() == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.get_structure_from_file(filepath)
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in get_structure_from_file)")

    def forward_parent_structure(self, parentpath, childpath, newname="POSCAR"):
        if self.keywords['program'].lower() == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            vasp_checker.forward_parent_structure(parentpath, childpath, newname)
            return None
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in forward_parent_structure)")
    
    def forward_parent_energy(self, parentpath, childpath, newname="OSZICAR"):
        if self.keywords['program'].lower() == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            vasp_checker.forward_parent_energy(parentpath, childpath, newname)
            return None
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in forward_parent_structure)")

    def forward_parent_dynmat(self, parentpath, childpath, newname="DYNMAT"):
        if self.keywords['program'].lower() == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            vasp_checker.forward_parent_dynmat(parentpath, childpath, newname)
            return None
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in forward_parent_structure)")
    def is_complete(self):
        '''Function to check if Ingredient is ready'''
        if self.keywords['program'].lower() == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            usepath = self.keywords['name']
            if 'images' in self.keywords['program_keys'].keys():
                mycomplete = vasp_checker.images_complete(self.keywords['name'],
                                self.keywords['program_keys']['images'])
                usepath = usepath + '/01'
            else:
                mycomplete = vasp_checker.is_complete(usepath)

            if mycomplete:
                self.metafile.write_data('Completed on', time.asctime())
                if 'vasprun.xml' in os.listdir(self.keywords['name']):
                    energy = vasp_checker.get_vasp_energy(self.keywords['name'])
                else:
                    energy = None
                self.metafile.write_data('Energy', energy)

                return mycomplete
            else:
                if not os.path.exists(usepath + '/OUTCAR'):
                    return False #hasn't started running yet.
                from MAST.ingredients.errorhandler import vasp_error
                if 'images' in self.keywords['program_keys'].keys():
                    errct = vasp_error.loop_through_errors(usepath, 1)
                else:
                    errct = vasp_error.loop_through_errors(usepath)
                if errct > 0:
                    pass #self.run() #Should try to rerun automatically or not?? NO.
                return False

        elif self.keywords['program'].lower() == 'phon':
            from MAST.ingredients.checker import phon_checker
            usepath = self.keywords['name']
            mycomplete = phon_checker.is_complete(usepath)
            self.metafile.write_data('Completed on', time.asctime())
            return mycomplete
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in is_complete)")

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
        if self.keywords['program'].lower() == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.is_ready_to_run(self.keywords['name'])
        elif self.keywords['program'].lower() == 'phon':
            from MAST.ingredients.checker import phon_checker
            return phon_checker.is_ready_to_run(self.keywords['name'])
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
    
    def run(self, mode='serial', curdir=os.getcwd()):
        from submit import queue_commands 
        
        curdir = os.getcwd()
        os.chdir(self.keywords['name'])

        if mode.lower() == 'noqsub':
            programpath = queue_commands.direct_shell_command()
            p = subprocess.Popen(programpath, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait()
            
        elif mode.lower() == 'serial':
            queuesub = queue_commands.queue_submission_command()
            runme = subprocess.Popen(queuesub, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            runme.wait()
            # for scheduling other jobs
            #runme.wait()
        os.chdir(curdir)
        self.metafile.write_data('Run start', time.asctime())

        return

    def set_up_program_input(self):
        if self.keywords['program'].lower() == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.set_up_program_input(self.keywords)
        elif self.keywords['program'].lower() == 'phon':
            from MAST.ingredients.checker import phon_checker
            return phon_checker.set_up_program_input(self.keywords)
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in set_up_program_input)")

    def get_path_to_write_neb_parent_energy(self, parent):
        """Get path to write the NEB's parent energy file.
            parent = 1 for initial, 2 for final
        """
        if self.keywords['program'].lower() == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.get_path_to_write_neb_parent_energy(self.keywords['name'], self.keywords['program_keys']['images'],parent)
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in get_path_to_write_neb_parent_energy)")

    def set_up_program_input_neb(self, image_structures):
        if self.keywords['program'].lower() == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.set_up_program_input_neb(self.keywords, image_structures)
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in set_up_neb)")

    def get_children(self):
        """Returns the children of this ingredient.
            If there are no children, it will return None instead.
        """
        if (self.keywords['child_dict']):
            return self.keywords['child_dict'].copy()
        else:
            return None

    @property
    def children(self):
        return self.get_children()

    def get_name(self):
        return self.keywords['name'].split('/')[-1]

    @property
    def name(self):
        return self.get_name()

    def get_keywords(self):
        return self.keywords.copy()

    def __repr__(self):
        return 'Ingredient %s of type %s' % (self.keywords['name'].split('/')[-1], self.__class__.__name__)

    def add_selective_dynamics_to_structure(self, sdarray):
        """Adds selective dynamics to a structure."""
        if self.keywords['program'].lower() == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.add_selective_dynamics_to_structure(self.keywords, sdarray)
        else:
            raise MASTError(self.__class__.__name__, 
                "Program not recognized (in add_selective_dynamics_to_structure)")

# The following functions need to be defined by the child class:
    def write_files(self):
        '''writes the files needed as input for the jobs'''
        raise NotImplementedError
    


