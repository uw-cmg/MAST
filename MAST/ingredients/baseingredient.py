import os
import subprocess
import time
import logging
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import Metadata
from MAST.ingredients.checker import BaseChecker
from MAST.ingredients.checker import VaspChecker
from MAST.ingredients.checker import VaspNEBChecker
from MAST.ingredients.checker import PhonChecker
from MAST.ingredients.errorhandler import BaseError
from MAST.ingredients.errorhandler import VaspError
from MAST.ingredients.errorhandler import PhonError
from MAST.ingredients.errorhandler import VaspNEBError

class BaseIngredient(MASTObj):
    """Base Ingredient class
        Attributes:
            self.meta_dict <dict>: Metadata dictionary
            self.metafile <Metadata>: Metadata file
            self.program <str>: program name, all lowercase,
                                from 'mast_program' in input
                                file
            self.checker <VaspChecker, PhonChecker, etc.>:
                    program-dependent checker object
            self.errhandler <VaspError, PhonError, etc.>:
                    program-dependent handler object
    """
    def __init__(self, allowed_keys, **kwargs):
        allowed_keys_base = dict()
        allowed_keys_base.update(allowed_keys) 
        MASTObj.__init__(self, allowed_keys_base, **kwargs)

        work_dir = '/'.join(self.keywords['name'].split('/')[:-1])
        topmeta = Metadata(metafile='%s/metadata.txt' % work_dir)
        data = topmeta.read_data(self.keywords['name'].split('/')[-1])

        self.meta_dict = dict()
        if data:
            for datum in data.split(','):
                self.meta_dict[datum.split(':')[0]] = datum.split(':')[1].strip()

        self.metafile = Metadata(metafile='%s/metadata.txt' % self.keywords['name'])

        self.program = self.keywords['program_keys']['mast_program'].lower()
        
        logging.basicConfig(filename="%s/mast.log" % os.getenv("MAST_CONTROL"), level=logging.DEBUG)
        self.logger = logging.getLogger(__name__)
        
        if self.program == 'vasp':
            self.checker = VaspChecker(name=self.keywords['name'],
            program_keys = self.keywords['program_keys'],
            structure = self.keywords['structure'])
            self.errhandler = VaspError(name=self.keywords['name'],
            program_keys = self.keywords['program_keys'],
            structure = self.keywords['structure'])
        elif self.program == 'vasp_neb':
            self.checker = VaspNEBChecker(name=self.keywords['name'],
            program_keys = self.keywords['program_keys'],
            structure = self.keywords['structure'])
            self.errhandler = VaspNEBError(name=self.keywords['name'],
            program_keys = self.keywords['program_keys'],
            structure = self.keywords['structure'])
        elif self.program == 'phon':
            self.checker = PhonChecker(name=self.keywords['name'],program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])
            self.errhandler = PhonError(name=self.keywords['name'],program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])
        else:
            allowed_keys={'name','program_keys','structure'}
            self.checker = BaseChecker(allowed_keys, name=self.keywords['name'],program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])
            self.errhandler = BaseError(allowed_keys, name=self.keywords['name'],program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])

    def write_directory(self):
        try:
            os.makedirs(self.keywords['name'])
            self.metafile.write_data('directory created', time.asctime())
            self.metafile.write_data('name', self.keywords['name'].split('/')[-1])
            self.metafile.write_data('program', self.program)
            self.metafile.write_data('ingredient type', self.__class__.__name__)
            #if 'mast_charge' in self.keywords['program_keys']:
            #    self.metafile.write_data('charge', self.keywords['program_keys']['mast_charge'])
            for key, value in self.meta_dict.items():
                self.metafile.write_data(key, value)
        except OSError:
            self.logger.info("Directory for %s already exists." % self.keywords['name'])
            return
        return

    def is_complete(self):
        '''Function to check if Ingredient is ready'''
        if self.checker.is_complete():
            self.metafile.write_data('completed on', time.asctime())
            if 'get_energy_from_energy_file' in dirutil.list_methods(self.checker):
                energy = self.checker.get_energy_from_energy_file(self.keywords['name'])
                self.metafile.write_data('energy', energy)
            return True
        else:
            if not self.checker.is_started():
                return False #hasn't started running yet.
            else:
                errct = self.errhandler.loop_through_errors()
                if errct > 0:
                    raise MASTError(self.__class__.__name__,"Error found for %s. Please correct this error or move this recipe out of $MAST_SCRATCH. Delete the $MAST_SCRATCH/mast.write_files.lock file if necessary." % self.keywords['name'])
                    #if not self.checker.is_on_queue():
                    #    self.run()
                return False

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
        return self.checker.is_ready_to_run()
        
    def getpath(self):
        '''getpath returns the directory of the ingredient'''
        return self.keywords['name']
    
    def write_submit_script(self):
        from submit import script_commands
        script_commands.write_submit_script(self.keywords)
        return
    
    def run(self, mode='serial', curdir=os.getcwd()):
        from submit import queue_commands 
        

        if mode.lower() == 'noqsub':
            curdir = os.getcwd()
            os.chdir(self.keywords['name'])
            programpath = queue_commands.direct_shell_command()
            p = subprocess.Popen(programpath, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait()
            os.chdir(curdir)
            self.metafile.write_data('run', time.asctime())
            
        elif mode.lower() == 'serial':
            queuesub = queue_commands.write_to_submit_list(self.keywords['name'])
            #runme = subprocess.Popen(queuesub, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #runme.wait()
            # for scheduling other jobs
            #runme.wait()
            self.metafile.write_data('queued', time.asctime())
        return

    def get_name(self):
        return self.keywords['name'].split('/')[-1]

    @property
    def name(self):
        return self.get_name()

    def get_keywords(self):
        return self.keywords.copy()

    def __repr__(self):
        return 'Ingredient %s of type %s' % (self.keywords['name'].split('/')[-1], self.__class__.__name__)


    
    def get_my_label(self, label):
        """Get the value of a label in the metadata file.
            Args:
                label <str>: Label to search for.
            Returns:
                <str>: Value of the label as a string, stripped.
        """
        myname = self.keywords['name']
        mymeta = Metadata(metafile=os.path.join(myname, "metadata.txt"))
        mylabel = mymeta.search_data(label)
        if mylabel == "":
            raise MASTError(self.__class__.__name__, 
                "No metadata for tag %s" % label)
        return mylabel[1]
