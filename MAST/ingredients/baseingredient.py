import os
import subprocess
import time

from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import Metadata
from MAST.ingredients.checker import BaseChecker
from MAST.ingredients.checker import VaspChecker
from MAST.ingredients.checker import PhonChecker

class BaseIngredient(MASTObj):
    """Base Ingredient class
        Attributes:
            self.meta_dict <dict>: Metadata dictionary
            self.metafile <Metadata>: Metadata file
            self.program <str>: program name, all lowercase,
                                from 'mast_program' in input
                                file
            self.checker <BaseChecker, PhonChecker, VaspChecker>:
                    program-dependent checker object
    """
    def __init__(self, allowed_keys, **kwargs):
        allowed_keys_base = dict()
        allowed_keys_base.update(allowed_keys) 
        MASTObj.__init__(self, allowed_keys_base, **kwargs)

        work_dir = '/'.join(self.keywords['name'].split('/')[:-1])
        topmeta = Metadata(metafile='%s/metadata.txt' % work_dir)
        data = topmeta.read_data(self.keywords['name'].split('/')[-1])
        #print 'GRJ DEBUG: name =', self.keywords['name'].split('/')[-1]
        #print 'GRJ DEBUG: data =', data
        #print 'GRJ DEBUG: topmeta =\n', topmeta

        self.meta_dict = dict()
        if data:
            for datum in data.split(','):
                self.meta_dict[datum.split(':')[0]] = datum.split(':')[1].strip()

        self.metafile = Metadata(metafile='%s/metadata.txt' % self.keywords['name'])

        self.program = self.keywords['program_keys']['mast_program'].lower()
        if self.program == 'vasp':
            self.checker = VaspChecker(name=self.keywords['name'],
            program_keys = self.keywords['program_keys'],
            structure = self.keywords['structure'])
        elif self.program == 'phon':
            self.checker = PhonChecker(name=self.keywords['name'],program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])
        else:
            allowed_keys={'name','program_keys','structure'}
            self.checker = BaseChecker(allowed_keys, name=self.keywords['name'],program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])

        #self.logger    = logger #keep this space
        #self.structure = dict() #TTM 2013-03-27 structure is in allowed_keys

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
            print "Directory exists."
            return
        return

    def is_complete(self):
        '''Function to check if Ingredient is ready'''
        if self.program == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            usepath = self.keywords['name']
            if vasp_checker._vasp_is_neb(self.keywords):
                mycomplete = vasp_checker.images_complete(self.keywords['name'],
                                self.keywords['program_keys']['images'])
                usepath = usepath + '/01'
            else:
                mycomplete = vasp_checker.is_complete(usepath)

            if mycomplete:
                self.metafile.write_data('completed on', time.asctime())
                if 'OSZICAR' in os.listdir(self.keywords['name']):
                    energy = self.checker.get_energy_from_energy_file(self.keywords['name'])
                else:
                    energy = None
                self.metafile.write_data('energy', energy)

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

        elif self.program == 'phon':
            from MAST.ingredients.checker import phon_checker
            usepath = self.keywords['name']
            mycomplete = phon_checker.is_complete(usepath)
            self.metafile.write_data('completed on', time.asctime())
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
        if self.program == 'vasp':
            from MAST.ingredients.checker import vasp_checker
            return vasp_checker.is_ready_to_run(self.keywords['name'])
        elif self.program == 'phon':
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
            queuesub = queue_commands.write_to_submit_list(self.keywords['name'])
            #runme = subprocess.Popen(queuesub, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #runme.wait()
            # for scheduling other jobs
            #runme.wait()
        os.chdir(curdir)
        self.metafile.write_data('run start', time.asctime())

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
