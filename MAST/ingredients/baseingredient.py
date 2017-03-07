##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import subprocess
import time
import logging
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import Metadata
from MAST.utility import MASTFile
from MAST.utility import loggerutils
from MAST.ingredients.checker import BaseChecker
from MAST.ingredients.checker import VaspChecker
from MAST.ingredients.checker import VaspNEBChecker
from MAST.ingredients.checker import PhonChecker
from MAST.ingredients.checker import GenericChecker
from MAST.ingredients.checker import LammpsChecker
from MAST.ingredients.errorhandler import BaseError
from MAST.ingredients.errorhandler import VaspError
from MAST.ingredients.errorhandler import PhonError
from MAST.ingredients.errorhandler import VaspNEBError
from MAST.ingredients.errorhandler import GenericError
from MAST.ingredients.pmgextend.atom_index import AtomIndex
from MAST.utility import loggerutils

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
            self.atomindex <AtomIndex>: atom index object
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
            for datum in data.split(';'):
                self.meta_dict[datum.split(':')[0]] = datum.split(':')[1].strip()

        self.metafile = Metadata(metafile='%s/metadata.txt' % self.keywords['name'])

        self.program = self.keywords['program_keys']['mast_program'].lower()
        
        self.logger = loggerutils.get_mast_logger(self.keywords['name'])
        
        sdir=os.path.join(os.path.dirname(self.keywords['name']),"structure_index_files")
        if os.path.exists(sdir):
            self.atomindex = AtomIndex(structure_index_directory=sdir)
        else:
            self.atomindex = None

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
        elif self.program == 'lammps':
            self.checker = LammpsChecker(name=self.keywords['name'],program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])
            self.errhandler = GenericError(name=self.keywords['name'],program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])
        elif self.program =='structopt':
            self.checker = StructoptChecker(name=self.keywords['name'],program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])
            self.errhandler = GenericError(name=self.keywords['name'],program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])
        else:
            allowed_keys={'name','program_keys','structure'}
            self.checker = GenericChecker(name=self.keywords['name'],program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])
            self.errhandler = GenericError(name=self.keywords['name'],program_keys=self.keywords['program_keys'],structure=self.keywords['structure'])

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
    
    def close_logger(self):
        """Close logger handlers 
            (prevents IOError from too many handlers being open)
        """
        handlerlist = list(self.logger.handlers) #TTM 428; also deleted return after OSError
        for myhandler in handlerlist:
            self.logger.removeHandler(myhandler)
            myhandler.flush()
            myhandler.close()
        return

    def is_complete(self):
        '''Function to check if Ingredient is ready'''
        if not self.checker.is_started():
            return False #hasn't started running yet
        complete = self.checker.is_complete()
        frozen = self.checker.is_frozen()
        if complete or frozen:
            errct = self.errhandler.loop_through_errors()
            if errct > 0:
                if os.path.isfile(os.path.join(self.keywords['name'],'error.5.tar.gz')):
                    self.logger.error("Ingredient directory already has 5 error zip files. A manual look is required.")
                    self.change_my_status("E")
                    return False
                if 'mast_auto_correct' in self.keywords['program_keys'].keys():
                    if str(self.keywords['program_keys']['mast_auto_correct']).strip()[0].lower() == 'f':
                        self.change_my_status("E")
                    else:
                        self.change_my_status("S")
                else:
                    self.change_my_status("S")
                self.errhandler.clean_up_directory()
                return False
            else:
                if complete:
                    self.metafile.write_data('completed on', time.asctime())
                    if 'get_energy_from_energy_file' in dirutil.list_methods(self.checker,0):
                        energy = self.checker.get_energy_from_energy_file()
                        self.metafile.write_data('energy', energy)
                    return True
                else:
                    return False
        else:
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
        from MAST.submit import script_commands
        script_commands.write_submit_script(self.keywords)
        return
    
    def run(self, mode='serial', curdir=os.getcwd()):
        if "dagman" in dirutil.get_mast_platform():
            return
        from MAST.submit import queue_commands 
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

    def change_my_status(self, newstatus):
        """Change an ingredient status by writing the new status to 
            change_status.txt in the ingredient folder, to get picked
            up by the recipe plan.
            Args:
                newstatus <str>: New status to which to change the ingredient.
        """
        ingdir = self.keywords['name']
        oneup = os.path.dirname(ingdir)
        tryrecipe = os.path.basename(oneup)
        statuspath = ""
        if dirutil.dir_is_in_scratch(tryrecipe):
            statuspath = "%s/change_status.txt" % ingdir
        else:
            twoup = os.path.dirname(oneup)
            tryrecipe = os.path.basename(twoup)
            if dirutil.dir_is_in_scratch(tryrecipe):
                statuspath = "%s/change_status.txt" % oneup
            else:
                raise MASTError(self.__class__.__name__, "Cannot change status of ingredient %s as recipe %s or %s is not found in $MAST_SCRATCH." % (self.keywords['name'],oneup, twoup))
        if os.path.isfile(statuspath):
            statusfile = MASTFile(statuspath)
        else:
            statusfile=MASTFile()
        statusfile.data.append("%s:recommend:%s" % (newstatus, time.asctime()))
        statusfile.to_file(statuspath)
        self.logger.info("Recommending status change to %s" % newstatus)

    def uses_atom_indexing(self):
        mysi = os.path.join(os.path.dirname(self.keywords['name']),
                "structure_index_files")
        if os.path.isdir(mysi):
            return True
        return False

    def create_atom_index_object(self):
        myai = AtomIndex(structure_index_directory=os.path.join(os.path.dirname(self.keywords['name']),'structure_index_files'))
        return myai
