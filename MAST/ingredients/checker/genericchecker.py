##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
from pymatgen.io.vasp import Poscar
from pymatgen.core.structure import Structure
from MAST.utility import dirutil
from MAST.utility.mastfile import MASTFile
from MAST.utility import MASTError
from MAST.ingredients.checker import BaseChecker
from MAST.ingredients.checker import VaspChecker
import os
import logging
import pymatgen
import numpy as np
import time
import shutil
class GenericChecker(BaseChecker):
    """Generic checker functions
        The generic program should accept an input file named "input.txt"
        with separate lines of the format:
            keyword <delimiter> value
        The delimiter should be set in mast_input_delimiter:
            mast_input_delimiter =
            The default is mast_input_delimiter None, which corresponds to a space.
        The generic program should be run as:
            <mast_exec value>
            for example, "python myscript.py input.txt"
        The generic starting signal should be given by the presence of a file:
            mast_started_file <filename>
        The generic completion signal should be given in the ingredient
        subsection in the $ingredients section of the input file as:
            mast_complete_file <filename>
            mast_complete_search <string to search for, for completion>
                              OR, use the value None in order to check
                              just the existence of the file.
                              Defaults to None.
    """
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseChecker.__init__(self, allowed_keys, **kwargs)

    def is_complete(self):
        """Check for the existence of mast_complete_file (default), or
            for a certain string, given by mast_complete_search"""
        if not 'mast_complete_file' in self.keywords['program_keys'].keys():
            raise MASTError(self.__class__.__name__, "No completion file indicated by mast_complete_file keyword for %s. Cannot determine whether run is complete." % self.keywords['name'])
            return False
        checkfile = os.path.join(self.keywords['name'],self.keywords['program_keys']['mast_complete_file'])
        searchstr = None
        if 'mast_complete_search' in self.keywords['program_keys'].keys():
            searchstr = self.keywords['program_keys']['mast_complete_search']
        if searchstr == "None":
            searchstr = None
        if os.path.isfile(checkfile):
            if searchstr == None:
                return True
            else:
                tempopen = MASTFile(checkfile)
                mymatch = tempopen.get_last_line_match(searchstr)
                if mymatch == None:
                    return False
                else:
                    return True
        else:
            return False

    def is_ready_to_run(self):
        """Generic program is ready to run if the input.txt file
            has been written, and if any files in 
            mast_copy_files are present.
            The mast_copy_files keyword must be a space-delimited list of full file paths.
        """
        dirname = self.keywords['name']
        notready=0
        if not(os.path.isfile(dirname + "/input.txt")):
            notready = notready + 1
        if 'mast_copy_files' in self.keywords['program_keys'].keys():
            myfiles = self.keywords['program_keys']['mast_copy_files'].split()
            for myfile in myfiles:
                fname = os.path.basename(myfile)
                if not os.path.isfile(os.path.join(dirname, fname)):
                    notready = notready + 1
        if notready > 0:
            return False
        else:
            return True

    def _copy_over_any_files(self):
        """
            The mast_copy_files ingredient keyword must be a 
            space-delimited list of full file paths:
            mast_copy_files //home/user/f1 //home/user/f2
        """
        dirname = self.keywords['name']
        if 'mast_copy_files' in self.keywords['program_keys'].keys():
            myfiles = self.keywords['program_keys']['mast_copy_files'].split()
            for myfile in myfiles:
                fname = os.path.basename(myfile)
                if not os.path.isfile(os.path.join(dirname, fname)):
                    shutil.copy(myfile, os.path.join(dirname, fname))

    def _phon_poscar_setup(self):
        """Set up a PHON POSCAR file. Strip out the "elements" line (that is,
            use VASP version 4 format. Also strip out anything beneath the atoms
            line.
        """
        name = self.keywords['name']
        pospath = os.path.join(name, "POSCAR")
        prepath = os.path.join(name, "POSCAR_prePHON")
        if os.path.isfile(pospath): #Already done. Return.
            return
        my_poscar = Poscar.from_file(prepath) 
        my_poscar.selective_dynamics=None #unset SD if it is set
        my_poscar.velocities=None #unset velocities
        dirutil.lock_directory(name)
        my_poscar.write_file(pospath)
        dirutil.unlock_directory(name)
        #pick up a copy and strip out the elements line.
        mypfile = MASTFile(pospath)
        myline6=mypfile.get_line_number(6)
        if myline6.strip().split()[0].isalpha:
            mypfile.modify_file_by_line_number(6,"D")
        mypfile.to_file(pospath)
        return



    def _input_get_non_mast_keywords(self):
        """Sort out the non-mast keywords and make a dictionary."""
        input_dict=dict()
        #allowedpath = os.path.join(dirutil.get_mast_install_path(), 'MAST',
        #                'ingredients','programkeys','phon_allowed_keywords.py')
        #allowed_list = self._phon_inphon_get_allowed_keywords(allowedpath)
        for key, value in self.keywords['program_keys'].iteritems():
            if not key[0:5] == "mast_":
                input_dict[key]=value
        return input_dict

    def _input_setup(self):
        """Set up the input.txt file. Note that this is not the recipe
            input file, it is the ingredient input file.
        """
        name=self.keywords['name']
        myd = dict()
        myd = self._input_get_non_mast_keywords()
        my_input = MASTFile()
        if not 'mast_delimiter' in self.keywords['program_keys'].keys():
            delim = " "
        else:
            delim = self.keywords['program_keys']['mast_delimiter']
            if delim == "None":
                delim = " "
        for key, value in myd.iteritems():
            my_input.data.append(str(key) + delim + str(value) + "\n")
        my_input.to_file(name + "/input.txt")
        self._copy_over_any_files()
        return 

    def set_up_program_input(self):
        """For the generic program, MAST will only set up the input.txt
            file. All other necessary files must be previously set up and
            located in the input.txt file, for example:
                 StartingFile = ~/structurebank/startingstructure.xyz 
            in input.txt. This file will not be created by MAST.
        """
        self._input_setup()
        return
    def is_started(self):
        """See if the ingredient has been started on
            the queue. For the generic program, a signal file must
            be specified in the ingredients keywords as mast_started_file
        """
        if not 'mast_started_file' in self.keywords['program_keys'].keys():
            raise MASTError(self.__class__.__name__, "No file indicated by mast_started_file keyword for %s. Cannot determine whether run has started." % self.keywords['name'])
            return False
        checkfile = os.path.join(self.keywords['name'],self.keywords['program_keys']['mast_started_file'])
        if os.path.isfile(checkfile):
            return True
        else:
            return False
