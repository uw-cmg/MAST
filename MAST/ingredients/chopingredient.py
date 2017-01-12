##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-06-13 by Zhewen Song
##############################################################
import os, re
import numpy as np
import logging
import shutil
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.structure import Lattice
from pymatgen.util.coord_utils import find_in_coord_list
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Chgcar
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import Metadata
from MAST.utility import MASTFile
from MAST.utility import dirutil
#from MAST.utility.defect_formation_energy import DefectFormationEnergy as DFE
from MAST.ingredients.baseingredient import BaseIngredient
from MAST.ingredients.pmgextend.structure_extensions import StructureExtensions
from MAST.ingredients.checker import VaspNEBChecker
from MAST.ingredients.checker import VaspChecker
from MAST.ingredients.pmgextend.atom_index import AtomIndex

class ChopIngredient(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of directory'),
            'program': (str, str(), 'Program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object'),
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)

    def _fullpath_childname(self, childname):
        """Get full path of the child directory.
            Args: 
                childname <str>: child directory
        """
        cshort = os.path.basename(childname)
        recipedir = os.path.dirname(self.keywords['name'])
        tryrecipe = os.path.basename(recipedir)
        if not (dirutil.dir_is_in_scratch(tryrecipe)):
            recipedir = os.path.dirname(recipedir)
            tryrecipe = os.path.basename(recipedir)
            if not (dirutil.dir_is_in_scratch(tryrecipe)):
                raise MASTError(self.__class__.__name__, "Could not find child directory %s in %s" % (childname, recipedir))
        return os.path.join(recipedir, cshort)

    def copy_file_with_prepend(self, copyfrom="", copyto="", childdir="", softlink=0):
        """Duplicate an ingredient file into the child
           directory, with the ingredient 
           name prepended
            e.g. "OSZICAR" becomes "defect_vac1_q=p2_stat_OSZICAR"
            Args:
                copyfrom <str>: File name to copy, e.g. OSZICAR
                copyto <str>: File name to copy to, e.g. OSZICAR
                childdir <str>: Child directory
                softlink <int>: 0 - copy (default)
                                1 - softlink
        """
        ingname = os.path.basename(self.keywords['name'])
        toname = "%s_%s" % (ingname, copyto)
        self.copy_file(copyfrom, toname, childdir, softlink)

    def copy_file_no_name_validation(self, copyfrom="", copyto="", childdir="", softlink=0):
        return self.copy_file(copyfrom, copyto, childdir, softlink, 0)

    def copy_file(self, copyfrom="", copyto="", childdir="", softlink=0, validate=1):
        """Copy a file.
            Args:
                copyfrom <str>: name to copy from, e.g. CONTCAR
                copyto <str>: name to copy to, e.g. POSCAR
                childdir <str>: child directory.
                    If not given, use the same directory.
                softlink <int>: 0 (default) - copy
                                1 - softlink
            Log an error if the file is already found
                and do not copy over.
            Do not copy over if the file is empty.
        """
        mydir = self.keywords['name']
        if copyfrom == "":
            raise MASTError(self.__class__.__name__, "No copy-from file given, ingredient %s" % mydir)
        if copyto == "":
            copyto = copyfrom
        if childdir == "":
            childdir = mydir
        if validate == 1:
            childdir = self._fullpath_childname(childdir)
        if not os.path.isdir(childdir):
            raise MASTError(self.__class__.__name__, "No directory for copying into, at %s" % childdir)
        if childdir == mydir and copyfrom == copyto:
            raise MASTError(self.__class__.__name__, "Copy-from and copy-to are the same file and directory for ingredient %s" % mydir)
        self.logger.info("Attempting to copy %s/%s to %s/%s" % (mydir, copyfrom, childdir, copyto))
        frompath = "%s/%s" % (mydir, copyfrom)
        topath = "%s/%s" % (childdir, copyto)
        if not os.path.isfile(frompath):
            self.logger.error("No file at %s. Skipping copy." % frompath)
        else:
            if os.path.isfile(topath):
                self.logger.error("File found at %s already. Skipping copy." % topath)
            else:
                myfile = MASTFile(frompath)
                if len(myfile.data) == 0:
                    self.logger.error("File at %s is empty. Skipping copy." % frompath)
                else:
                    if softlink == 1:
                        curpath = os.getcwd()
                        os.chdir(childdir)
                        mylink=subprocess.Popen("ln -s %s/%s %s" % (mydir, copyfrom, copyto), shell=True)
                        mylink.wait()
                        os.chdir(curpath)
                        self.logger.info("Softlinked file from %s to %s" % (frompath, topath))
                    else:
                        shutil.copy(frompath, topath)
                        self.logger.info("Copied file from %s to %s" % (frompath, topath))
        return

    def get_saddle_dir(self):
        """Get the saddle directory of an NEB calculation.
            VASP only supported.
        """
        if not 'vasp' in self.program:
            raise MASTError(self.__class__.__name__, "Program %s not supported" % self.program)
        myname = self.keywords['name']
        subdirs = dirutil.immediate_subdirs(myname)
        highenergystr="Not evaluated"
        highenergy=-1000000000000.0
        for subdir in subdirs:
            fullsub = os.path.join(myname, subdir)
            singlechecker = VaspChecker(name=fullsub, program_keys = dict(self.checker.keywords['program_keys']), structure = self.checker.keywords['structure'].copy())
            singlechecker.keywords['program_keys']['mast_program'] = "vasp"
            myenergy = singlechecker.get_energy_from_energy_file()
            if myenergy > highenergy:
                highenergystr = subdir
                highenergy = myenergy
        self.logger.info("Saddle directory: %s" % highenergystr)
        return highenergystr

    def copy_saddle_file(self, oldfname, newfname, childname=""):
        """Forward a file from the saddle of a NEB calculation.
            VASP calculations only.
            Args:
                oldfname <str>: Old file name
                newfname <str>: New file name
                childname <str>: Child directory
        """
        saddledir = self.get_saddle_dir()
        self.copy_file("%s/%s" % (saddledir, oldfname), newfname, childname)
        return
    def copy_saddle_file_with_prepend(self, oldfname, newfname, childname=""):
        """Forward a file from the saddle of a NEB calculation with
            the original NEB prepend, but not the inage name.
            VASP calculations only.
            Args:
                oldfname <str>: Old file name
                newfname <str>: New file name
                childname <str>: Child directory
            e.g. 
            self.copy_saddle_file_with_prepend("CONTCAR","POSCAR", childdir)
            neb_vac1-vac2_stat/03/CONTCAR will be copied into childdir as:
            neb_vac1-vac2_stat_POSCAR
        """
        saddledir = self.get_saddle_dir()
        ingname = os.path.basename(self.keywords['name'])
        prependedname = "%s_%s" % (ingname, newfname)
        self.copy_file("%s/%s" % (saddledir, oldfname), prependedname, childname)
        return
    
    def softlink_file(self, linkfrom="", linkto="", childdir=""):
        """Softlink a file.
            Args:
                linkfrom <str>: name to link from, e.g. CONTCAR
                linkto <str>: name to link to, e.g. POSCAR
                childdir <str>: child directory.
            Log an error if the file is already found
                and do not copy over.
        """
        return self.copy_file(linkfrom, linkto, childdir, 1)

    def file_exists(self, filename=""):
        """Check for the existence of a file in the ingredient
            folder.
            Args:
                filename <str>: File name (e.g. OUTCAR)
            Returns:
                True if the file exists; False otherwise
        """
        checkfile = os.path.join(self.keywords['name'], filename)
        filefound = os.path.isfile(checkfile)
        if not filefound:
            self.logger.warning("File %s not found." % checkfile)
        return filefound

    def file_has_string(self, filename="", searchstring="", last_x_lines=0):
        """Check for a string in the file.
            Args:
                filename <str>: File name (e.g. OUTCAR)
                searchstring <str>: Search string
                last_x_lines <str, will be converted to int>:
                    Optional: search only the last X lines.
                    0 (default) - search all lines. 
            Returns:
                True if the string is found; False otherwise
        """
        if not (self.file_exists(filename)):
            return False
        if searchstring == "":
            self.logger.error("No search string given. Returning None.")
            return None
        checkfile = os.path.join(self.keywords['name'],filename)
        tempopen = MASTFile(checkfile)
        mymatch = tempopen.get_last_x_lines_line_match(searchstring, last_x_lines)
        if mymatch == None:
            return False
        else:
            return True

    def copy_fullpath_file(self, copyfromfullpath="", copyto="", childdir="", softlink=0):
        """Copy a file.
            Args:
                copyfromfullpath <str>: full path of the file
                    to copy from, e.g.
                    //home/user/POSCAR_allstart
                copyto <str>: name to copy to, e.g. POSCAR
                childdir <str>: child directory.
                    If not given, use the ingredient directory.
                softlink <int>: 0 (default) - copy
                                1 - softlink
            Log an error if the file is already found
                and do not copy over.
            Do not copy over if the file is empty.
        """
        mydir = self.keywords['name']
        if copyfromfullpath == "":
            raise MASTError(self.__class__.__name__, "No copy-from file given, ingredient %s" % mydir)
        if copyto == "":
            copyto = os.path.basename(copyfromfullpath)
        if childdir == "":
            childdir = mydir
        childdir = self._fullpath_childname(childdir)
        if not os.path.isdir(childdir):
            raise MASTError(self.__class__.__name__, "No directory for copying into, at %s" % childdir)
        self.logger.info("Attempting to copy %s to %s/%s" % (copyfromfullpath, childdir, copyto))
        topath = "%s/%s" % (childdir, copyto)
        if not os.path.isfile(copyfromfullpath):
            self.logger.error("No file at %s. Skipping copy." % copyfromfullpath)
        else:
            if os.path.isfile(topath):
                self.logger.error("File found at %s already. Skipping copy." % topath)
            else:
                myfile = MASTFile(copyfromfullpath)
                if len(myfile.data) == 0:
                    self.logger.error("File at %s is empty. Skipping copy." % copyfromfullpath)
                else:
                    if softlink == 1:
                        curpath = os.getcwd()
                        os.chdir(childdir)
                        mylink=subprocess.Popen("ln -s %s %s" % (copyfromfullpath, copyto), shell=True)
                        mylink.wait()
                        os.chdir(curpath)
                        self.logger.info("Softlinked file from %s to %s" % (copyfromfullpath, topath))
                    else:
                        shutil.copy(copyfromfullpath, topath)
                        self.logger.info("Copied file from %s to %s" % (copyfromfullpath, topath))
        return

    def write_ingred_input_file(self, fname="", allowed_file="all", not_allowed_file="none", upperkey=1, delimiter=" "):
        """Write an input file.
            Args: 
                fname <str>: File name for the ingredient input 
                            file to be created, e.g. INCAR
                allowed_file <str>: File name for the list
                    of allowed keywords. Use "all" to allow
                    all non-mast keywords.
                upperkey <str or int>: 
                        1 - uppercase keywords (default)
                        0 - leave keywords their own case
                delimiter <str>: Delimiter to place between
                    keywords and values in the input file.
                    Default is space.
                    Omit this parameter to use a space.
                    Neither semicolon nor comma may be used
                    as a delimiter, as these are special
                    delimiters in the MAST input file already.
        """
        if delimiter == "":
            delimiter = " "
        if allowed_file.lower() == "all":
            okay_keys = self._get_allowed_non_mast_keywords("", upperkey)
        else:
            okay_keys = self._get_allowed_non_mast_keywords(allowed_file, upperkey)
        if not_allowed_file.lower() == "none": bad_keys = list()
        else:
            path = os.path.join(dirutil.get_mast_install_path(),'ingredients','programkeys',not_allowed_file)
            bad_keys = list()
            fp = open(path)
            lines = fp.readlines()
            for i in range(1,len(lines)-1):
                if int(upperkey)==1:
                    bad_keys.append(lines[i].strip().upper())
                else:
                    bad_keys.append(lines[i].strip().lower())
                
        my_input = MASTFile()
        for key, value in okay_keys.iteritems():
            if not key in bad_keys:
                my_input.data.append(str(key) + delimiter + str(value) + "\n")
        inputpath = os.path.join(self.keywords['name'],fname)
        if os.path.isfile(inputpath):
            self.logger.error("File already exists at %s. Skipping input file writing." % inputpath)
        else:
            my_input.to_file(inputpath)
        return 

    def write_ordered_ingred_input_file(self, fname="", allowed_file="all", upperkey=1, delimiter=" "):
        """Write an input file.
            Assumes that keyword dictionary is given as
            ["#.keyword"] = value
            For example ["2.fix"]=["nvt"]
            Creates an input with keywords in order.
            Args: 
                fname <str>: File name for the ingredient input 
                            file to be created, e.g. INCAR
                allowed_file <str>: File name for the list
                    of allowed keywords. Use "all" to allow
                    all non-mast keywords.
                upperkey <str or int>: 
                        1 - uppercase keywords (default)
                        0 - leave keywords their own case
                delimiter <str>: Delimiter to place between
                    keywords and values in the input file.
                    Default is space.
                    Omit this parameter to use a space.
                    Neither semicolon nor comma may be used
                    as a delimiter, as these are special
                    delimiters in the MAST input file already.
        """
        if delimiter == "":
            delimiter = " "
        if allowed_file.lower() == "all":
            okay_keys = self._get_allowed_non_mast_keywords("", upperkey)
        else:
            okay_keys = self._get_allowed_non_mast_keywords(allowed_file, upperkey)
        keylist=list()
        for key, value in okay_keys.iteritems():
            key = str(key)
            keynum = int(key.split(".")[0])
            keywordval = key.split(".")[1]
            value = str(value)
            keylist.append([keynum, keywordval, value])
        keylist.sort()
        my_input = MASTFile()
        for keytriplet in keylist:
            my_input.data.append(keytriplet[1] + delimiter + keytriplet[2] + "\n")
        inputpath = os.path.join(self.keywords['name'],fname)
        if os.path.isfile(inputpath):
            self.logger.error("File already exists at %s. Skipping input file writing." % inputpath)
        else:
            my_input.to_file(inputpath)
        return 

    def _get_allowed_non_mast_keywords(self, allowed_file="", upperkey=1):
        """Get the non-mast keywords and make a dictionary.
            Args:
                allowed_file <str>: File name containing
                    allowed keywords.
                    If this argument is not entered,
                    then any non-mast keyword will be allowed.
                upperkey <str or int>: 
                        1 - uppercase keywords (default)
                        0 - leave keywords their own case
            Returns:
                my_dict <dict>: Dictionary of keyword and
                    value entries
        """
        my_dict=dict()
        allowed_list = list()
        if not(allowed_file == ""):
            if os.path.isfile(allowed_file):
                allowed_list = self._get_allowed_keywords(allowed_file, upperkey)
            else:
                allowedpath = os.path.join(dirutil.get_mast_install_path(),
                            'ingredients','programkeys',allowed_file)
                allowed_list = self._get_allowed_keywords(allowedpath, upperkey)
        for key, value in self.keywords['program_keys'].iteritems():
            if not key[0:5] == "mast_":
                if int(upperkey) == 1:
                    keytry = key.upper()
                else:
                    keytry = key
                keyokay = None
                if len(allowed_list) > 0:
                    if keytry in allowed_list:
                        keyokay = True
                    else:
                        keyokay = False
                        self.logger.warning("Ignoring program key %s for INCAR. To allow this keyword, add it to %s" % (keytry, allowedpath))
                else:
                    keyokay = True
                if keyokay:
                    if type(value)==str and value.lower() in ['false','true']:
                        my_dict[keytry]=value.capitalize() #First letter cap
                    else:
                        my_dict[keytry]=value
        return my_dict
    
    def _get_allowed_keywords(self, allowedpath, upperkey=1):
        """Get allowed keywords.
            Args:
                allowedpath <str>: file path for allowed keywords
                upperkey <int or str>: uppercase allowed keywords
        """
        if not os.path.isfile(allowedpath):
            self.logger.error("No file at %s for allowed keywords. Returning empty list." % allowedpath)
            return list()
        allowed = MASTFile(allowedpath)
        allowed_list=list()
        for entry in allowed.data:
            if upperkey == 1:
                allowed_list.append(entry.strip().upper())
            else:
                allowed_list.append(entry.strip())
        return allowed_list

    def no_setup(self):
        """No setup is needed."""
        return
    def no_update(self, childname=""):
        """No updating is needed."""
        return
    
    def run_command(self, commandstr=""):
        """Run a command. Command needs to be of the form
            "<str>.py arg1 arg2 ...." and cannot contain
            the following characters: *, |, &, <, >
        """
        oopslist = ["*","|","&","<",">"]
        for oopchar in oopslist:
            if oopchar in commandstr:
                self.logger.error("Cannot run command specified by command string (see input file) because of bad character.")
                return 
        stringsplit = commandstr.split()
        if not stringsplit[0].strip()[-3:] == ".py":
            self.logger.error("Cannot run command specified by command string (see input file) because of non-python file.")
            return 
        else:
            import subprocess
            myproc = subprocess.Popen(commandstr, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            myresults = myproc.communicate()[0]
            return myresults

    def write_neb(self):
        """Get the parent structures, sort and match atoms, and interpolate.
            Write images to the appropriate folders.
        """
        parentstructures = self.get_parent_structures()
        parentimagestructures = self.get_parent_image_structures()
        image_structures = list()
        if len(parentimagestructures) == 0:
            sxtend = StructureExtensions(struc_work1=parentstructures[0], struc_work2=parentstructures[1],name=self.keywords['name'])
            image_structures = sxtend.do_interpolation(self.keywords['program_keys']['mast_neb_settings']['images'])
            if os.path.exists(os.path.join(os.path.dirname(self.keywords['name']),'structure_index_files')):
                image_structures_raw = list(image_structures)
                image_structures = list()
                self.logger.info("Attempt to unsort.")
                myai = AtomIndex(structure_index_directory=os.path.join(os.path.dirname(self.keywords['name']),'structure_index_files'))
                manifestep1=myai.guess_manifest_from_ingredient_metadata(self.keywords['name'],0)
                manifestep2=myai.guess_manifest_from_ingredient_metadata(self.keywords['name'],1)
                myai.make_temp_manifest_from_scrambled_structure(self.keywords['name'],image_structures_raw[0],os.path.join(self.keywords['name'],'scrambledep1'))
                myai.make_temp_manifest_from_scrambled_structure(self.keywords['name'],image_structures_raw[-1],os.path.join(self.keywords['name'],'scrambledep2'))
                for sidx in range(0,len(image_structures_raw)-1):
                    onestruc = image_structures_raw[sidx]
                    newstruc = myai.unscramble_a_scrambled_structure(self.keywords['name'], onestruc, manifestep1, os.path.join(self.keywords['name'],"scrambledep1"))
                    image_structures.append(newstruc)
                newstruc = myai.unscramble_a_scrambled_structure(self.keywords['name'], image_structures_raw[-1], manifestep2, os.path.join(self.keywords['name'],"scrambledep2"))
                image_structures.append(newstruc)


        else:
            image_structures.append(parentstructures[0])
            image_structures.extend(parentimagestructures)
            image_structures.append(parentstructures[1])
        if image_structures == None:
            raise MASTError(self.__class__.__name__,"Bad number of images")
        self.checker.keywords['program_keys']['mast_neb_settings']['image_structures']=image_structures
        self.checker.set_up_program_input()
        self.place_parent_energy_files()
        self.write_submit_script()
        return


    def write_pathfinder_neb(self, chgcarfolder, addlsites=0):
        """Get the parent structures, sort and match atoms, and interpolate
            using Daniil Kitchaev's NEB Pathfinder
            Write images to the appropriate folders.
            Args:
                chgcarfolder <str>: CHGCAR folder, e.g. a perfect static
                addlsites <int>: number of additional sites at the end of
                            the NEB manifest to pick up for needing
                            path finding, e.g. vacancy positions
        """
        parentstructures = self.get_parent_structures()
        parentimagestructures = self.get_parent_image_structures()
        image_structures = list()
        if len(parentimagestructures) > 0:
            raise MASTError(self.__class__.__name__, "Can only use write_pathfinder_neb without mast_coordinates, or on initial NEB setup.")
        


        #sxtend = StructureExtensions(struc_work1=parentstructures[0], struc_work2=parentstructures[1],name=self.keywords['name'])
        #image_structures = sxtend.do_interpolation(self.keywords['program_keys']['mast_neb_settings']['images'])

        s1 = parentstructures[0]
        s2 = parentstructures[1]
        #s1 = Poscar.from_file(args.s1).structure
        #s2 = Poscar.from_file(args.s2).structure
        
        #Must be an interstitial
        mysi = os.path.join(os.path.dirname(self.keywords['name']),
                "structure_index_files")
        if not os.path.isdir(mysi):
            raise MASTError(self.__class__.__name__, "Can only use write_pathfinder_neb with structure indexing turned on.")
        
        #get aidxs from ;int lines of defect manifests
        #match aidx to position in sorted NEB manifests
        #get those sites

        # Find diffusing species site indices
        #mg_sites = []
        #for site_i, site in enumerate(s1.sites):
        #    if site.specie == Element("Mg"):
        #        mg_sites.append(site_i)
        sxtend = StructureExtensions(struc_work1=parentstructures[0], struc_work2=parentstructures[1],name=self.keywords['name'])
        pmg_image_structures = sxtend.do_interpolation(self.keywords['program_keys']['mast_neb_settings']['images'])
        
        myai=AtomIndex(structure_index_directory=mysi)
        manifestep1=myai.guess_manifest_from_ingredient_metadata(self.keywords['name'],0)
        manifestep2=myai.guess_manifest_from_ingredient_metadata(self.keywords['name'],1)
        myai.make_temp_manifest_from_scrambled_structure(self.keywords['name'],pmg_image_structures[0],os.path.join(self.keywords['name'],'pmg_scrambledep1'))
        myai.make_temp_manifest_from_scrambled_structure(self.keywords['name'],pmg_image_structures[-1],os.path.join(self.keywords['name'],'pmg_scrambledep2'))
        
        nebman=myai.guess_manifest_from_ingredient_metadata(self.keywords['name'],0)
        nebsplit=nebman.rsplit("_", 1)
        defectman="%s_" % nebsplit[0]
        intlist=list()
        defectmanlist=myai.read_manifest_file(defectman)
        for defectidx in defectmanlist:
            if "int" in defectidx:
                intlist.append(defectidx.split(';')[0])
        
        nebmanlist=myai.read_manifest_file(nebman)
        pmg_nebmanlist=myai.read_manifest_file(os.path.join(self.keywords['name'],'pmg_scrambledep1'))
        addlsites=int(addlsites)
        if addlsites > 0:
            self.logger.info("additional sites: %s" % addlsites)
            intlist.extend(nebmanlist[-1*addlsites:])

        intsites=list()
        for intitem in intlist:
            intidx = pmg_nebmanlist.index(intitem)
            intsites.append(intidx)

        print "ADDLSITES: %s, INTSITES: %s" % (addlsites, intsites)
        print "S1:"
        print s1
        print "S2:"
        print s2
        self.logger.info("intsites: %s" % intsites)
        # Interpolate
        #use a param for perf_stat or othe ing label for chgcar
        #print("Using CHGCAR potential mode.")
        #chg = Chgcar.from_file(args.chg)
        #pf = NEBPathfinder(s1, s2, relax_sites=mg_sites, v=ChgcarPotential(chg).get_v(), n_images=10)
        
        if not os.path.isfile("%s/CHGCAR" % chgcarfolder):
            chgcarfolder = self._fullpath_childname(chgcarfolder)
            if not os.path.isfile("%s/CHGCAR" % chgcarfolder):
                raise MASTError(self.__class__.__name__, "No CHGCAR in %s " % chgcarfolder)
        if not (os.stat("%s/CHGCAR" % chgcarfolder).st_size > 1):
            raise MASTError(self.__class__.__name__, "CHGCAR in %s appears to be empty" % chgcarfolder)
        chg = Chgcar.from_file("%s/CHGCAR" % chgcarfolder)
        
        numim = self.keywords['program_keys']['mast_neb_settings']['images']
        
        from MAST.utility.daniil_pathfinder import NEBPathfinder
        from MAST.utility.daniil_pathfinder import ChgcarPotential
        pf = NEBPathfinder(s1, s2, relax_sites=intsites, v=ChgcarPotential(chg).get_v(), n_images=numim+1)

        image_structures=pf.images
        
        image_structures_raw = list(image_structures)
        image_structures = list()
        self.logger.info("Attempt to unsort.")
        myai = AtomIndex(structure_index_directory=os.path.join(os.path.dirname(self.keywords['name']),'structure_index_files'))
        manifestep1=myai.guess_manifest_from_ingredient_metadata(self.keywords['name'],0)
        manifestep2=myai.guess_manifest_from_ingredient_metadata(self.keywords['name'],1)
        for sidx in range(0,len(image_structures_raw)-1):
            onestruc = image_structures_raw[sidx]
            newstruc = myai.unscramble_a_scrambled_structure(self.keywords['name'], onestruc, manifestep1, os.path.join(self.keywords['name'],"pmg_scrambledep1"))
            image_structures.append(newstruc)
        newstruc = myai.unscramble_a_scrambled_structure(self.keywords['name'], image_structures_raw[-1], manifestep2, os.path.join(self.keywords['name'],"pmg_scrambledep2"))
        image_structures.append(newstruc)

        if image_structures == None:
            raise MASTError(self.__class__.__name__,"Bad number of images")
        self.checker.keywords['program_keys']['mast_neb_settings']['image_structures']=image_structures
        self.checker.set_up_program_input()
        self.place_parent_energy_files()
        self.write_submit_script()
        return



        
    def get_parent_structures(self):
        """Assume that parents have written files
            named 'parent_structure_<N>'. 
            For VASP these are CONTCAR-type files.
            Returns:
                [struct_init, struct_fin]: pymatgen Structure objects
        """
        header = os.path.join(self.keywords['name'], "parent_structure_")
        mylabel = BaseIngredient.get_my_label(self, "neb_label").split("-")
        pfpath_init = header + mylabel[0]
        pfpath_fin = header + mylabel[1]
        if not os.path.isfile(pfpath_init):
            raise MASTError(self.__class__.__name__,
                "Error: no parent file at" + pfpath_init)
        if not os.path.isfile(pfpath_fin):
            raise MASTError(self.__class__.__name__,
                "Error: no parent file at" + pfpath_fin)
        struct_init = self.checker.get_structure_from_file(pfpath_init)
        struct_fin = self.checker.get_structure_from_file(pfpath_fin)
        if self.uses_atom_indexing():
            myai = self.create_atom_index_object()
            iname = self.keywords['name']
            nebpcs=[0,1]
            nebstructs=[struct_init, struct_fin]
            scramblenames=['scrambledep_init','scrambledep_fin']
            for nidx in [0,1]:
                sname = os.path.join(iname, scramblenames[nidx])
                manifestep=myai.guess_manifest_from_ingredient_metadata(iname,
                                nebpcs[nidx])
                myai.make_temp_manifest_from_scrambled_structure(iname, 
                                nebstructs[nidx], sname)
                onestruc = nebstructs[nidx]
                newstruc = myai.unscramble_a_scrambled_structure(iname,
                            onestruc, manifestep, sname)
                if nidx == 0:
                    sorted_init = newstruc
                elif nidx == 1:
                    sorted_fin = newstruc
        else:
            base_struct = self.keywords['structure']
            myimages = self.keywords['program_keys']['mast_neb_settings']['images']
            neblines = self.keywords['program_keys']['mast_neb_settings']['lines']
            sxtendi = StructureExtensions(struc_work1 = struct_init, struc_init = base_struct, name=self.keywords['name'])
            sorted_init = sxtendi.sort_structure_and_neb_lines(neblines, '00', myimages) 
            sxtendf = StructureExtensions(struc_work1 = struct_fin, struc_init = base_struct, name=self.keywords['name'])
            sorted_fin = sxtendf.sort_structure_and_neb_lines(neblines, str(myimages + 1).zfill(2), myimages)
        return [sorted_init, sorted_fin]

    def get_parent_image_structures(self):
        """A low-mesh NEB may have written files
            named 'parent_structure_<N-N>_0N'. 
            For VASP these are CONTCAR-type files.
            Returns:
                list of <Structure>: list of pymatgen Structure objects
        """
        header = "parent_structure_"
        numim = self.keywords['program_keys']['mast_neb_settings']['images']
        imct = 1
        imstrs=list()
        if self.uses_atom_indexing():
            myai = self.create_atom_index_object()
            iname = self.keywords['name']
            manifestep=myai.guess_manifest_from_ingredient_metadata(iname,0) #images always use initial endpoint manifest
        while imct <= numim:
            pfpath=""
            for myfile in os.listdir(self.keywords['name']):
                if (header in myfile) and (str(imct).zfill(2) in myfile):
                    pfpath = os.path.join(self.keywords['name'],myfile)
            if pfpath == "":
                pass
            else:
                struct_im = self.checker.get_structure_from_file(pfpath)
                if self.uses_atom_indexing():
                    sname = os.path.join(iname, "scrambledim_%s" % str(imct.zfill(2))) 
                    myai.make_temp_manifest_from_scrambled_structure(iname, 
                                        struct_im, sname)
                    onestruc = struct_im
                    newstruc = myai.unscramble_a_scrambled_structure(iname,
                                    onestruc, manifestep, sname)
                    imstrs.append(newstruc)
                else:
                    base_struct = self.keywords['structure']
                    neblines = self.keywords['program_keys']['mast_neb_settings']['lines']
                    sxtend = StructureExtensions(struc_work1 = struct_im, struc_init = base_struct, name=self.keywords['name'])
                    sorted_im = sxtend.sort_structure_and_neb_lines(neblines, str(imct).zfill(2), self.keywords['program_keys']['mast_neb_settings']['images']) 
                    imstrs.append(sorted_im)
            imct = imct + 1
        if len(imstrs) > 0 and not (len(imstrs) == numim):
            raise MASTError(self.__class__.__name__, "Incomplete number of forwared images found!")
        return imstrs

    def place_parent_energy_files(self):
        """Assume parents have written files parent_energy_<N>.
            Copy these files into the 00 and 0N directories.
        """
        header = os.path.join(self.keywords['name'], "parent_energy_")
        mylabel=BaseIngredient.get_my_label(self, "neb_label").split("-")
        pfpath1= header + mylabel[0]
        pfpath2= header + mylabel[1]
        pffile1=MASTFile(pfpath1)
        pffile2=MASTFile(pfpath2)
        pffile1.to_file(self.checker.get_path_to_write_neb_parent_energy(1)) #MASTFile contains directory locking.
        pffile2.to_file(self.checker.get_path_to_write_neb_parent_energy(2))
        return



    def write_neb_subfolders(self):
        """Write the static runs to each subfolder.
        """
        myname=self.keywords['name']
        mystr=self.keywords['structure']
        numim = int(self.keywords['program_keys']['mast_neb_settings']['images'])
        singlechecker = self.checker
        if self.program == 'vasp':
            nebchecker = VaspNEBChecker(name=self.checker.keywords['name'], program_keys = dict(self.checker.keywords['program_keys']), structure = self.checker.keywords['structure'].copy())
        else:
            raise MASTError(self.__class__.__name__, "NEB checker not supported for program %s") % self.program
        self.checker = nebchecker
        self.write_neb() #Write the POSCAR and OSZICAR files from existing parent-given structures
        self.checker = singlechecker
        for imct in range(1, numim+1):
            subdir = str(imct).zfill(2)
            newname = os.path.join(myname, subdir)
            self.checker.keywords['name']=newname
            self.checker.set_up_program_input()
            self.keywords['name'] = newname
            self.write_submit_script()
        self.keywords['name'] = myname
        self.checker.keywords['name']=myname
        self.checker.keywords['structure']=mystr
        return
    def write_singlerun(self):
        self.checker.set_up_program_input()
        self.write_submit_script()
    def write_singlerun_automesh(self):
        if not self.checker.has_starting_structure_file():
            raise MASTError(self.__class__.__name__,"No starting structure file for ingredient %s; write_singlerun_automesh cannot be used." % self.keywords['name']) 
        if not 'mast_kpoint_density' in self.keywords['program_keys'].keys():
            raise MASTError(self.__class__.__name__,"No mast_kpoint_density ingredient keyword found for ingredient %s; write_singlerun_automesh cannot be used." % self.keywords['name'])
        self.checker.scale_mesh(self.checker.get_initial_structure_from_directory(self.keywords['name']),int(self.keywords['program_keys']['mast_kpoint_density']))
        self.checker.set_up_program_input()
        self.write_submit_script()
    def write_phonon_multiple(self):
        """Write the multiple phonon files, one for each atom and each direction.
            #TTM add atom index
        """
        self.checker.set_up_program_input()
        self.write_submit_script()
        mystructure = self.checker.get_initial_structure_from_directory()
        if self.atomindex:
            sdarrlist=self.atomindex.get_sd_array(self.keywords['name'], True)
        else:
            [pcs,pcr,thresh] = self.get_my_phonon_params()
            sxtend = StructureExtensions(struc_work1 = mystructure, name=self.keywords['name'])
            sdarrlist = sxtend.get_multiple_sd_array(pcs, pcr, thresh)
        if sdarrlist == None:
            raise MASTError(self.__class__.__name__, "No phonons to run!")
        sct=1
        myname=self.keywords['name']
        for sdarr in sdarrlist:
            newname = os.path.join(myname,"phon_%s" % str(sct).zfill(2))
            try:
                os.mkdir(newname)
            except OSError:
                pass
            self.checker.keywords['name']=newname
            self.checker.keywords['structure']=mystructure
            self.checker.set_up_program_input()
            self.checker.add_selective_dynamics_to_structure_file(sdarr)
            self.keywords['name'] = newname
            self.write_submit_script()
            self.checker.keywords['name']=myname
            self.checker.softlink_charge_density_file(newname)
            self.checker.softlink_wavefunction_file(newname)
            sct = sct + 1
        self.checker.keywords['name'] = myname
        self.keywords['name']=myname
        

    def write_phonon_single(self):
        """Write the phonon files to a directory.
            #TTM add atom index
        """
        self.checker.set_up_program_input()
        self.write_submit_script()
        if self.atomindex:
            sdarr=self.atomindex.get_sd_array(self.keywords['name'])
        else:
            mystructure = self.checker.get_initial_structure_from_directory()
            [pcs,pcr,thresh] = self.get_my_phonon_params()
            sxtend = StructureExtensions(struc_work1 = mystructure, name=self.keywords['name'])
            sdarr = sxtend.get_sd_array(pcs, pcr,thresh)
            if sum(sum(sdarr)) == 0:
                raise MASTError(self.__class__.__name__,"Selective dynamics array was all false; no selective dynamics added for %s" % self.keywords['name'])
                return
        self.checker.add_selective_dynamics_to_structure_file(sdarr)

    def get_my_phonon_params(self):
        """Get phonon parameters from 
            ['program_keys']['mast_phonon_settings']['phonon_center_site'] and
            ['program_keys']['mast_phonon_settings']['phonon_center_radius'] 
            ['program_keys']['mast_phonon_settings']['threshold'] 
            Returns:
                [phonon_center_site, phonon_center_radius, threshold]
        """
        myphdict = dict(self.keywords['program_keys']['mast_phonon_settings'])
        return [myphdict['phonon_center_site'],myphdict['phonon_center_radius'],myphdict['threshold']]


    def ready_singlerun(self):
        return BaseIngredient.is_ready_to_run(self)
    def ready_structure(self):
        if self.directory_is_locked():
            return False
        return self.checker.has_starting_structure_file()

    def ready_defect(self):
        self.logger.warning("ready_defect is deprecated. use ready_structure instead")
        return self.ready_structure()

    def ready_neb_subfolders(self):
        """Make sure all subfolders are ready to run."""
        myname=self.keywords['name']
        notready=0
        numim = int(self.keywords['program_keys']['mast_neb_settings']['images'])
        for imct in range(0,numim+2):
            subdir = str(imct).zfill(2)
            newname = os.path.join(myname, subdir)
            self.checker.keywords['name']=newname
            if imct == 0 or imct > numim:
                pass
            else:
                if dirutil.directory_is_locked(newname):
                    notready = notready + 1
                if not self.checker.is_ready_to_run():
                    notready = notready + 1
            imct = imct + 1
        self.checker.keywords['name']=myname
        if notready == 0:
            return True
        else:
            return False

    def ready_subfolders(self):
        """Make sure all subfolders are ready to run."""
        myname=self.keywords['name']
        subdirs = dirutil.walkdirs(myname,1,1)
        notready=0
        if len(subdirs) == 0:
            return False
        for subdir in subdirs:
            newname = os.path.join(myname, subdir)
            self.checker.keywords['name']=newname
            if dirutil.directory_is_locked(newname):
                notready = notready + 1
            if not self.checker.is_ready_to_run():
                notready = notready + 1
        self.checker.keywords['name']=myname
        if notready == 0:
            return True
        else:
            return False

    def run_singlerun(self, mode='serial'):
        return BaseIngredient.run(self, mode)
    def run_neb_subfolders(self):
        """Run all image subfolders."""
        myname=self.keywords['name']
        subdirs = dirutil.walkdirs(myname,1,1)
        imct = 0
        numim = int(self.keywords['program_keys']['mast_neb_settings']['images'])
        for subdir in subdirs:
            if imct == 0 or imct > numim:
                pass
            else:
                newname = os.path.join(myname, subdir)
                self.keywords['name']=newname
                BaseIngredient.run(self,"serial")
            imct = imct + 1
        self.keywords['name']=myname
        return
    def run_subfolders(self):
        """Run all subfolders."""
        myname=self.keywords['name']
        phondirs = dirutil.walkdirs(myname,1,1)
        for phondir in phondirs:
            newname = os.path.join(myname, phondir)
            self.keywords['name']=newname
            BaseIngredient.run(self,"serial")
        self.keywords['name']=myname
        return

    def run_strain(self):
        """Strain the lattice.
            Args:
                Looks for mast_strain in input file, for 
                percent strains.
                mast_strain 1.01 1.01 -0.98
            Returns:
                Creates structure file in directory 
        """
        mystructure = self.checker.get_initial_structure_from_directory()
        mystrain = self.keywords['program_keys']['mast_strain']
        sxtend = StructureExtensions(struc_work1 = mystructure, name=self.keywords['name'])
        strained_structure = sxtend.strain_lattice(mystrain)
        self.checker.write_final_structure_file(strained_structure)
        return

    def run_scale(self):
        """
        """
        try:
            base_structure = self.checker.get_initial_structure_from_directory()
        except: #no initial structure
            base_structure = self.keywords['structure'].copy()
            self.logger.warning("No parent structure detected for induce defect ingredient %s. Using initial structure of the recipe." % self.keywords['name'])
        scalingsize = self.metafile.read_data('scaling_size')
        workdir=os.path.dirname(self.keywords['name'])
        if scalingsize == None: scalingsize = '1 1 1'
        elif '[' and ']' in scalingsize:
            scalingsize = scalingsize.split('[')[1].split(']')[0]
        else: raise MASTError(self.__class__.__name__, "Error in scaling size for the ingredient %s" % self.keywords['name'])
        scalextend = StructureExtensions(struc_work1=base_structure, scaling_size=scalingsize, name=self.keywords['name'])
        scaled = scalextend.scale_structure()
        if self.atomindex:
            mymeta=Metadata(metafile="%s/metadata.txt" % self.keywords['name'])
            scaling_label=mymeta.read_data("scaling_label")
            if scaling_label == None:
                scaling_label = ""
            defect_label=mymeta.read_data("defect_label")
            if defect_label == None:
                defect_label = ""
            manname="manifest_%s_%s_" % (scaling_label, defect_label)
            scaled=self.atomindex.graft_new_coordinates_from_manifest(scaled, manname, "")
            self.logger.info("Getting coordinates from manifest.")
        self.checker.write_final_structure_file(scaled)
        return

    def run_defect(self):
        try:
            base_structure = self.checker.get_initial_structure_from_directory()
        except: #no initial structure
            base_structure = self.keywords['structure'].copy()
            self.logger.warning("No parent structure detected for induce defect ingredient %s. Using initial structure of the recipe." % self.keywords['name'])
        scalingsize = self.metafile.read_data('scaling_size')
        if scalingsize == None: scalingsize = '1,1,1'
        elif '[' and ']' in scalingsize:
            scalingsize = scalingsize.split('[')[1].split(']')[0]        
        else: raise MASTError(self.__class__.__name__, "Error in scaling size for the ingredient %s"%self.keywords['name'])        
        defect = self.keywords['program_keys']['mast_defect_settings']
        scaled = base_structure.copy()
        if self.atomindex:
            mymeta=Metadata(metafile="%s/metadata.txt" % self.keywords['name'])
            scaling_label=mymeta.read_data("scaling_label")
            if scaling_label == None:
                scaling_label = ""
            defect_label=mymeta.read_data("defect_label")
            if defect_label == None:
                raise MASTError(self.__class__.__name__,"Ingredient %s has no defect_label in metadata. Cannot get manifest." % self.keywords['name'])
            parent=mymeta.read_data("parent") 
            manname="manifest_%s_%s_" % (scaling_label, defect_label)
            sxtend = StructureExtensions(struc_work1=scaled, scaling_size=scalingsize, name=self.keywords['name'])
            scaled = sxtend.scale_structure()
            scaled=self.atomindex.graft_new_coordinates_from_manifest(scaled, manname, parent)
            self.logger.info("Getting coordinates from manifest.")
            self.checker.write_final_structure_file(scaled)
            return
        for key in defect:
            if 'subdefect' in key:
                subdefect = defect[key]
                sxtend = StructureExtensions(struc_work1=scaled, scaling_size=scalingsize, name=self.keywords['name'])
                scaled = sxtend.scale_defect(subdefect, defect['coord_type'], defect['threshold'])
            else:
                pass
        self.checker.write_final_structure_file(scaled)
        return

    def run_scale_defect(self):
        try:
            base_structure = self.checker.get_initial_structure_from_directory() 
        except: #no initial structure
            base_structure = self.keywords['structure'].copy()
            self.logger.warning("No parent structure detected for induce defect ingredient %s. Using initial structure of the recipe." % self.keywords['name'])
        scalextend = StructureExtensions(struc_work1=base_structure, name=self.keywords['name'])
        if not 'mast_scale' in self.keywords['program_keys'].keys():
            raise MASTError(self.__class__.__name__,"No mast_scale ingredient keyword for scaling ingredient %s." % self.keywords['name'])
        scaled = scalextend.scale_structure(self.keywords['program_keys']['mast_scale'])
        defect = self.keywords['program_keys']['mast_defect_settings']
        for key in defect:
            if 'subdefect' in key:
                subdefect = defect[key]
                sxtend = StructureExtensions(struc_work1=scaled, struc_work2=base_structure, name=self.keywords['name'])
                scaled = sxtend.scale_defect(subdefect, defect['coord_type'], defect['threshold'])
            else:
                pass
        self.checker.write_final_structure_file(scaled)
        return

    def run_supercell_defect_set(self, parentdir=""):
        """For finite size scaling.
            Given a text file named supercell_list.txt:

            /some/directory/of/parent/ingredient
            ----------------------------------
            ScalingLMN     V_M     Kpoint mesh
            ----------------------------------
            [4, 4, 4] 2.55350475682 [2, 2, 2]
            [3, 4, 3] 3.00247070991 [2, 2, 2]
            [4, 3, 2] 3.05854316622 [2, 2, 4]
            [3, 3, 2] 3.67455399526 [2, 2, 4]
            [2, 2, 2] 5.10700950784 [4, 4, 4]

            There are folders in the parent ingredient
                4x1x1
                3x2x3 etc. within which are a perfect POSCAR
                file and a correspondingly scaled KPOINTS file.
            Copy each folder with POSCAR and KPOINTS.
            Induce the correct defect group for each new folder
                in the current ingredient, scaling by the
                appropriate amount.
            Args:
                parentdir <str>: Parent directory. If not
                        supplied, method will look for
                        supercell_list.txt in the current
                        directory and use the path 
                        specified within the file.
        """
        myname = self.keywords['name']
        supercell_info = MASTFile(os.path.join(myname, "supercell_list.txt"))
        if parentdir == "":
            scaleset_dir = supercell_info.data[0].strip()
        else:
            scaleset_dir = parentdir
        scale_folders = dirutil.immediate_subdirs(scaleset_dir)
        for scale_folder in scale_folders:
            shutil.copytree(os.path.join(scaleset_dir,scale_folder), os.path.join(myname, scale_folder))
        for scale_folder in scale_folders:
            myfolder = os.path.join(myname, scale_folder)
            singlechecker = VaspChecker(name=myfolder,program_keys = dict(self.checker.keywords['program_keys']))
            scaled = singlechecker.get_initial_structure_from_directory() 
            defect = self.keywords['program_keys']['mast_defect_settings']
            scalestring = scale_folder
            for key in defect:
                if 'subdefect' in key:
                    subdefect = defect[key]
                    sxtend = StructureExtensions(struc_work1=scaled, name=self.keywords['name'])
                    scaled = sxtend.scale_defect_by_LMN(scalestring, subdefect, defect['coord_type'], defect['threshold'])
                else:
                    pass
            singlechecker.write_final_structure_file(scaled)
        return True

    def complete_structure(self):
        if self.directory_is_locked():
            return False
        iscomplete = self.checker.has_ending_structure_file()
        if iscomplete:
            if 'update_atom_index_for_complete' in dirutil.list_methods(self.checker,0):
                self.checker.update_atom_index_for_complete()
        return iscomplete
    def complete_singlerun(self):
        iscomplete = BaseIngredient.is_complete(self)
        if iscomplete:
            if 'update_atom_index_for_complete' in dirutil.list_methods(self.checker,0):
                self.checker.update_atom_index_for_complete()
        return iscomplete
    def complete_neb_subfolders(self):
        """Make sure all subfolders are complete."""
        myname=self.keywords['name']
        subdirs = dirutil.walkdirs(myname,1,1)
        if len(subdirs) == 0:
            return False
        notready=0
        imct = 0
        numim = int(self.keywords['program_keys']['mast_neb_settings']['images'])
        for subdir in subdirs:
            newname = os.path.join(myname, subdir)
            self.keywords['name']=newname
            self.checker.keywords['name']=newname
            self.errhandler.keywords['name']=newname
            if imct == 0 or imct > numim:
                pass
            elif not self.is_complete():
                notready = notready + 1
            imct = imct + 1
        self.keywords['name']=myname
        self.checker.keywords['name']=myname
        self.errhandler.keywords['name']=myname
        if notready == 0:
            return True
        else:
            return False

    def complete_subfolders(self):
        """Make sure all subfolders are complete."""
        myname=self.keywords['name']
        phondirs = dirutil.walkdirs(myname,1,1)
        notready=0
        if len(phondirs) == 0:
            return False
        for phondir in phondirs:
            newname = os.path.join(myname, phondir)
            self.keywords['name']=newname
            self.checker.keywords['name']=newname
            self.errhandler.keywords['name']=newname
            if not self.is_complete():
                notready = notready + 1
        self.keywords['name']=myname
        self.checker.keywords['name']=myname
        self.errhandler.keywords['name']=myname
        if notready == 0:
            return True
        else:
            return False
    def complete_supercell_defect_set(self):
        """Check that a supercell defect set is complete
            by making sure that there is an ending structure
            (CONTCAR for VASP) in every subfolder.
        """
        myname = self.keywords['name']
        mysubdirs = dirutil.immediate_subdirs(myname)
        oklist=list()
        if len(mysubdirs) == 0: #No subdirectories set yet
            return False
        for subdir in mysubdirs:
            self.checker.keywords['name']=os.path.join(myname, subdir)
            if self.checker.has_ending_structure_file():
                oklist.append(True)
            else:
                oklist.append(False)
        self.checker.keywords['name']=myname
        okarr = np.array(oklist)
        if okarr.all():
            return True
        else:
            return False

    def give_structure(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)

    def give_neb_structures_to_neb(self, childname):
        """Update to ANOTHER NEB."""
        childname = self._fullpath_childname(childname)
        myct=1
        while myct <= self.keywords['program_keys']['mast_neb_settings']['images']:
            imno = str(myct).zfill(2)
            impath = os.path.join(self.keywords['name'], imno)
            self.checker.keywords['name'] = impath
            self.checker.forward_final_structure_file(childname,"parent_structure_" + BaseIngredient.get_my_label(self, "neb_label") + '_' + imno)
            myct = myct + 1
        return

    def give_supercell_subfolder_file(self, oldfname, newfname, childname):
        """Give each CONTCAR to a corresponding scale1 through scale5
            child folder.
            Args:
                oldfname <str>: Old file name
                newfname <str>: New file name
                childname <str>: Child directory name (fullpath)
        """
        myname = self.keywords['name']
        childbase = os.path.basename(childname)
        if not "scale" in childbase:
            raise MASTError(self.__class__.__name__,"Child directory %s needs 'scale' in its name" % childname)
        namesplit = childbase.split("_")
        for nspl in namesplit:
            if "scale" in nspl:
                scalestr = nspl
                break
        scaleidx = int(scalestr.strip()[5:]) - 1 #e.g. scale3 = index 2
        subdirs = dirutil.immediate_subdirs(myname)
        mysubdir = subdirs[scaleidx]
        self.copy_file(os.path.join(mysubdir, oldfname),newfname,childname)
        


    def give_supercell_defect_kpoints(self, childname):
        """Give each KPOINTS to a corresponding scale1 through scale5
            child folder.
        """
        myname = self.keywords['name']
        if not "scale" in childname:
            raise MASTError(self.__class__.__name__,"Child directory %s needs 'scale' in its name" % childname)
        subdirs = dirutil.immediate_subdirs(myname)





    def give_saddle_structure(self, childname):
        """Forward the middle image structure."""
        childname = self._fullpath_childname(childname)
        myname=self.keywords['name']
        subdirs = dirutil.walkdirs(myname,1,1)
        highenergyimct=0
        imct=0
        highenergy=-1000000000000.0
        #numim = int(self.keywords['program_keys']['images'])
        #middleim = int(numim/2)+1 #returns 1 for 1, 2 for 3 im, 3 for 5 im, etc
        subdirs.sort()
        if self.program == 'vasp_neb' or self.program == 'vasp':
            singlechecker = VaspChecker(name=self.checker.keywords['name'],program_keys = dict(self.checker.keywords['program_keys']),structure = self.checker.keywords['structure'].copy())
        else:
            raise MASTError(self.__class__.__name__, "Program %s not supported" % self.program)
        for subdir in subdirs:
            singlechecker.keywords['name'] = subdir
            myenergy = singlechecker.get_energy_from_energy_file()
            if myenergy > highenergy:
                highenergyimct = imct
                highenergy = myenergy
            imct = imct + 1
        imno = str(highenergyimct).zfill(2)
        impath = os.path.join(self.keywords['name'], imno)
        self.checker.keywords['name'] = impath
        self.checker.forward_final_structure_file(childname)
        self.checker.softlink_charge_density_file(childname)
        self.checker.softlink_wavefunction_file(childname)
        return
    def give_phonon_multiple_forces_and_displacements(self,childname):
        self.checker.combine_dynamical_matrix_files(self.keywords['name'])
        self.checker.combine_displacement_files(self.keywords['name'])
        self.give_phonon_single_forces_and_displacements(childname)
    def give_phonon_single_forces_and_displacements(self, childname):
        #Do NOT forward the CONTCAR structure, since the ending CONTCAR contains a displacement in it. Instead, forward the POSCAR
        childname = self._fullpath_childname(childname)
        self.checker.forward_dynamical_matrix_file(childname)
        self.checker.forward_displacement_file(childname)
        self.checker.forward_initial_structure_file(childname, "POSCAR_prePHON")

    def give_structure_and_energy_to_neb(self, childname):
        childname = self._fullpath_childname(childname)
        label = BaseIngredient.get_my_label(self, "defect_label")
        self.checker.forward_final_structure_file(childname,"parent_structure_" + label)
        self.checker.forward_energy_file(childname, "parent_energy_" + label)
    def give_structure_and_restart_files_softlinks(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)
        self.checker.softlink_charge_density_file(childname)
        self.checker.softlink_wavefunction_file(childname)
    
    def give_structure_and_restart_files(self, childname):
        self.give_structure_and_restart_files_softlinks(childname)

    def give_structure_and_restart_files_full_copies(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)
        self.checker.forward_charge_density_file(childname)
        self.checker.forward_wavefunction_file(childname)
   
    def give_structure_and_charge_density_full_copy(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)
        self.checker.forward_charge_density_file(childname)
    
    def give_structure_and_wavefunction_full_copy(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)
        self.checker.forward_wavefunction_file(childname)

    def give_structure_and_charge_density_softlink(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)
        self.checker.softlink_charge_density_file(childname)
    
    def give_structure_and_wavefunction_softlink(self, childname):
        childname = self._fullpath_childname(childname)
        self.checker.forward_final_structure_file(childname)
        self.checker.softlink_wavefunction_file(childname)

    def give_energy_to_dfe(self, childname):
        childpath = self._fullpath_childname(childname)
        shutil.copy(self.keywords['name']+'/OSZICAR', childpath+'/'+self.keywords['name'].split('/')[-1]+'_OSZICAR')
        shutil.copy(self.keywords['name']+'/CONTCAR', childpath+'/'+self.keywords['name'].split('/')[-1]+'_CONTCAR')
        shutil.copy(self.keywords['name']+'/OUTCAR', childpath+'/'+self.keywords['name'].split('/')[-1]+'_OUTCAR')
        #fp = open(childpath+"/dfe.in", "a")
        #fp.write(self.keywords['name'].split('/')[-1]+":"+open(self.keywords['name']+'/OSZICAR').readlines()[-1].split('E0=')[1].split()[0]+'\n')
        #dfe = DFE(os.path.dirname(childpath))
        #scalingsize = self.metafile.read_data('scaling_size')
        #Ef = dfe._calculate_defect_formation_energies(scalingsize)

    def give_doscar_to_dfe(self, childname):
        childpath = self._fullpath_childname(childname)
        shutil.copy(self.keywords['name']+'/DOSCAR', childpath+'/'+self.keywords['name'].split('/')[-1]+'_DOSCAR')

    def give_outcar_to_dfe(self, childname):
        childpath = self._fullpath_childname(childname)            
        shutil.copy(self.keywords['name']+'/OUTCAR', childpath+'/'+self.keywords['name'].split('/')[-1]+'_OUTCAR')



    
    def give_structure_w_random_displacements(self, childname, disp=0.01):
        import pymatgen as mg
        import random
        disp = float(disp)
        childpath = self._fullpath_childname(childname)
        print childpath,self.keywords['name']
        struct=mg.read_structure(self.keywords['name']+'/CONTCAR')
        for i in range(len(struct)):
            coords = struct[i].coords
            x1 = random.uniform(-1,1)
            x2 = random.uniform(-1,1)
            while (x1**2+x2**2>=1):
                x1=random.uniform(-1,1)
                x2=random.uniform(-1,1)
            coords[0]+=(2*x1*np.sqrt(1-x1**2-x2**2)*disp)
            coords[1]+=(2*x2*np.sqrt(1-x1**2-x2**2)*disp)
            coords[2]+=((1-2*(x1**2+x2**2))*disp)
            struct.replace(i, struct[i].specie, coords, True)
        mg.write_structure(struct,childpath+'/POSCAR')
