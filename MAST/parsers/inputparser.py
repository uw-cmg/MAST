##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-06-13 by Zhewen Song
##############################################################
import os, re
import time
import fnmatch
import logging

import numpy as np
import pymatgen as pmg
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifParser

from MAST.utility import InputOptions
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import MAST2Structure
from MAST.utility import dirutil
from MAST.utility import Metadata
from MAST.utility import loggerutils
ALLOWED_KEYS = {\
                 'inputfile'    : (str, 'mast.inp', 'Input file name'),\
               }
MAST_KEYWORDS = {
                 'system_name': 'mast',
                 'mast_auto_correct': 'True'
                }

STRUCTURE_KEYWORDS = {'posfile': None,
                      'spacegroup': None,
                      'symmetry_only': False,
                      'coord_type': 'cartesian',
                      'atom_list': None,
                      'coordinates': None,
                      'lattice': None,
                      'primitive': False,
                      'structure': None,
                      'use_structure_index': "False"
                     }

DEFECTS_KEYWORDS = {'coord_type': 'cartesian',
                    'vacancy': list(),
                    'interstial': list(),
                    'antisite': list(),
                    'substitution': list(),
                   }

INGREDIENTS_KEYWORDS = ['singlepoint',
                        'optimization',
                        'neb',
                   ]

RECIPE_KEYWORDS = {'recipe_file': None,
                  }

PERS_RECP_KEYWORDS = {'personal_recipe_list': None,
                      }

class InputParser(MASTObj):
    """Class for parsing a *.inp format input file.
        Attributes:
            self.section_end <str>: delimiter for the end of a section
            self.delimiter <str>: character(s) for breaking up each line
            self.section_parsers <dict>: supported sections
            self.logger <logging Logger>: logger, either recipe-level or MAST input level
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.section_end = '$end'
        self.delimiter = ' ' # how we're breaking up each line
        self.section_parsers = {\
                'mast'     : self.parse_mast_section,
                'structure' : self.parse_structure_section,
                'scaling' : self.parse_scaling_section,
                'ingredients' : self.parse_ingredients_section,
                'defects'  : self.parse_defects_section,
                'recipe'   : self.parse_recipe_section,
                'neb'      : self.parse_neb_section,
                'chemical_potentials' : self.parse_chemical_potentials_section,
                'summary' : self.parse_summary_section,
                'personal_recipe' : self.parse_personal_recipe_section
                               }
        scratchpath = dirutil.get_mast_scratch_path()
        inputlocation = os.path.dirname(self.keywords['inputfile'])
        if (inputlocation == ""):
            self.logger = loggerutils.get_mast_logger("mast input parser")
        elif scratchpath not in inputlocation:
            self.logger = loggerutils.get_mast_logger("mast input parser")
        else:
            self.logger = loggerutils.get_mast_logger("mast input parser %s" % inputlocation)



    def parse(self):
        """Parses information from the *.inp style input file.
            Returns:
                options <InputOptions>: a dictionary of sections and values,
                                        where each section is its own
                                        dictionary:
                                        options.options['section'] = dict()
        """
        options   = InputOptions()
        infile    = file(self.keywords['inputfile'])
        section_name=""
        eof = False
        while (not eof):
            line = infile.readline() #.lower() #Remove the lowercasing
            # This works because the blank lines in the file are really 
            # newline (\n)!
            # This is why strip() didn't work, cause it removed that \n, 
            # creating an empty string.
            if (not line):
                eof = True
            elif (line.startswith('#')) or (line.startswith('!')):
            # Lines that start with a # or a ! are treated as a comment and 
            # are skipped
                pass
            elif (line.startswith('$')) and (self.section_end not in line):
                section_name = line[1:].strip()

                if section_name not in self.section_parsers:
                    error = 'Section %s not recognized' % section_name
                    MASTError(self.__class__.__name__, error)
                    return

                section_content = list()
                self.logger.debug('Found section %s.  Reading in options.' % section_name)
            elif (self.section_end in line) and (section_name):
                self.section_parsers[section_name](section_name, 
                        section_content, options)
                self.logger.debug('Finished parsing the %s section.' % section_name)
            else:
                if (section_name != 'recipe') and (section_name != 'personal_recipe'): 
                        line = line.strip()
                if (line):
                    section_content.append(line)
        infile.close()

        self.perform_element_mapping(options)
        self.set_structure_from_inputs(options)
        self.logger.debug(options) # change to debug level
        self.validate_execs(options)
        return options

    def parse_mast_section(self, section_name, section_content, options):
        """Parses the mast section and populate the options."""
        mast_dict = MAST_KEYWORDS.copy()

        for line in section_content:
            line = line.split()
            if (line[0] not in mast_dict):
                error = 'Section keyword %s not recognized' % line[0]
                MASTError(self.__class__.__name__, error)
                return
            else:
                mast_dict[line[0]] = line[1]

        for key, value in mast_dict.items():
            options.set_item(section_name, key, value)
        
    def parse_scaling_section(self, section_name, section_content, options):
        """Parses the scaling section and populate the options."""
        scaling_dict = dict()
        for line in section_content:
            if (line.startswith('begin')):
                size_name = ' '.join(line.split()[1:])
                size_dict = dict()
            elif ('end' not in line):
                opt = line.split()
                if (opt[0] == 'mast_size'):
                    matrix = line.split('[')[1].split(']')[0]
                    size_dict['mast_size'] = matrix
                elif (opt[0] == 'mast_kpoints'):
                    size_dict['mast_kpoints'] = ' '.join(opt[1:])
            elif ('end' in line):
                scaling_dict[size_name] = size_dict
        for key,value in scaling_dict.items():
            options.set_item(section_name, key, value)


    def parse_structure_section(self, section_name, section_content, options):
        """Parses the structure section and populate the options.
            Does not create the structure.

            Format is along the lines of:
                coord_type fractional

                begin coordinates
                Ga 0.000000 0.000000 0.000000
                Ga 0.500000 0.500000 0.000000
                Ga 0.000000 0.500000 0.500000
                Ga 0.500000 0.000000 0.500000
                As 0.250000 0.250000 0.250000
                As 0.750000 0.750000 0.250000
                As 0.250000 0.750000 0.750000
                As 0.750000 0.250000 0.750000
                end

                begin lattice
                2.0 0.0 0.0
                0.0 2.0 0.0
                0.0 0.0 2.0
                end

            Note that coord_type will default to "cartesian" if not specified.

        """
        # Initialize with default values
        structure_dict = STRUCTURE_KEYWORDS.copy() 

        subsection_dict = dict()
        for myline in section_content:
            line = myline.split()

            if (line[0] in structure_dict):
                structure_dict[line[0]] = line[1]
            elif ('begin' in line[0]):
                subsection = line[1]
                subsection_list = list()
            elif ('end' not in line):
                lsplit = myline.split()
                lineval = list()
                for lval in lsplit:
                    lval.strip()
                    if len(lval) > 0:
                        lineval.append(lval)
                subsection_list.append(lineval)
            elif ('end' in line):
                subsection_dict[subsection] = subsection_list

        # GRJ: Since we lowercase the whole input file, and filenames may not
        # conform to this, we examine the current directory for any files that
        # may match and use that.  If there are multiple matches, we'll throw
        # an error.
        if (structure_dict['posfile'] is not None): # Do we have a geometry file?
            # TTM for issue #423
            posfiletry = os.path.join(os.getcwd(), structure_dict['posfile'])
            if os.path.isfile(posfiletry):
                self.logger.info("Using posfile %s in directory %s" % (structure_dict['posfile'], os.getcwd()))
            else:
                # First build a list of likely files
                origindir = os.path.dirname(self.keywords['inputfile'])
                if origindir == "":
                    origindir = os.getcwd()
                metatry = os.path.join(os.getcwd(), 'metadata.txt')
                if os.path.isfile(metatry):
                    myrecipemeta = Metadata(metafile=metatry)
                    origindir = myrecipemeta.search_data('origin_dir')[1]
                myfiles = dirutil.walkfiles(origindir)
                file_list=list()
                for myfile in myfiles:
                    if structure_dict['posfile'] in myfile:
                        file_list.append(myfile)
                if (len(file_list) > 1):
                    # If we have multiple files with the same name, but different capitalization, throw an error here
                    self.logger.warning('Found multiple files with the name %s' % structure_dict['posfile'])
                    self.logger.warning('Found the files:')
                    for file in file_list:
                        self.logger.warning(file)
                    error='Found ambiguous file names'
                    raise MASTError(self.__class__.__name__, error)
                elif len(file_list) == 0:
                    raise MASTError(self.__class__.__name__, "No structure file %s found in %s" % (structure_dict['posfile'], origindir))
                else:
                    structure_dict['posfile'] = file_list[0]
        # print 'in InputParser.parse_structure_section:', subsection_dict
        # TM
        element_map = dict()
        atom_list = list()
        for key, value in subsection_dict.items():
            if (key == 'coordinates'):
                value = np.array(value)
                # Here we use .title() to re-capitalize the first letter of all 
                # the atomic symbols to comply with what
                # pymatgen needs
                atom_list = [val.title() for val in value[:, 0]]
                structure_dict['atom_list'] = atom_list
                coordinates = np.array(value[:, 1:], dtype='float')
                structure_dict['coordinates'] = coordinates
            if (key == 'lattice'):
                #print "VALUE: ", value
                lattice = np.array(value, dtype='float')
                structure_dict['lattice'] = lattice
            if (key == 'elementmap'):
                for elline in value:
                    elkey = elline[0].strip().upper() #all caps
                    elname = elline[1].strip().title() #Title case
                    element_map[elkey]=elname
                structure_dict['element_map'] = element_map
        if len(element_map) > 0 and len(atom_list) > 0:
            new_atom_list = list()
            for atomval in atom_list:
                if atomval.upper() in element_map.keys():
                    new_atom_list.append(element_map[atomval])
                else:
                    new_atom_list.append(atomval)
            structure_dict['atom_list'] = new_atom_list

        for key, value in structure_dict.items():
            options.set_item(section_name, key, value)
    def parse_defects_section(self, section_name, section_content, options):
        """Parses the defects section and populates the options.
        """
        defect_list = ['antisite', 'vacancy', 'substitution', 'interstitial']
        defect_types = dict()
        multidefect = False
        charge = [0]
        count = 1
        coord_type = 'cartesian'
        threshold = 1.e-4

        for line in section_content:
            #line = line.split(self.delimiter)
            line = line.split()
            if (line[0] == 'coord_type'):
                coord_type = line[1]
            elif (line[0] == 'threshold'):
                threshold = float(line[1])
            elif (line[0] in defect_list) and (not multidefect):
                type_dict = dict()
                label = None
                charge = [0]

                if (len(line) < 5):
                    error = 'Defect specification requires at least 5 arguments.'
                    MASTError(self.__class__.__name__, error)

                # Check for static options
                type_dict['type'] = line[0]
                coord = line[1:4]
                type_dict['coordinates'] = np.array(coord, dtype='float')
                type_dict['symbol'] = line[4]

                if (len(line) > 5):
                    for lin in line[5:]:
                        lin = lin.split('=')
                        if (lin[0] == 'charge'):
                            charge_range = lin[1].split(',')
                            charge = range(int(charge_range[0]),
                                           int(charge_range[1])+1)
                        elif (lin[0] == 'label'):
                            label = lin[1]

                if (not label):
                    label = 'defect_%s' % str(count)

                defect = {'charge': charge,
                          'subdefect_1': type_dict,
                          'coord_type': coord_type,
                          'threshold': threshold,
                          'phonon': dict()}

                #defect_types['defect_%s' % label] = defect
                defect_types[label] = defect

                count += 1
            elif (line[0] == 'begin'):
                defect = dict()
                multidefect = True
                subcount = 1
                defect['charge'] = charge
                defect['coord_type'] = coord_type
                defect['threshold'] = threshold
                defect['phonon'] = dict()

                try:
                    label = line[1]
                except IndexError:
                    label = None
            elif (line[0] == 'end'):
                if (not label):
                    label = 'defect_%s' % str(count)

                defect_types[label] = defect
                count += 1
                multidefect = False
            elif ('charge' in line[0]):
                charge_range = line[0].split('=')[-1].split(',')
                defect['charge'] = range(int(charge_range[0]),
                                         int(charge_range[1])+1)
            elif ('phonon' in line[0]):
                plabel = line[1]
                p_center = ' '.join(line[2:5])
                p_radius = float(line[5])
                p_thresh = 0.1
                if len(line) > 6:
                    p_thresh = float(line[6])
                defect['phonon'][plabel]=dict()
                defect['phonon'][plabel]['phonon_center_site'] = p_center
                defect['phonon'][plabel]['phonon_center_radius'] = p_radius
                defect['phonon'][plabel]['threshold'] = p_thresh
            else:
                type_dict = dict()

                if (len(line) < 5):
                    error = 'Defect specification requires at least 5 arguments.'
                    MASTError(self.__class__.__name__, error)

                # Check for static options
                type_dict['type'] = line[0]
                coord = line[1:4]
                type_dict['coordinates'] = np.array(coord, dtype='float')
                type_dict['symbol'] = line[4]

                defect['subdefect_%i' % subcount] = type_dict
                # print 'Rawr!', defect
                subcount += 1
        options.set_item(section_name, 'num_defects', count-1)
        options.set_item(section_name, 'defects', defect_types)
        options.set_item(section_name, 'coord_type', coord_type)

    def parse_recipe_section(self, section_name, section_content, options):
        """Parses the recipe section and populates the options."""
        recipe_dict = RECIPE_KEYWORDS.copy()

        if not section_content:
            error = 'Recipe section is not specified'
            MASTError(self.__class__.__name__, error)
            return	
        else:
            recipe_dict['recipe_file'] = section_content

        for key, value in recipe_dict.items():
            options.set_item(section_name, key, value)

    def parse_personal_recipe_section(self, section_name, section_content, options):
        """Parses the recipe section and populates the options."""
        personal_recipe_dict = PERS_RECP_KEYWORDS.copy()

        if not section_content:
            error = 'Personal Recipe section is not specified or is empty'
            MASTError(self.__class__.__name__, error)
            return
        else:
            personal_recipe_dict['personal_recipe_list'] = section_content

        for key, value in personal_recipe_dict.items():
            options.set_item(section_name, key, value)

    def parse_ingredients_section(self, section_name, section_content, options):
        """Parses the ingredients section and populates the options.
            Section takes the form of:
                $ingredients
                begin ingredients_global
                mast_kpoints 3x3x3
                mast_xc pbe
                end

                begin singlepoint
                encut 400
                end

                begin optimize
                encut 300
                ibrion 2
                end

                $end

            mast_kpoints are parsed out as a 3 index list of integers. 
            Everything else is parsed out as a string.

            Anything in ingredients_global is then appended onto each 
            individual ingredient.
        """

        global_dict = dict()
        ingredients_dict = dict()

        for line in section_content:
            if (line.startswith('begin')):
                # Each ingredient section starts with "begin".
                # Check for this line, and initialize the individual
                # ingredient dictionary
                ingredient_name = line.split()[1]
                ingredient_dict = dict()
            elif (not (line == 'end')):
                opt = line.split()
                # print opt
                if (opt[0] == 'mast_kpoints'):
                    try: 
                        kpts = map(int, opt[1].split('x'))
                        if len(opt) > 2:
                            kpts.append(opt[2])
                    except ValueError: 
                        kpts = opt[1:]
                    # Second option after mast_kpoints tells where to 
                    # center the k-point mesh.
                    # If it's there, we append it onto the k-point grid list.
                    ingredient_dict[opt[0]] = kpts
                elif (opt[0] == 'mast_pp_setup'):
                    psp_dict = dict()
                    for key in opt[1:]:
                        key = key.split('=')
                        ref = key[0].title()
                        val = str().join(key[1][0].upper() + key[1][1:])
                        psp_dict[ref] = val
                    ingredient_dict[opt[0]] = psp_dict
                #elif (opt[0] == 'mast_exec'):
                #    ingredient_dict[opt[0]] = ' '.join(opt[1:]) #preserve whole line
                #elif (opt[0] == 'mast_strain'):
                #    ingredient_dict[opt[0]] = ' '.join(opt[1:]) #preserve whole line
                #elif (opt[0] == 'ptemp'):
                #    ingredient_dict[opt[0]] = ' '.join(opt[1:]) #preserve whole line 
                #elif (opt[0] == 'rwigs'):
                #    ingredient_dict[opt[0]] = ' '.join(opt[1:]) #preserve whole line
                #elif (opt[0] == 'mast_setmagmom'):
                #    ingredient_dict[opt[0]] = ' '.join(opt[1:]) #preserve whole line
                elif (opt[0] == 'mast_coordinates'):
                    shortsplit = opt[1].split(",")
                    filesplit=list()
                    origindir = os.getcwd()
                    metatry = os.path.join(os.getcwd(), 'metadata.txt')
                    if os.path.isfile(metatry):
                        myrecipemeta = Metadata(metafile=metatry)
                        origindir = myrecipemeta.search_data('origin_dir')[1]
                    myfiles = dirutil.walkfiles(origindir)
                    for shortname in shortsplit:
                        for fullfile in myfiles:
                            if shortname.strip() in os.path.basename(fullfile):
                                filesplit.append(fullfile)
                    if not (len(filesplit) == len(shortsplit)):
                        raise MASTError(self.__class__.__name__, "Not all files given by %s were found in %s." % (shortsplit, origindir))
                    ingredient_dict[opt[0]] = filesplit
                else:
                    ingredient_dict[opt[0]] = ' '.join(opt[1:]) #preserve whole line
                    #ingredient_dict[opt[0]] = opt[1] #old behavior took only the first part
                if (opt[0] == 'mast_program') and (opt[1] == 'vasp' or opt[1] == 'vasp_neb'):
                    if os.getenv('VASP_PSP_DIR') == None:
                        raise MASTError(self.__class__.__name__, "Input file specifies program vasp, but no POTCAR directory is set in environment variable VASP_PSP_DIR")
            elif ('end' in line):
                # Each ingredient section ends with "end", if present finish 
                # that current section and assign
                # the neccessary element in the ingredients dictionary and 
                # create the global dictionary
                if (ingredient_name == 'ingredients_global'):
                    global_dict = ingredient_dict
                else:
                    ingredients_dict[ingredient_name] = ingredient_dict

        # Each value in ingredients_dict is a dictionary containing the relevant
        # ingredient and option(s).
        for ing_key, ing_value in ingredients_dict.items():
            options.set_item(section_name, ing_key, ing_value)

        options.set_item(section_name, 'global', global_dict)

    def parse_neb_section(self, section_name, section_content, options):
        """Parses the neb section and populates the options.
            
            The number of images is specified and must be the same
            throughout the recipe.

            Each neb gets a label which should correspond to a defect group
            or defect group number.

            Within each subsection, each line specifies which atom is moving,
            giving its element name and its approximate coordinates from
            the starting unrelaxed structure.

            $neb

            images 3

            begin crowdion1-crowdion2
            Cr, 0 0.3 0, 0 0 0
            Ni, 0 0.6 0, 0 0.3 0
            end

            begin int1-int2
            Mg, 0 0 0.5, 0.3 0.1 2
            end

            $end

            In the crowdion example, the starting undefected structure is
            defected three times in each defect group: 1 vacancy and two
            antisites. The neb section allows us to unambiguously match up 
            the atoms and determine which atoms are moving, and where. 
        """
        nebs = dict()
        images = 0
        neblabel=""
        for line in section_content:
            line = line.strip()
            #if 'images' in line:
            #    line = line.split()
            #    images = int(line[1])
            if 'begin' in line:
                line = line.split()
                neblabel = line[1].strip()
                nebs[neblabel]=dict()
                nebs[neblabel]['lines']=list()
                nebs[neblabel]['phonon']=dict()
            elif 'end' in line:
                neblabel = ""
            else:
                line = line.split(',') #use commas for element lines
                if 'phonon' in line[0]:
                    sline = line[0].split()
                    plabel = sline[1]
                    p_center = ' '.join(sline[2:5])
                    p_radius = float(sline[5])
                    p_thresh = 0.1
                    if len(sline) > 6:
                        p_thresh = float(sline[6])
                    nebs[neblabel]['phonon'][plabel]=dict()
                    nebs[neblabel]['phonon'][plabel]['phonon_center_site'] = p_center
                    nebs[neblabel]['phonon'][plabel]['phonon_center_radius'] = p_radius
                    nebs[neblabel]['phonon'][plabel]['threshold'] = p_thresh
                elif 'images' in line[0]:
                    nebs[neblabel]['images'] = int(line[0].split()[1])
                else:
                    nebs[neblabel]['lines'].append(line)

        #options.set_item(section_name, 'images', images)
        options.set_item(section_name, 'nebs', nebs)

    def parse_chemical_potentials_section(self, section_name, section_content, options):
        """Parses the chemical_potentials section and populates the options.
            Section uses the standard begin...end subsection structure, but with
            a modification: instead of strict subsection titles (i.e. structure,
            lattice etc.), subsection titles are the conditions under which the
            chemical potentials are for.  Any combination of white spacing is
            allowed, however note that all conditions will be converted to lower
            case first!

            Using a GaAs example, hypothetically we could have a 
            chemical_potential section like the following:
                $chemical_potentials
                begin As rich
                Ga 4.5
                As 3.5
                end

                begin Ga rich
                Ga 3.5
                As 4.5
                end
                $end

            For charged defects, note that the chemical potential of the
            electron will be calculated automatically from the Fermi level and
            appropiate potential shift correction.  However, if one desires to
            manually give an electron chemical potential (for whatever reason!)
            this can be done by specifiying \'electron\' in the appropiate
            condition.
        """
        chempot_dict = dict()

        for line in section_content:
            if (line.startswith('begin')):
                condition_name = ' '.join(line.split()[1:])
                condition_dict = dict()
            elif ('end' not in line):
                opt = line.split()
                condition_dict[opt[0]] = float(opt[1])
            elif ('end' in line):
                chempot_dict[condition_name] = condition_dict

        for key, value in chempot_dict.items():
            options.set_item(section_name, key, value)


    def parse_phonon_section(self, section_name, section_content, options):
        """Parses the phonon section and populates the options.
            
            The phonon_center_site is a coordinate (as a string, like "0 0 0")
                for the center of the phonon calculations. Give a coordinate
                of an atom for that atom.
                Do not set this value in order to calculate phonons for all
                atoms.
            The phonon_center_radius is a radius in ANGSTROMS to search for
                neighbors to the phonon_center_site, and also take into account
                their vibrations.
                Do not set this value or set to 0 to only calculate phonons
                for the atom given in phonon_center_site.
            All phonons must be enclosed in a begin and end tag within the 
            phonon section. Phonons should be labeled
            xxxx_phonon_label and the label should correspond to the 
            auto-generation of the
            recipe (e.g. corresponding exactly to defect group names like vac1
            or to neb image names like vac1-vac2_02)

            $phonon

            begin crowdion1
            phonon_center_site 0 0.3 0
            phonon_center_radius 2
            end
            
            begin crowdion1-crowdion2_02
            phonon_center_site 0 0.15 0
            phonon_center_radius 2
            end

            $end
        """
        raise MASTError(self.__class__.__name__,"Phonon section in input files *.inp is obsolete. Phonons now go in the defects and neb sections.")
        phonon = dict()
        label=""
        for line in section_content:
            line = line.strip()
            if 'begin' in line:
                line = line.split()
                label = line[1].strip()
                phonon[label]=dict()
            elif 'end' in line:
                pass
            else:
                line = line.split()
                phonon[label][line[0]] = line[1]

        options.set_item(section_name, 'phonon', phonon)


    def perform_element_mapping(self, input_options):
        """Perform element mapping if there is an element map specified,
            so that defects and NEB lines get the correct element name.
            Args:
                input_options <InputOptions>
            Returns:
                modifies input_options dictionary entries in place.
        """
        #print self.input_options.get_section_keys('structure')
        eldict=dict()
        if 'element_map' in input_options.get_section_keys('structure'):
            eldict = input_options.get_item('structure','element_map')
        if len(eldict) == 0:
            return
        if 'defects' in input_options.get_sections():
            defdict = input_options.get_item('defects','defects')
            for dkey in defdict.keys():
                for sdkey in defdict[dkey].keys():
                    if 'subdefect' in sdkey:
                        symbol = defdict[dkey][sdkey]['symbol'].upper()
                        if symbol in eldict.keys():
                            defdict[dkey][sdkey]['symbol'] = eldict[symbol]
        if 'neb' in input_options.get_sections():
            nebdict = input_options.get_item('neb','nebs')
            for nebkey in nebdict.keys():
                nlinenum = len(nebdict[nebkey]['lines'])
                nlinect=0
                while nlinect < nlinenum:
                    symbol = nebdict[nebkey]['lines'][nlinect][0].upper()
                    if symbol in eldict.keys():
                        nebdict[nebkey]['lines'][nlinect][0] = eldict[symbol]
                    nlinect = nlinect + 1
        return



    def validate_execs(self, input_options):
        """Make sure each ingredient has a mast_exec line.
            Args:
                input_options <InputOptions>
        """
        have_exec = False
        for ingredient, options in input_options.get_item('ingredients').items():
            if 'mast_exec' in options:
                have_exec = True
                break
            else:
                have_exec = False

        if (not have_exec):
            error = 'mast_exec keyword not found in the $ingredients section'
            raise MASTError(self.__class__.__name__, error)
    def set_structure_from_inputs(self, input_options):
        """Make a pymatgen structure and update the
            structure key.
            Args:
                input_options <InputOptions>
        """
        strposfile = input_options.get_item('structure','posfile')
        if strposfile is None:
            iopscoords=input_options.get_item('structure','coordinates')
            iopslatt=input_options.get_item('structure','lattice')
            iopsatoms=input_options.get_item('structure','atom_list')
            iopsctype=input_options.get_item('structure','coord_type')
            structure = MAST2Structure(lattice=iopslatt,
                coordinates=iopscoords, atom_list=iopsatoms,
                coord_type=iopsctype)
        elif ('poscar' in strposfile.lower()):
            structure = Poscar.from_file(strposfile).structure
        elif ('cif' in strposfile.lower()):
            structure = CifParser(strposfile).get_structures()[0]
        else:
            error = 'Cannot build structure from file %s' % strposfile
            raise MASTError(self.__class__.__name__, error)
        input_options.update_item('structure','structure',structure)
        if input_options.get_item('structure','use_structure_index') in ['True','true','T','t']:
            self.do_structure_indexing(input_options)
        return

    def do_structure_indexing(self, input_options):
        """Index the structure into a dictionary.
        """
        if os.path.exists("structure_index_files"):
            return
        from MAST.ingredients.pmgextend.atom_index import AtomIndex
        myindex = AtomIndex(input_options = input_options)
        myindex.set_up_initial_index()
        return

    def parse_summary_section(self, section_name, section_content, options):
        """Parses the summary section and populates the options."""
        mast_dict = dict()
        for line in section_content:
            line = line.split()
            mast_dict.setdefault(line[0], []).append(line[1].strip())

        for key, value in mast_dict.items():
            options.set_item(section_name, key, value)
