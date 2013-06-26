############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
import os
import fnmatch

import numpy as np
import pymatgen as pmg

from MAST.utility import InputOptions
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import MAST2Structure
from MAST.utility import PickleManager

ALLOWED_KEYS = {\
                 'inputfile'    : (str, 'mast.inp', 'Input file name'),\
               }

MAST_KEYWORDS = {'program': 'vasp',
                 'system_name': 'mast',
                 'scratch_directory': os.path.expanduser(os.environ['MAST_SCRATCH']),
                }

STRUCTURE_KEYWORDS = {'posfile': None,
                      'spacegroup': None,
                      'symmetry_only': False,
                      'coord_type': 'cartesian',
                      'atom_list': None,
                      'coordinates': None,
                      'lattice': None,
                      'primitive': False,
                      'structure': None
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

class InputParser(MASTObj):
    """Class for parsing a *.inp format input file.
        Attributes:
            self.section_end <str>: delimiter for the end of a section
            self.delimiter <str>: character(s) for breaking up each line
            self.section_parsers <dict>: supported sections
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.section_end = '$end'
        self.delimiter = ' ' # how we're breaking up each line
        self.section_parsers = {\
                'mast'     : self.parse_mast_section,
                'structure' : self.parse_structure_section,
                'ingredients' : self.parse_ingredients_section,
                'defects'  : self.parse_defects_section,
                'recipe'   : self.parse_recipe_section,
                'neb'      : self.parse_neb_section,
                'chemical_potentials' : self.parse_chemical_potentials_section,
                               }

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

        eof = False
        while (not eof):
            line = infile.readline().lower()
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
                print '\nFound section %s.  Reading in options.' % section_name
            elif (self.section_end in line) and (section_name):
                self.section_parsers[section_name](section_name, 
                        section_content, options)
                print 'Finished parsing the %s section.' % section_name
            else:
                line = line.strip()
                if (line):
                    section_content.append(line)
        infile.close()

        return options

    def parse_mast_section(self, section_name, section_content, options):
        """Parses the mast section and populate the options."""
        mast_dict = MAST_KEYWORDS.copy()

        for line in section_content:
            line = line.split(self.delimiter)
            if (line[0] not in mast_dict):
                error = 'Section keyword %s not recognized' % line[0]
                MASTError(self.__class__.__name__, error)
                return
            else:
                mast_dict[line[0]] = line[1]

        for key, value in mast_dict.items():
            options.set_item(section_name, key, value)
        
        if options.get_item(section_name,'program') == 'vasp':
            if os.getenv('VASP_PSP_DIR') == None:
                raise MASTError(self.__class__.__name__, "Input file specifies program vasp, but no POTCAR directory is set in environment variable VASP_PSP_DIR")

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
        for line in section_content:
            line = line.split(self.delimiter)

            if (line[0] in structure_dict):
                structure_dict[line[0]] = line[1]
            elif ('begin' in line[0]):
                subsection = line[1]
                subsection_list = list()
            elif ('end' not in line):
                subsection_list.append(line)
            elif ('end' in line):
                subsection_dict[subsection] = subsection_list

        # GRJ: Since we lowercase the whole input file, and filenames may not
        # conform to this, we examine the current directory for any files that
        # may match and use that.  If there are multiple matches, we'll throw
        # an error.
        if (structure_dict['posfile'] is not None): # Do we have a geometry file?
            # First build a list of likely files
            file_list = [file for file in os.listdir('.') if fnmatch.fnmatch(file.lower(), structure_dict['posfile'])]
            if (len(file_list) > 1):
                # If we have multiple files with the same name, but different capitalization, throw an error here
                print 'Found mutliple files with the name %s' % structure_dict['posfile']
                print 'Found the files:'
                for file in file_list:
                    print file
                error = 'Found ambiguous file names'
                MASTError(self.__class__.__name__, error)
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
                lattice = np.array(value, dtype='float')
                structure_dict['lattice'] = lattice
            if (key == 'elementmap'):
                for elline in value:
                    elkey = elline[0].strip().upper() #all caps
                    elname = elline[1].strip().title() #Title case
                    element_map[elkey] = elname

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

        for line in section_content:
            line = line.split(self.delimiter)

            if (line[0] == 'coord_type'):
                coord_type = line[1]
            elif (line[0] in defect_list) and (not multidefect):
                type_dict = dict()
                label = None

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
                    label = str(count)

                defect = {'charge': charge,
                          'subdefect_1': type_dict}

                defect_types['defect_%s' % label] = defect

                count += 1
            elif (line[0] == 'begin'):
                defect = dict()
                multidefect = True
                subcount = 1
                defect['charge'] = charge

                try:
                    label = line[1]
                except IndexError:
                    label = None
            elif (line[0] == 'end'):
                if (not label):
                    label = str(count)

                defect_types['defect_%s' % label] = defect
                count += 1
                multidefect = False
            elif ('charge' in line[0]):
                charge_range = line[0].split('=')[-1].split(',')
                defect['charge'] = range(int(charge_range[0]),
                                         int(charge_range[1])+1)
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

        for line in section_content:
            line = line.split(self.delimiter)
            if (line[0] not in recipe_dict):
                error = 'Section keyword %s not recognized' % line[0]
                MASTError(self.__class__.__name__, error)
                return
            elif (line[0] == 'recipe_file'):
                try:
                    recipe_path = os.environ['MAST_RECIPE_PATH']
                except KeyError:
                    error = 'MAST_RECIPE_PATH environment variable not set'
                    MASTError(self.__class__.__name__, error) 
                recipe_dict['recipe_file'] = '%s/%s' % (recipe_path, line[1])
            else:
                recipe_dict[line[0]] = line[1]
        for key, value in recipe_dict.items():
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
            elif ('end' not in line):
                opt = line.split()
                # print opt
                if (opt[0] == 'mast_kpoints'):
                    kpts = map(int, opt[1].split('x'))
                    # Second option after mast_kpoints tells where to 
                    # center the k-point mesh.
                    # If it's there, we append it onto the k-point grid list.
                    if len(opt) > 2:
                        kpts.append(opt[2])
                    ingredient_dict[opt[0]] = kpts
                else:
                    ingredient_dict[opt[0]] = opt[1]
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
        neblines = dict()
        images = 0

        neblabel=""
        for line in section_content:
            type_dict = dict()
            line = line.strip()
            if 'images' in line:
                line = line.split(self.delimiter)
                images = int(line[1])
            elif 'begin' in line:
                line = line.split(self.delimiter)
                neblabel = line[1].strip()
                neblines[neblabel]=list()
            elif 'end' in line:
                neblabel = ""
            else:
                line = line.split(',') #use commas for elementl ines
                neblines[neblabel].append(line)

        options.set_item(section_name, 'images', images)
        options.set_item(section_name, 'neblines', neblines)

    def parse_chemical_potentials_section(self, section_name, section_content,
                                          options):
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


