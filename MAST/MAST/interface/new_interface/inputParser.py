############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
import numpy as np
from MAST.utility import MASTObj
from MAST.utility import MASTError
from inputOptions import InputOptions

ALLOWED_KEYS     = {\
                       'inputfile'    : (str, 'mast.inp', 'Input file name'),\
                   }


class InputParser(MASTObj):
    """Parses input file and returns the options.
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.section_end     = '$end'
        self.section_parsers = {\
                                    'mast'     : self.parse_mast_section,\
                                    'geometry' : self.parse_geometry_section,\
                                    'unitcell' : self.parse_unit_cell_section,\
                                    'program'  : self.parse_program_section,\
                                    'defects'  : self.parse_defects_section,\
                                    'neb'      : self.parse_neb_section,\
                               }

    def parse(self):
        """Parses information from the input file"""
        options   = InputOptions()
        infile    = file(self.keywords['inputfile'])
        contents  = infile.read()
        infile.close()

        sections  = contents.strip().split(self.section_end)

        for section_content in sections:
            section_content = section_content.strip()
            if not section_content:
                continue

            section_name    = section_content[0][1:]
            if section_name not in self.section_parsers:
                error = 'Section %s not recognized' % section_name
                MASTError(self.__class__.__name__, error)
                return
          
            self.section_parsers[section_name](section_name, section_content[1:], options)

    def parse_mast_section(self, section_name, section_content, options):
        """Parse the mast section and populate the options"""
        raise NotimplementedError

    def parse_geometry_section(self, section_name, section_content, options):
        """Parse the geometry section and populate the options

        Every line is expected to have a list of coordinates
        """
        coordinates = list()
        for con in section_contents:
            coordinates.append(con.split())
        coordinates = np.array(coordinates, dtype='float')
        options.set_item(section_name, 'coordinates', coordinates)


    def parse_unitcell_section(self, section_name, section_content, options):
        """Parse the unit cell section and populate the options

        Every line is expected to have a list of cells
        """
        cell = list()
        for con in section_contents:
            cell.append(con.split())
        cell = np.array(cell, dtype='float')
        options.set_item(section_name, 'cell', cell)

    def parse_program_section(self, section_name, section_content, options):
        """Parse the program section and populate the options"""
        raise NotImplementedError

    def parse_defects_section(self, section_name, section_content, options):
        """Parse the defects section and populate the options"""
        raise NotImplementedError

    def parse_neb_section(self, section_name, section_content, options):
        """Parse the neb section and populate the options"""
        raise NotImplementedError
