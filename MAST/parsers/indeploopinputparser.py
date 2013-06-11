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

import numpy as np
import pymatgen as pmg

from MAST.utility import InputOptions
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import MAST2Structure
from MAST.utility.mastfile import MASTFile

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

class IndepLoopInputParser(MASTObj):

    """Scans an input file for "indeploop" keyword and copies it into
        many input files.
        Attributes:
            self.indeploop <str>: flag for independent looping
            self.loop_delim <str>: character for delimiting loops
            self.loop_start <str>: character indicating start of loop
            self.loop_end <str>: character indicating end of loop
            self.baseinput <MASTFile>: MASTFile created from *.inp input file

        Looping must be indicated at the beginning of the line, and the
        text to be looped must be complete:
        indeploop mast_kpoints (3x3x3 G, 5x5x5 G, 2x2x2 M, 4x4x4 M) 
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.indeploop = "indeploop"
        self.loop_delim = ","
        self.loop_start = "("
        self.loop_end = ")"
        self.baseinput = MASTFile(self.keywords['inputfile'])
    
    def scan_for_indep_loop(self):
        """Scan for independent loops.
            Return:
                loopdict <dict>: Dictionary where loop line indices in the file
                                are keys, with values
                                'looplist', containing a list of strings, and
                                'prepend', containing any text to prepend
                                'append', containing any text to append
                    e.g. loopdict[12]['looplist']=['3x3x3 G','5x5x5 G']
                                     ['prepend']='mast_kpoints '
                                     ['append']=''
        """
        loopdict=dict()
        realline=""
        split1=""
        split2=""
        looplist=""
        numlines = len(self.baseinput.data)
        lidx =0
        while lidx < numlines:
            dataline = self.baseinput.data[lidx].strip()
            if self.indeploop in dataline:
                loopdict[lidx] = dict()
                realline=dataline.split(' ',1)[1] #strip off indeploop
                split1 = realline.split(self.loop_start)
                split2 = split1[1].split(self.loop_end)
                loopdict[lidx]['prepend'] = split1[0]
                loopdict[lidx]['append'] = split2[1]
                loopdict[lidx]['looplist'] = split2[0]
            lidx = lidx + 1
        return loopdict
        
    def prepare_looped_lines(self, loopdict):
        """Prepare looped lines from looping dictionary.
            Args:
                loopdict <dict>: dictionary of looping items from
                                    scan_for_indep_loop
            Returns:
                loopline_dict <dict of dict>: dictionary of different lines for
                                        input files, with keys being
                                        the looped line index
                
        """
        loopline_dict=dict()
        for lidx in loopdict.keys():
            loopline_dict[lidx] = dict()
            lct=0
            for loopitem in loopdict[lidx]['looplist'].split(self.loop_delim):
                newline = loopdict[lidx]['prepend'] + loopitem + loopdict[lidx]['append'] + '\n'
                loopline_dict[lidx][str(lidx) + '-' + str(lct)] = newline
                lct = lct + 1
        return loopline_dict

    def prepare_looped_datasets(self, loopline_dict):
        """Prepare looped datasets from looping lines.
            Args:
                loopline_dict <dict of list>: lines to loop over
            Returns:
                datasets_dict <dict of list>: full datasets for new input files
        """
        datasets_dict=dict()
        flat_dict=dict()
        loopidx_list = loopline_dict.keys()
        prod_list=list()
        for loopidx in loopline_dict.keys():
            prod_list.append(loopline_dict[loopidx].keys())
            for loopvalidx in loopline_dict[loopidx].keys():
                flat_dict[loopvalidx] = loopline_dict[loopidx][loopvalidx]
        import itertools
        loopcombs = itertools.product(*prod_list) #asterisk flattens the list
        stopiter = 0
        while not stopiter:
            try:
                mycomb = loopcombs.next()
            except StopIteration:
                stopiter = 1
            if stopiter == 0:
                newdata = list(self.baseinput.data)
                for lvalidx in list(mycomb):
                    lidx = int(lvalidx.split('-')[0])
                    newdata[lidx] = flat_dict[lvalidx]
                datasets_dict[mycomb] = newdata
        return datasets_dict

    def create_input_files(self, datasets_dict):
        """Create independently looped input files.
            Args:
                datasets_dict <dict of list>

            Returns:
                createdfiles <list of str>: list of newly-created file names
            Creates an input file for each entry.


        """
        ifname = self.keywords['inputfile']
        dirstem = os.path.dirname(ifname)
        basename = os.path.basename(ifname).split('.')[0]
        createdfiles=list()
        if dirstem == "":
            dirstem = os.getcwd()
        for didx in datasets_dict.keys():
            newfile = MASTFile()
            newfile.data = list(datasets_dict[didx])
            newname = newfile.to_unique_file(dirstem, 
                        basename + '_loop', '.inp', 10000)
            createdfiles.append(newname)
        return createdfiles

