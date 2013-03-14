############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
from creator import Creator
from MAST.utility.mastobj import MASTObj

ALLOWED_KEYS     = {\
                       'inputfile'    : (str, 'mast.inp', 'Input file name'),\
                       'ouputfile'    : (str, 'mast.out', 'Output file name'),\
                   } 

class interface(MASTObj):
    """User interface to set up a calculation group.

    Each instance of Interface sets up one calculation group.

    Attributes:
        options <InputOptions object>: used to store the options
                                       parsed from input file
    """

    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.options = None
    
    def start(self):
        """Starts the calculation group from interface.

        This function parses the Input file and fetches the options.
        Then passes these options to the creator class instance.
        """
        parser_obj   = InputParser(inputfile=self.keywords['inputfile'])
        self.options = parser_obj.parse()

        creator_obj  = Creator(self.options)
        creator_obj.start() 
         
