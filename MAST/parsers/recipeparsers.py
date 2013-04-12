############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
from MAST.utility import MASTObj
from MAST.utility import InputOptions

ALLOWED_KEYS = {\
                 'templateFile'    : (str, 'recipe.template', 'template file name'),\
                 'inputOptions'    : (InputOptions, InputOptions(), 'input options parsed using input parser'),\
                 'personalRecipe'  : (str, 'sic.recipe', 'personalized recipe file'),\
               }

class RecipeParser(MASTObj):
    """Parses the input file and produces the personalized recipe file
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.input_options   = self.keywords['inputOptions']
        self.ingredient_list = list()

    def parse(self):
        ''' Parses the template recipe file and creates
            the personalized recipe file
        '''
        #fetch required paramaters
        f_ptr           = open(self.keywords['templateFile'], "r")
        o_ptr           = open(self.keywords['personalRecipe'], "w")
        system_name     = self.input_options.get_item("mast", "system_name", "sys_")
        n_param         = self.input_options.get_item("defects", "num_defects", 0)
        n_images        = self.input_options.get_item("neb", "images", 0)
        n_hops_dict     = self.input_options.get_item("neb", "hopfrom_dict", {})

        print system_name, self.input_options.get_item('mast', 'system_name')

        for line in f_ptr.readlines():
            line             = line.strip()
            line             = line.lower()
            processing_lines = []
            #validate the input line
            if not line or line.startswith('#'):
                continue

            #collect ingredents
            line2 = line.split()
            if (line2[0].lower() == 'ingredient'):
                self.ingredient_list.append(line2[2])


            #replace line params<N>
            #step 1 : replace <sys> with system_name
            #step 2 : replace <n-n> with appropriate hop combinations
            #step 3 : replace <img-n> with appropriate image numbers
            #step 4 : replace <n> with appropriate defect numbers
            processing_lines.append(line)
            #step 1
            processing_lines = self.process_system_name(processing_lines, system_name)
            #step 2
            processing_lines = self.process_hop_combinations(processing_lines, n_hops_dict)
            #step 3
            processing_lines = self.process_images(processing_lines, n_images)
            #step 4
            processing_lines = self.process_defects(processing_lines, n_param)

            #dump the processed lines to file
            output_str = "\n".join(processing_lines)
            o_ptr.write("%s\n" % output_str)

        f_ptr.close()
        o_ptr.close()
#        print 'in RecipeParser.parse():', list(set(self.ingredient_list))

    def process_system_name(self, processing_lines, system_name):
        '''replace <sys> with the system name from the input options
        '''
        for index in xrange(len(processing_lines)):
            processing_lines[index] = processing_lines[index].replace('<sys>', system_name)
        return processing_lines

    def process_hop_combinations(self, processing_lines, n_hops_dict):
        '''replace <n-n> with valid combinations mentioned in the
           hopfrom_dict  of the input options
        '''
        new_lines = []
        if not n_hops_dict:
            return processing_lines
        for line in processing_lines:
            if "<n-n>" in line:
                for hop_start, hop_end_list in n_hops_dict.iteritems():
                    for hop_end in hop_end_list:
                        r_str = "%d-%d" % (hop_start, hop_end)
                        new_lines.append(line.replace('<n-n>', r_str))
            else:
                new_lines.append(line)
        return new_lines

    def process_images(self, processing_lines, n_images):
        '''replace <img-n> with the equivalent number of 
           lines based on the number of images found in the
           input_options
        '''
        new_lines = []
        if not n_images:
            return processing_lines
        for line in processing_lines:
            if '<img-n>' in line:
                for index in xrange(n_images):
                    new_lines.append(line.replace('<img-n>', str(index+1)))
            else:
                 new_lines.append(line)
        return new_lines

    def process_defects(self, processing_lines, n_defects):
        '''replace <N> with the equivalent number of lines
           based on the number of defects given in the
           input options
        '''
        new_lines = []
        if not n_defects:
            return processing_lines
        for line in processing_lines:
            if '<n>' in line:
                for index in xrange(n_defects):
                    new_lines.append(line.replace("<n>", str(index+1)))
            else:
                new_lines.append(line)
        return new_lines


    def get_unique_ingredients(self):
        '''fetches the ingredients names'''
        return list(set(self.ingredient_list))
