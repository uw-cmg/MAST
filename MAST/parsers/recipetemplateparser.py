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

from MAST.utility import MASTObj
from MAST.utility import InputOptions
from MAST.utility import MASTError

ALLOWED_KEYS = {\
                 'templateFile'    : (str, None, 'template file name'),\
                 'inputOptions'    : (InputOptions, None, 'input options parsed using input parser'),\
                 'personalRecipe'  : (str, None, 'personalized recipe file'),\
               }

class RecipeTemplateParser(MASTObj):
    """Class for parsing the template recipe file in the input file and
        creating a "personalized" recipe file where any <sys> and <N> from
        the recipe template have been filled in according to information
        from the input file.
        Attributes:
            self.input_options <InputOptions>: Input options from input file.
            self.template_file <str>: Path to the recipe template file
            self.personal_recipe <str>: Path to the personalized recipe.
            self.ingredient_list <list of str>: List of ingredients mentioned
                                                in the recipe template file.
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.input_options   = self.keywords['inputOptions']
        self.template_file   = self.keywords['templateFile']
        self.personal_recipe = self.keywords['personalRecipe']
        self.ingredient_list = list()

    def parse(self):
        """ Parses the template recipe file and creates
            the personalized recipe file
        """
        if self.template_file is None:
            raise MASTError(self.__class__.__name__, "Template file not provided!")

        if not os.path.exists(self.template_file):
            raise MASTError(self.__class__.__name__, "Template file not found!")

        if self.input_options is None:
            raise MASTError(self.__class__.__name__, "Input Options not provided!")

        if self.personal_recipe is None:
            raise MASTError(self.__class__.__name__, "Personal recipe file not provided!")

        #fetch required paramaters
        f_ptr           = open(self.template_file, "r")
        o_ptr           = open(self.personal_recipe, "w")
        system_name     = self.input_options.get_item("mast", "system_name", "sys")
        n_defects       = self.input_options.get_item("defects", "num_defects", 0)
        d_defects       = self.input_options.get_item("defects","defects")
        n_images        = self.input_options.get_item("neb", "images", 0)
        d_neblines      = self.input_options.get_item("neb", "neblines", {})
        recipe_name     = None

#        print system_name, self.input_options.get_item('mast', 'system_name')

        for line in f_ptr.readlines():
            line             = line.strip()
            line             = line.lower()
            processing_lines = []
            #validate the input line
            if not line or line.startswith('#'):
                continue

            #collect recipe name
            line2 = line.split()
            if (line2[0].lower() == 'recipe'):
                recipe_name = line2[1] 

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
            #step 4
            processing_lines = self.process_defects(processing_lines, n_defects, d_defects)
            #step 2
            processing_lines = self.process_hop_combinations(processing_lines, d_neblines)
            #step 3
            processing_lines = self.process_images(processing_lines, n_images)

            #dump the processed lines to file
            output_str = "\n".join(processing_lines)
            o_ptr.write("%s\n" % output_str)

        f_ptr.close()
        o_ptr.close()
#        print 'in RecipeParser.parse():', list(set(self.ingredient_list))
        return recipe_name

    def process_system_name(self, processing_lines, system_name):
        """replace <sys> with the system name from the input options
        """
        for index in xrange(len(processing_lines)):
            processing_lines[index] = processing_lines[index].replace('<sys>', system_name)
        return processing_lines

    def process_hop_combinations(self, processing_lines, d_neblines):
        """replace <n-n> with neb labels which are keys of the
           neblines dict of the input options.

            Args:
                processing_lines <list of str>: recipe lines to process.
                d_neblines <dict of str>: dictionary of NEB lines.
        """
        new_lines = []
        eval_lines = []
        if not d_neblines:
            return processing_lines
        line=""
        for line in processing_lines:
            if "<n-n>" in line:
                for neblabel in d_neblines.keys():
                    n_line = line.replace('<n-n>', neblabel)
                    eval_lines.append(n_line)
            else:
                eval_lines.append(line)
        line=""
        for line in eval_lines:
            #print "TTM DEBUG line: ", line
            if not 'neb' in line.split('_'):
                #print "TTM DEBUG line safe"
                new_lines.append(line)
            else:
                evalsplit = line.split('child')
                if len(evalsplit) == 1:
                    new_lines.append(line)
                else:
                    childsplit=evalsplit[1].split('_')
                    parsplit=evalsplit[0].split('_')
                    okay=1
                    for neblabel in d_neblines.keys():
                        if okay == 1 and (neblabel in childsplit):
                            parlabels = neblabel.split('-')
                            #print "TTM DEBUG: parlabels: ",parlabels
                            #print "TTM DEBUG: parsplit: ",parsplit
                            if not parlabels[0] in parsplit and not (parlabels[1] in parsplit):
                                if not neblabel in parsplit: #image static
                                    okay=0
                    if okay == 0:
                        pass
                    else:
                        new_lines.append(line)
        return new_lines

    def process_images(self, processing_lines, n_images):
        """replace <img-n> with the equivalent number of 
           lines based on the number of images found in the
           input_options
        """
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

    def process_defects(self, processing_lines, n_defects, d_defects):
        """replace <N> with the equivalent number of lines
           based on the number of defects given in the
           input options

           Args:
            processing_lines <list of str>: recipe lines to proceess
            n_defects <int>: number of defected systems
            d_defects <dict>: dictionary of defects, including labels and 
                                positions.
        """
        import inspect
        print 'GRJ DEBUG: %s.%s' % (self.__class__.__name__, inspect.stack()[0][3])

        new_lines = list()

        if not n_defects:
            return processing_lines

        for line in processing_lines:
            if '<n>' in line:
                for defect_key in d_defects.keys():
                    defect_label = defect_key.split('_')[1] #defect_1, etc.
                    new_lines.append(line.replace("<n>", defect_label))
            else:
                new_lines.append(line)
        return new_lines

    def get_unique_ingredients(self):
        """fetches the ingredients names"""
        return list(set(self.ingredient_list))
