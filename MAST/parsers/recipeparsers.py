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

    def parse(self):
        ''' Parses the template recipe file and creates
            the personalized recipe file
        '''
        f_ptr           = open(self.keywords['templateFile'], "r")
        o_ptr           = open(self.keywords['personalRecipe'], "w")
        system_name     = self.input_options.get_item("mast", "system_name", "sys_")
        n_param         = self.input_options.get_item("defects", "num_defects", 0)

        for line in f_ptr.readlines():
            line = line.strip()
            #validate the input line
            if not line or line.startswith('#'):
                continue

            #replace line params
            line = line.lower()
            line = line.replace("<sys>", system_name)
            if "<n>" in line:
                for index in xrange(n_param):
                    new_line = line.replace("<n>", str(index))
                    o_ptr.write("%s\n" % new_line)
            else:
                o_ptr.write("%s\n" % line)

        f_ptr.close()
        o_ptr.close()


def main(template_file, personal_recipe_file):
    input_options = InputOptions()
    input_options.set_item("mast", "system_name", "SiC")
    input_options.set_item("defects", "num_defects", 2)
    rp_obj        = RecipeParser(templateFile=template_file, inputOptions=input_options, personalRecipe=personal_recipe_file)
    rp_obj.parse()


if __name__ == "__main__":
    template_file          = "test/recipetest/recipe_template.txt"
    personal_recipe_file   = "test/recipetest/sic_recipe.txt"
    main(template_file, personal_recipe_file)
