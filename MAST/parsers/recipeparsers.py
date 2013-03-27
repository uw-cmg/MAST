import re

class RecipeParser:
    def __init__(self):
        pass

    def parse(self, template_file, input_options, personal_recipe_file):
        ''' Parses the template recipe file and creates
            the personalized recipe file
        '''
        print 'template_file =', template_file
        f_ptr           = open(template_file, "r")
        o_ptr           = open(personal_recipe_file, "w")
#        system_name     = input_options.get("mast", dict()).get("system_name", "sys_")
#        n_param         = len(input_options.get("defects", dict()).get("num_defects", 0))
        system_name = input_options.get_item('mast', 'system_name')
        n_param = input_options.get_item('defects', 'num_defects')

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
    input_options = {"mast" : {"system_name" : "SiC"}, "defects" : {"num_defects" : 2}}
    rp_obj        = RecipeParser()
    rp_obj.parse(template_file, input_options, personal_recipe_file)


if __name__ == "__main__":
    template_file          = "recipe_template.txt"
    personal_recipe_file   = "sic_recipe.txt"
    main(template_file, personal_recipe_file)
