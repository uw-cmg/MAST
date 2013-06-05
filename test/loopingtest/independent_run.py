from MAST.recipe.recipeinput import RecipeInput
from MAST.buffet.buffetmanager import Buffet

####################################
# input parameters
####################################

#mast info
PROGRAM     = "VASP"
SYSTEM_NAME = ["MgAl", "Fe"]

#structure info
COORD_TYPE  = "fractional"
COORDINATES = [\
                 {\
                    'Mg' : '0.000000 0.000000 0.000000',\
                    'Al' : '0.500000 0.500000 0.500000',\
                 },\
                 {\
                    'Fe' : '0.300000 0.300000 0.300000'\
                 }\
              ]
LATTICE     = [\
                 "3.0 0.0 0.0",\
                 "0.0 3.0 0.0",\
                 "0.0 0.0 3.0"\
              ]

#ingredients info
ING_GLOBAL_PARAM = {\
                       'mast_kpoints' : '3x3x3',\
                       'xc'           : 'pbe',\
                   }

#optimize ingredient info
OPTIMIZE_ING  = "optimize"
OPTIMIZE_INFO = {\
                    'encut'  : '300',\
                    'ibrion' : '2',\
                }

#recipe info
RECIPE_FILE = 'recipe_test.txt'

####################################
#independent run
for run_idx in xrange(len(SYSTEM_NAME)):
    recipe_input = RecipeInput(recipe_name="independent_recipe%s" % run_idx)

    #add input parameters
    recipe_input.addMastInfo(PROGRAM, SYSTEM_NAME[run_idx])
    recipe_input.addStructureInfoFromCoords(COORD_TYPE, COORDINATES[run_idx], LATTICE)
    recipe_input.addGlobalIngredientsInfo(ING_GLOBAL_PARAM)
    recipe_input.addIngredientsInfo(OPTIMIZE_ING, OPTIMIZE_INFO)
    recipe_input.addRecipeInfo(RECIPE_FILE)

    #add recipe input to buffet and cook it
    buffet_obj = Buffet(name="independent_recipe_%s" % SYSTEM_NAME[run_idx])
    buffet_obj.addRecipe(recipe_input)
    buffet_obj.cook()


