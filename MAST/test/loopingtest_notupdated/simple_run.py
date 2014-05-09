from MAST.recipe.recipeinput import RecipeInput
from MAST.buffet.buffetmanager import Buffet

####################################
# input parameters
####################################

#mast info
PROGRAM     = "VASP"
SYSTEM_NAME = "MgAl"

#structure info
COORD_TYPE  = "fractional"
COORDINATES = {\
                 'Mg' : '0.000000 0.000000 0.000000',\
                 'Al' : '0.500000 0.500000 0.500000',\
              }
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
#simple run
recipe_input = RecipeInput(recipe_name="simple_recipe")

#add input parameters
recipe_input.addMastInfo(PROGRAM, SYSTEM_NAME)
recipe_input.addStructureInfoFromCoords(COORD_TYPE, COORDINATES, LATTICE)
recipe_input.addGlobalIngredientsInfo(ING_GLOBAL_PARAM)
recipe_input.addIngredientsInfo(OPTIMIZE_ING, OPTIMIZE_INFO)
recipe_input.addRecipeInfo(RECIPE_FILE)

#add recipe input to buffet and cook it
buffet_obj = Buffet(name="simple_recipe")
buffet_obj.addRecipe(recipe_input)
buffet_obj.cook()


